using namespace seqan;

struct SHit{
    Pair <uint64_t, uint64_t> occ;
//     uint8_t errors;
//     DnaString read;
    template <typename TOcc>
    SHit(TOcc & inOcc):
        occ(inOcc)
    {}
};

//why cant i sort occ in a custom way?

bool sHit_smaller(const SHit & x, const SHit & y)
{
    if(getSeqNo(x.occ) == getSeqNo(y.occ))
        return getSeqOffset(x.occ) < getSeqOffset(y.occ);
    else
        return getSeqNo(x.occ) < getSeqNo(y.occ);

}

template<int disT>
bool sHit_similar(const SHit & x, const SHit & y)
{
    return(getSeqNo(x.occ) == getSeqNo(y.occ) && getSeqOffset(x.occ) + disT >= getSeqOffset(y.occ) && getSeqOffset(x.occ) <= getSeqOffset(y.occ) + disT);
}

template<typename TOcc>
bool occs_smaller(TOcc & x, TOcc & y)
{
    if(getSeqNo(x) == getSeqNo(y))
        return getSeqOffset(x) < getSeqOffset(y);
    else
        return getSeqNo(x) < getSeqNo(y);
}

template <unsigned errors, typename TIndex, typename TText, typename TContainer, typename TOptions, typename TDistanceTag>
inline void runAlgo4(TIndex & index, TText const & text, TContainer & c, TOptions const & opt, TDistanceTag const &)
{
    typedef typename TContainer::value_type value_type;

    constexpr uint64_t max_val = std::numeric_limits<value_type>::max();
    typedef Iter<TIndex, VSTree<TopDown<> > > TIter;

    auto const & limits = stringSetLimits(text);

    uint64_t const textLength = length(text); // lengthSum() forwards to length() for a single string

    const uint64_t max_i = textLength - opt.k_length + 1;
    const uint64_t step_size = opt.k_length - opt.overlap + 1;
    #pragma omp parallel for schedule(dynamic, std::max(1ul, max_i/(step_size*opt.threads*50000))) num_threads(opt.threads)
    for (uint64_t i = 0; i < max_i; i += step_size)
    {
        uint64_t max_pos = std::min(i + opt.k_length - opt.overlap, textLength - opt.k_length) + 1;

        // overlap is the length of the infix!
        uint64_t leading = 0, trailing = 0;
        for (uint64_t xx = i; c[xx] != 0 && xx < max_pos; ++xx)
            ++leading;
        for (uint64_t xx = max_pos - 1; c[xx] != 0 && xx >= i; --xx) // TODO: i could theoretically be 0 ... overflow because of unsigned value! but c[...] will be zero for i=0 (unless some really weired scheduling happens)
            ++trailing;

        if (trailing != max_pos - i)
        // doesn't work either: trailing != length - overlap + 1 because last interval might be smaller than length - overlap + 1
        // doesn't word: leading != max_pos - i. trailing is computed last and might have found the full range (while another thread writing) to be non-zero while leading didn't find a full range before!
        {
            uint64_t begin_pos = i + leading;
            uint64_t end_pos = max_pos - trailing; // excluding
            uint64_t new_overlap = opt.k_length - (end_pos - begin_pos) + 1;

            auto scheme = OptimalSearchSchemes<0, errors>::VALUE; // TODO: move out as array
            _optimalSearchSchemeComputeFixedBlocklength(scheme, new_overlap); // only do when new_overlap != overlap

            TIter it_zero_errors[end_pos - begin_pos];
//             value_type hits[end_pos - begin_pos] = {};
            std::set<uint64_t> hits[end_pos - begin_pos] = {};

            auto delegate = [&hits, &it_zero_errors, begin_pos, &opt, textLength, new_overlap, &text](
                TIter it, auto const & /*read*/, unsigned const errors_spent)
            {
                uint64_t const bb = std::min(textLength - 1, begin_pos + opt.k_length - 1 + opt.k_length - new_overlap);
                if (errors_spent == 0)
                {
                    extend3<errors>(it, hits, it_zero_errors, errors - errors_spent, text, opt.k_length,
                        begin_pos + opt.k_length - new_overlap, begin_pos + opt.k_length - 1, // searched interval
                        begin_pos, bb, TDistanceTag() // entire interval
                    );
                }
                else
                {
                    extend(it, hits, errors - errors_spent, text, opt.k_length,
                        begin_pos + opt.k_length - new_overlap, begin_pos + opt.k_length - 1, // searched interval
                        begin_pos, bb, TDistanceTag() // entire interval
                    );
                }
            };

            auto const & needle = infix(text, begin_pos + opt.k_length - new_overlap, begin_pos + opt.k_length);
            TIter it(index);
            _optimalSearchScheme(delegate, it, needle, scheme, TDistanceTag());
            for (uint64_t j = begin_pos; j < end_pos; ++j)
            {
                value_type cValue = 0;
                uint64_t dist = 3 * errors;

                if(std::is_same<TDistanceTag, EditDistance>::value){
                    std::set<uint64_t> & chits = hits[j - begin_pos];
                    if(chits.empty())
                        continue;
                    //locate
                    typedef Pair <uint32_t, uint64_t>       TOcc;
                    std::vector<TOcc > occs;
                    occs.reserve(chits.size());
                    for(std::set<uint64_t>::iterator cit = chits.begin(); cit != chits.end(); ++cit){
                        occs.push_back(index.fwd.sa[*cit]);
                    }


                    std::sort(occs.begin(), occs.end(),
                              [](TOcc & x, TOcc & y){
                            if(getSeqNo(x) == getSeqNo(y))
                                return getSeqOffset(x) < getSeqOffset(y);
                            else
                                return getSeqNo(x) < getSeqNo(y);
                    });
/*
                    // check for similar occ
                    occs.erase(std::unique(occs.begin(), occs.end(),
                                           [&dist](TOcc & x, TOcc & y){
                                               return (getSeqNo(x) == getSeqNo(y) &&
                                               getSeqOffset(x) + dist >= getSeqOffset(y) &&
                                               getSeqOffset(x) <= getSeqOffset(y) + dist);
                                        }), occs.end()); //TODO use 3*errors instead*/


                    ++cValue;
                    uint64_t prev = 0;
                    for(uint64_t i = 1; i < occs.size(); ++i){
//                         std::cout << occs[i] << "\n";
                        if(! (getSeqNo(occs[i]) == getSeqNo(occs[prev]) &&
                           getSeqOffset(occs[i]) + dist >= getSeqOffset(occs[prev]) &&
                           getSeqOffset(occs[i]) <= getSeqOffset(occs[prev]) + dist))
                        {
                            prev = i;
                            ++cValue;
//                             std::cout << "prev: " << prev << "\n";
                        }
                    }
//                     std::cout << "\n";

                    cValue = std::min((uint64_t) cValue, max_val);


//                     std::vector<SHit> & chits = hits[j - begin_pos];
//                     std::sort(chits.begin(), chits.end(), sHit_smaller);
//                     chits.erase(std::unique(chits.begin(), chits.end(), sHit_similar<3 * errors>), chits.end()); //TODO use 3*errors instead
//                     cValue = std::min((uint64_t) chits.size(), max_val);
                }
                else
                {
                    cValue = hits[j - begin_pos].size();
                }

                if (countOccurrences(it_zero_errors[j - begin_pos]) > 1) // guaranteed to exist, since there has to be at least one match!
                {
                    for (auto const & occ : getOccurrences(it_zero_errors[j-begin_pos], Fwd()))
                    {
                        auto const occ_pos = posGlobalize(occ, limits);
                        c[occ_pos] = cValue;//hits[j - begin_pos];
                    }
                }
                else
                {
                    c[j] = cValue;//hits[j - begin_pos];
                }
            }
        }
    }

    resetLimits(indexText(index), c, opt.k_length);
}
