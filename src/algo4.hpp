using namespace seqan;
/*
struct smallHit{
    Pair <uint64_t, uint64_t> occ;
    uint8_t errors;
//     DnaString read;
};

bool occ_s(const smallHit & x, const smallHit & y)
{
    if(x.occ.i1 == y.occ.i1){
        if(x.occ.i2 == y.occ.i2){
             return x.errors < y.errors;
        }else{
            return x.occ.i2 < y.occ.i2;
        }
    }else{
        return x.occ.i1 < y.occ.i1;
    }
}

template<int disT>
bool occ_sim(const smallHit & x, const smallHit & y)
{
    return(x.occ.i1 == y.occ.i1 && x.occ.i2 + disT >= y.occ.i2 && x.occ.i2 <= y.occ.i2 + disT);
}

template<typename TOcc>
bool occs_smaller(TOcc & x, TOcc & y)
{
    if(getSeqNo(x) == getSeqNo(y))
        return getSeqOffset(x) < getSeqOffset(y);
    else
        return getSeqNo(x) < getSeqNo(y);
}*/

template <unsigned errors, typename TIndex, typename TText, typename TSeqLengths, typename TContainer, typename TOptions, typename TDistanceTag>
inline void runAlgo4(TIndex & index, TText const & text, TSeqLengths const & sequenceLengths, TContainer & c, TOptions const & opt, TDistanceTag const &)
{
    typedef typename TContainer::value_type value_type;

    constexpr uint64_t max_val = std::numeric_limits<value_type>::max();
    typedef Iter<TIndex, VSTree<TopDown<> > > TIter;

    uint64_t const textLength = length(text); // lengthSum() forwards to length() for a single string

    //     auto const & limits = stringSetLimits(text);
//     std::vector<uint64_t> sequenceLengths = getSeqLengths<uint64_t, uint64_t>(text); //has to be calc before since concat(text)
    /*
    std::cout << "Seqlengths\n";
    for(int i = 0; i < sequenceLengths.size(); ++i)
        std::cout << (int)sequenceLengths[i] << "\n";
    std::cout << "Finished Seq len\n";*/

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
//             std::vector<smallHit> hits[end_pos - begin_pos] = {};
            std::set<uint64_t> hits[end_pos - begin_pos] = {};
            /*
            auto delegate = [&hits, &it_zero_errors, begin_pos, &opt, textLength, new_overlap, &text](
                TIter it, auto const & hits, auto errors_spent)
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
             */
            auto const & needle = infix(text, begin_pos + opt.k_length - new_overlap, begin_pos + opt.k_length);
            TIter it(index);
            //_optimalSearchScheme(delegate, it, needle, scheme, TDistanceTag());
            for (uint64_t j = begin_pos; j < end_pos; ++j)
            {
                value_type cValue = 0;
                uint64_t dist = 3 * errors;

                std::set<uint64_t> & chits = hits[j - begin_pos];
                if(chits.empty()){
                    c[j] = 0;
                    continue;
                }

//                     std::sort(occs.begin(), occs.end(), occ_s);
//                     occs.erase(std::unique(occs.begin(), occs.end(), occ_sim<3 * errors>), occs.end());
//                     cValue = std::min((uint64_t) occs.size(), max_val);

                if(std::is_same<TDistanceTag, EditDistance>::value){
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

                    ++cValue;
                    uint64_t prev = 0;
                    for(uint64_t i = 1; i < occs.size(); ++i){
//                         std::cout << occs[i] << "\n";
                        if(!(getSeqNo(occs[i]) == getSeqNo(occs[prev]) &&
                           getSeqOffset(occs[i]) + dist >= getSeqOffset(occs[prev]) &&
                           getSeqOffset(occs[i]) <= getSeqOffset(occs[prev]) + dist))
                        {
                            prev = i;
                            ++cValue;
                        }
                    }
//                     std::cout << "\n";

                    cValue = std::min((uint64_t) cValue, max_val);


                }
                else
                {
                    cValue = hits[j - begin_pos].size();
                }

//                 std::sort(chits.begin(), chits.end(), occ_s);
//                 chits.erase(std::unique(chits.begin(), chits.end(), occ_sim<3 * errors>), chits.end());
//                 cValue = std::min((uint64_t) chits.size(), max_val);


                if (countOccurrences(it_zero_errors[j - begin_pos]) > 1) // guaranteed to exist, since there has to be at least one match!
                {
                    for (auto const & occ : getOccurrences(it_zero_errors[j-begin_pos], Fwd()))
                    {
//                         auto const occ_pos = posGlobalize(occ, limits);
//                         std::cout << "Occ: "<< occ << "\t" << (int)getSeqOffset(occ) << "\t" << (int)sequenceLengths[getSeqNo(occ)] << "\n";
                        uint64_t const occ_pos = getSeqOffset(occ) + sequenceLengths[getSeqNo(occ)];
//                         std::cout << (int)occ_pos << "\n";
                        c[occ_pos] = cValue;//hits[j - begin_pos];
                    }
//                     std::cout << "filled all positions\n";
                }
                else
                {
                    c[j] = cValue;//hits[j - begin_pos];
                }

            }
        }
    }

    resetLimits(text, sequenceLengths, c, opt.k_length);
}
