using namespace seqan;

struct smallHit{
    Pair <uint32_t, uint32_t> occ;
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

template <unsigned errors, typename TIndex, typename TText, typename TContainer, typename TOptions>
inline void runAlgoTrivial(TIndex & index, TText const & text, TContainer & c, TOptions const & opt)
{
    typedef typename TContainer::value_type value_type;

    constexpr uint64_t max_val = std::numeric_limits<typename TContainer::value_type>::max();
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, opt.k_length);
//     _optimalSearchSchemeComputeChronBlocklength(scheme);

    uint64_t textLength = seqan::length(text);

//1000000
    //textLength/(opt.threads*1000)
    #pragma omp parallel for schedule(dynamic, 500) num_threads(opt.threads)
    for (uint64_t i = 0; i < textLength - opt.k_length + 1; ++i)
    {
        value_type hits = 0;

        std::vector<smallHit> myhits;
/*
        auto delegate = [&hits](auto const &it, auto const & , unsigned const ) {
            if ((uint64_t) hits + countOccurrences(it) <= max_val)
                hits += countOccurrences(it);
            else
                hits = max_val;
        };*/

        auto delegate = [&myhits](auto const &it, auto const & /*needle*/, unsigned const occErrors) {
           for (auto occ : getOccurrences(it.fwdIter)){
                smallHit me;
                me.occ = occ;
                me.errors = occErrors;
                myhits.push_back(me);
            }
        };

        auto const & needle = infix(text, i, i + opt.k_length);
        Iter<TIndex, VSTree<TopDown<> > > it(index);

        if(opt.indels)
            _optimalSearchScheme(delegate, it, needle, scheme, EditDistance());
        else
            _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());

        std::sort(myhits.begin(), myhits.end(), occ_s);
        myhits.erase(std::unique(myhits.begin(), myhits.end(), occ_sim<15>), myhits.end());

        if(myhits.size() < max_val)
            c[i] = myhits.size();
        else
            c[i] = max_val;
//         c[i] = hits;
    }

    resetLimits(indexText(index), c, opt.k_length);
}
