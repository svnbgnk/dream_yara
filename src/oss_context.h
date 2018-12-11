#ifndef OSSBESTX_H_
#define OSSBESTX_H_

#include "common.h"

using namespace seqan;

struct SparseVDesc{
    seqan::Pair<uint32_t, uint32_t> range;
    uint32_t smaller;
    uint32_t repLen;
    Dna lastChar;
};


struct SparseIter{
    seqan::Pair<uint32_t, uint32_t> fwdRange;
    uint32_t revRangeStart;
    uint32_t repLen;
    /*
    SparseVDesc fwd;
    SparseVDesc rev;
    SparseVDesc fwdp;
    SparseVDesc revp;*/
};

template<typename TIter2>
struct State{
    TIter2 it;
    uint32_t nlp;
    uint32_t nrp;
    uint8_t sId;
    uint8_t blockIndex;
    bool fwdDirection;


    //TODO use this only when TIter2 == MyIter
    template<typename TIter>
    State(TIter inIt,
          uint32_t nlp,
          uint32_t nrp,
          uint8_t sId,
          uint8_t blockIndex,
          bool fwdDirection) :
        it(inIt),
        nlp(nlp),
        nrp(nrp),
        sId(sId),
        blockIndex(blockIndex),
        fwdDirection(fwdDirection)
    {}
    /*
    State(TIter inIt,
          uint32_t nlp,
          uint32_t nrp,
          uint8_t sId,
          uint8_t blockIndex,
          bool fwdDirection) :
//         it(inIt),
        nlp(nlp),
        nrp(nrp),
        sId(sId),
        blockIndex(blockIndex),
        fwdDirection(fwdDirection)
    {



        it.fwdRange = inIt.fwdIter.vDesc.range;
        it.revRangeStart = inIt.revIter.vDesc.range.i1;
        it.repLen = inIt.fwdIter.vDesc.repLen;
    }*/
};

/*
        it.fwd.range = inIt.fwdIter.vDesc.range;
        it.fwd.smaller = inIt.fwdIter.vDesc.smaller;
        it.fwd.repLen = inIt.fwdIter.vDesc.repLen;
        it.fwd.lastChar = inIt.fwdIter.vDesc.lastChar;

        it.rev.range = inIt.revIter.vDesc.range;
        it.rev.smaller = inIt.revIter.vDesc.smaller;
        it.rev.repLen = inIt.revIter.vDesc.repLen;
        it.rev.lastChar = inIt.revIter.vDesc.lastChar;

        it.fwdp.range = inIt.fwdIter._parentDesc.range;
        it.fwdp.smaller = inIt.fwdIter._parentDesc.smaller;
        it.fwdp.repLen = inIt.fwdIter._parentDesc.repLen;
        it.fwdp.lastChar = inIt.fwdIter._parentDesc.lastChar;

        it.revp.range = inIt.revIter._parentDesc.range;
        it.revp.smaller = inIt.revIter._parentDesc.smaller;
        it.revp.repLen = inIt.revIter._parentDesc.repLen;
        it.revp.lastChar = inIt.revIter._parentDesc.lastChar;*/


template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TVector, typename TVSupport,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 int strata,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 std::vector<std::pair<TVector, TVSupport>> & bitvectors, // cant be const since TVSupport.set_vector(&TVector
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 TDistanceTag const & )
{
    if(strata == 99 || ossContext.oneSSBestXMapper)
        strata = maxErrors;

    switch (maxErrors)
    {
        case 1:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<1, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;

        }
        case 2:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<1, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 2 :
                {
                    find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;

        }
        case 3:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<1, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 2 :
                {
                    find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 3 :
                {
                    find<0, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;

        }
        case 4:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<1, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 2 :
                {
                    find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 3 :
                {
                    find<0, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 4 :
                {
                    find<0, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;
        }
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                 exit(1);
    }
}

//no bitvectors so substitute them
template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 const int strata,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 TDistanceTag const & )
{
    std::vector<std::pair<TBitvector, TSupport>> emtpy_bitvectors;
    find(minErrors, maxErrors, strata, ossContext, delegate, delegateDirect, index, emtpy_bitvectors, needles, TDistanceTag());
}

struct modusParameters{
public:
    bool nomappability;
    bool directsearch;
    bool compmappable;
    bool suspectunidirectional;

    bool testflipdensity;
    uint32_t step;
    uint32_t distancetoblockend;
    uint32_t directsearch_th;
    uint32_t directsearchblockoffset;
    float filter_th;
    float invflipdensity;
    uint32_t intervalsize;

    modusParameters(){
        setdefault();
    }

    void setdefault(){
        nomappability = true;
        directsearch = true;
        compmappable = true;
        suspectunidirectional = true;

        testflipdensity = true;
        //binaryNumber //has to be 2^x - 1 for fast modulo calculation
        step = 0b11;
        distancetoblockend = 2;

        directsearchblockoffset = 0;
        directsearch_th = 2;
        filter_th = 0.5;

        invflipdensity = 0.5;

        intervalsize = 3;
    }

    void print(){
        std::cout << "Cases Enabled: " << "\n";
        std::cout << nomappability << " " << directsearch << " " << compmappable << " " << suspectunidirectional << "\n";
        std::cout << "Params: " << "\n";

        std::cout << "step: " << step << "\n";
        std::cout << "distancetoblockend: " << distancetoblockend << "\n";
        std::cout << "directsearchblockoffset: " << directsearchblockoffset << "\n";
        std::cout << "directsearch_th: " << directsearch_th << "\n";
        std::cout << "filter_th: " << filter_th << "\n";
        std::cout << "invflipdensity: " << invflipdensity << "\n";
        std::cout << "intervalsize: " << intervalsize << "\n";
    }
};

// this struct carries Traits of Mapper, the reference to matches, text, read context and will
// contain conditions and parameters for OSS
template <typename TSpec, typename TConfig>
struct OSSContext
{
    //Parameters
    modusParameters normal;
    modusParameters comp;
    modusParameters uni;

    bool filterDelegate = true;
    bool trackReadCount = false;
    bool itv = true;
    bool bestXMapper = false; //still needed multiple searches
    bool oneSSBestXMapper = false;

    uint8_t maxError;
    uint8_t strata;
    uint32_t readLength;

    typedef MapperTraits<TSpec, TConfig>        TTraits;
    typedef typename TTraits::TReadsContext     TReadsContext;
    typedef typename TTraits::TMatchesAppender  TMatches;
    typedef typename TTraits::TMatch            TMatch;
    typedef typename TTraits::TContigSeqs       TContigSeqs;
    typedef typename TTraits::TReadSeqs         TReadSeqs;
    typedef typename TTraits::TSparseIter       TSparseIter;
//     typedef typename TTraits::TBiIter           MySparseIter;
    typedef State<TSparseIter>                  TTState;
    typedef typename TTraits::TBitvectorsMeta   TBitvectorsMeta;


    // Shared-memory read-write data.
    TReadsContext &     ctx;
    TMatches &          matches;
    ReadsContextOSS<TSpec, TConfig>     ctxOSS;
    std::vector<std::vector<TTState> > states;

    //track occ count per read
    std::vector<uint32_t> /*&*/ readOccCount;

    // Shared-memory read-only data.
    TReadSeqs & readSeqs;
    TContigSeqs const & contigSeqs;
    TBitvectorsMeta bitvectorsMeta;
//     Options const &     options;

    OSSContext(TReadsContext & ctx,
                 TMatches & matches,
                 TReadSeqs & readSeqs,
                 TContigSeqs const & contigSeqs) :
        ctx(ctx),
        matches(matches),
        readSeqs(readSeqs),
        contigSeqs(contigSeqs)
    {
        ;
    }

    void setdefault(){
        normal.setdefault();
        comp.setdefault();
        uni.setdefault();
    }

    void print(){
        std::cout << "Normal: ";
        normal.print();
        std::cout << "Comp: ";
        comp.print();
        std::cout << "Uni: ";
        uni.print();
    }

    void setReadContextOSS(uint8_t nerrors, uint8_t instrata, bool mScheme = false){
        maxError = nerrors;
        strata = instrata;

//         initReadsContext(ctx, readCount);
        std::cout << "maxError: " << (int)nerrors << "\tStrata: " << (int)strata << "\n";
        if(!mScheme){
            clear(ctxOSS);
            resize(ctxOSS, readSeqs);
            std::cout << "Using one Scheme" << "\n";
            oneSSBestXMapper = true;

            std::vector<TTState> v;
            for(int i = 0; i < maxError + 1; ++i)
                states.push_back(v);
    //         states.reserve(maxError);
        }else{
            std::cout << "Using multiple Schemes" << "\n";
            bestXMapper = true;
        }
    }

    template <typename TIter>
    inline void saveState(TIter & iter, uint32_t nlp, uint32_t nrp, uint8_t sid, uint8_t blockIndex, bool right, uint8_t errors){
        TTState state(iter, nlp, nrp, sid, blockIndex, right);
        states[errors].push_back(state);
//         int sizee = states[errors].size() - 1;
//         std::cout << (int)states[errors][sizee].blockIndex << "\n";
//         std::cout << "saved State" << "\n";
    }

    template <size_t nbrBlocks>
    bool itvCondition(OptimalSearch<nbrBlocks> const & s,
                      uint8_t const blockIndex,
                      uint32_t ivalOne)
    {
        return(itv && ivalOne < (static_cast<int>(s.pi.size()) - blockIndex - 1 + normal.directsearchblockoffset) * normal.directsearch_th);
    }


    template <typename TIter,
              size_t nbrBlocks>
    bool itvConditionComp(TIter iter,
                          uint32_t const needleLeftPos,
                          uint32_t const needleRightPos,
                          uint8_t const errors,
                          OptimalSearch<nbrBlocks> const & s,
                          uint8_t const blockIndex)
    {
        return(itv && iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1 < (s.pi.size() - blockIndex - 1 + comp.directsearchblockoffset) * comp.directsearch_th);
    }

    bool itvConditionUni(uint8_t const blockSize,
                         uint8_t const blockIndex,
                         uint32_t ivalOne)
    {
        return(itv && ivalOne < (static_cast<int>(blockSize) - blockIndex - 1 + uni.directsearchblockoffset) * uni.directsearch_th);
    }

    template<size_t nbrBlocks>
    bool inBlockCheckMappabilityCondition(uint32_t needleLeftPos,
                                          uint32_t needleRightPos,
                                           OptimalSearch<nbrBlocks> const & s,
                                          uint8_t blockIndex)
    {
        uint32_t prevBlocklength = (blockIndex > 0) ? s.blocklength[blockIndex - 1] : 0;
        uint32_t nextBlocklength = s.blocklength[blockIndex];
        uint32_t step = (needleRightPos - needleLeftPos - 1);


        bool enoughDistanceToBlockEnds = step + normal.distancetoblockend < nextBlocklength && step - normal.distancetoblockend > prevBlocklength;
        return(((step & normal.step) == 0) && enoughDistanceToBlockEnds);
    }

    void delegate(auto const & iter, uint32_t const needleId, uint8_t errors, bool const rev){
        std::cout << "Trying to report occ" << "\n";
        uint32_t occLength = repLength(iter);
        for (auto occ : getOccurrences(iter)){
            TMatch hit;
            setContigPosition(hit, occ, posAdd(occ, occLength));
            hit.errors = errors;
            setReadId(hit, readSeqs, needleId);
            setMapped(ctx, needleId);
//             setMinErrors(ctx, needleId, errors);
    //         hit.readId = needleId;
//TODO to du maybe use this form
    //         THit hit = { range(indexIt), (TSeedId)position(seedsIt), errors };

            appendValue(matches, hit, Generous(), typename TTraits::TThreading());
        }
    }
};

#endif
