#ifndef OSSBESTX_H_
#define OSSBESTX_H_

#include "common.h"

using namespace seqan;

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
    if(strata == 99)
        strata = maxErrors;

    switch (maxErrors)
    {
        case 0:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;

        }
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


template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TVector, typename TVSupport,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void findChris(const int minErrors,
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
    if(strata == 99)
        strata = maxErrors;

    switch (maxErrors)
    {
        case 0:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;

        }
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
                    find<2, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag()); // MinError = 1
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
                    find<2, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 2 :
                {
                    find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
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
        case 5:
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
                    find<5, 5>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 5>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 2 :
                {
                    find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 5>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 3 :
                {
                    find<0, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 5>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 4 :
                {
                    find<0, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<5, 5>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 5 :
                {
                    find<0, 5>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;
        }
        case 6:
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
                    find<5, 5>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<6, 6>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 5>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<6, 6>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 2 :
                {
                    find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 5>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<6, 6>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 3 :
                {
                    find<0, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 6>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 5>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 4 :
                {
                    find<0, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<5, 6>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<5, 5>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 5 :
                {
                    find<0, 5>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<6, 6>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 6 :
                {
                    find<0, 6>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
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

        directsearchblockoffset = 2;
        directsearch_th = 4;

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
//     bool bestXMapper = false; //still needed multiple searches
    bool delayITV = false; //TODO add as option
    bool saFilter = false;
//     bool oneSSBestXMapper = false;
    uint32_t itvOccThreshold = 10;

    uint8_t maxError;
    uint8_t strata;
    uint32_t readLength;
    uint32_t numberOfSequences;

    typedef MapperTraits<TSpec, TConfig>        TTraits;
    typedef typename TTraits::TReadsContext     TReadsContext;
    typedef typename TTraits::TMatchesAppender  TMatches;
    typedef typename TTraits::TMatch            TMatch;
    typedef typename TTraits::TContigSeqs       TContigSeqs;
    typedef typename TTraits::TReadSeqs         TReadSeqs;
//     typedef typename TTraits::TBiIter           MySparseIter;
    typedef typename TTraits::TBitvectorsMeta   TBitvectorsMeta;


    // Shared-memory read-write data.
    TReadsContext &     ctx;
    TMatches &          matches;
    ReadsContextOSS<TSpec, TConfig>     ctxOSS;

    //track occ count per read
//     std::vector<uint32_t> /*&*/ readOccCount;

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

    void loadInputParameters(uint8_t inMaxError, uint8_t inStrata, uint32_t inReadLength, uint32_t inNumberOfSequences){
        maxError = inMaxError;
        strata = inStrata;
        readLength = inReadLength;
        numberOfSequences = inNumberOfSequences;
    }

    template <size_t nbrBlocks, typename TSALength>
    bool itvCondition(OptimalSearch<nbrBlocks> const & s,
                      uint8_t const blockIndex,
                      TSALength ivalOne)
    {
        return(itv && ivalOne < itvOccThreshold/*(static_cast<int>(s.pi.size()) - blockIndex - 1 + normal.directsearchblockoffset) * normal.directsearch_th*/);
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
        return(itv && iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1 < itvOccThreshold/*< (s.pi.size() - blockIndex - 1 + comp.directsearchblockoffset) * comp.directsearch_th*/);
    }

    template <typename TSALength>
    bool itvConditionUni(uint8_t const blockSize,
                         uint8_t const blockIndex,
                         TSALength ivalOne)
    {
        return(itv && ivalOne < itvOccThreshold/*(static_cast<int>(blockSize) - blockIndex - 1 + uni.directsearchblockoffset) * uni.directsearch_th*/);
    }

    template<size_t nbrBlocks>
    bool inBlockCheckMappabilityCondition(uint32_t needleLeftPos,
                                          uint32_t needleRightPos,
                                           OptimalSearch<nbrBlocks> const & s,
                                          uint8_t blockIndex)
    {
        uint32_t step = (needleRightPos - needleLeftPos - 1);
        if(normal.distancetoblockend > step)
            return false;
        uint32_t prevBlocklength = (blockIndex > 0) ? s.blocklength[blockIndex - 1] : 0;
        uint32_t nextBlocklength = s.blocklength[blockIndex];

        bool enoughDistanceToBlockEnds = step + normal.distancetoblockend < nextBlocklength && step - normal.distancetoblockend > prevBlocklength;
        return(((step & normal.step) == 0) && enoughDistanceToBlockEnds);
    }
};

#endif
