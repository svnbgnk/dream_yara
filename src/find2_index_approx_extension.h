#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_

#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <numeric>
#include <limits>
// #include "find2_index_approx_unidirectional.h"




namespace seqan{

template<typename TContigsSum>
struct SARange
{
    Pair<TContigsSum, TContigsSum> range;
    Pair<TContigsSum, TContigsSum> revRange;
    Pair<TContigsSum, TContigsSum> parentRange;
    Pair<TContigsSum, TContigsSum> parentRevRange;
    Pair<int8_t, int8_t> limOffsets;
    Dna lastChar;
    uint32_t repLength;
    uint8_t errors;
    bool goRight;
//     uint32_t shift = 0; //for itv
};



template<typename TSpec, typename TConfig,
         typename TDelegate,
         typename TContigsSum,
         typename TUniIter>
inline void fm_tree(OSSContext<TSpec, TConfig> & ossContext,
                    TDelegate & delegate,
                    uint32_t const needleId,
                    seqan::SARange<TContigsSum> & saRange,
                    TUniIter it,
                    bool const goRightDir,
                    uint32_t offset,
                    uint8_t const skipLastLayer,
                    uint32_t & counter)
{
    typedef MapperTraits<TSpec, TConfig>                     TTraits;
    typedef typename TTraits::TContigsPos                    TContigsPos;

    if((it.vDesc.range.i2 - it.vDesc.range.i1) < ossContext.fmTreeBreak){
//         ossContext.fmtreeBreakLocates += (it.vDesc.range.i2 - it.vDesc.range.i1) * (ossContext.samplingRate - offset - 1) / 2;
        for(TContigsSum i = it.vDesc.range.i1; i < it.vDesc.range.i2; ++i)
        {
            bool found = false;
            uint32_t max = ossContext.samplingRate - offset - skipLastLayer;
            TContigsPos pos = it.index->sa(i, max, found);
            if(found){
                setSeqOffset(pos, getSeqOffset(pos) + offset);
                if(goRightDir){
                    setSeqOffset(pos, ossContext.sequenceLengths[getSeqNo(pos)] - getSeqOffset(pos) - saRange.repLength);
//                     std::cout << "FM-tree revt: " << pos << "\n";
                }
                else{
//                      std::cout << "FM-treet: " << pos << "\n";
                }
//                 ossContext.allPos += getSeqOffset(pos);
                delegate(ossContext, needleId, saRange, pos);
            }
        }
    }
    else
    {
        TContigsSum ssp = getRank(it.index->sa.sparseString.indicators, it.vDesc.range.i1 - 1);
        TContigsSum sep = getRank(it.index->sa.sparseString.indicators, it.vDesc.range.i2 - 1);
        for(TContigsSum i = ssp; i < sep; ++i)
        {
            TContigsPos pos = posAdd(getValue(it.index->sa.sparseString.values, i), offset);
//             ossContext.allPos += getSeqOffset(pos);
            if(goRightDir){
            //    std::cout << "rev: " << pos  << "\n";
                setSeqOffset(pos, ossContext.sequenceLengths[getSeqNo(pos)] - getSeqOffset(pos) - saRange.repLength);
//                 std::cout << "FM-tree rev: " << pos << "\toffset:" << offset << "\t" << i << "\n";
            }
            else
            {
//                 std::cout << "FM-tree: " << pos << "\toffset:" << offset << "\t" << i << "\n";
            }
            delegate(ossContext, needleId, saRange, pos);
        }

        if(offset + skipLastLayer + 1 < ossContext.samplingRate && goDown(it)){
//             counter += 4;

            do{
                fm_tree(ossContext, delegate, needleId, saRange, it, goRightDir, offset + 1, skipLastLayer, counter);
            }while(goRight(it));
        }
    }
}

template<typename TSpec, typename TConfig,
         typename TContigSeqs,
         typename TDelegate,
         typename TContigsSum,
         typename TIndex>
inline void locate(OSSContext<TSpec, TConfig> & ossContext,
                   TContigSeqs const & genome,
                   TDelegate & delegate,
                   uint32_t const needleId,
                   seqan::SARange<TContigsSum> & saRange,
                   Iter<TIndex, VSTree<TopDown<> > > iter)
{
    typedef MapperTraits<TSpec, TConfig>                     TTraits;
    typedef typename TTraits::TContigsPos                    TContigsPos;

    if(ossContext.fmTreeThreshold < saRange.range.i2 - saRange.range.i1 && ossContext.samplingRate > 1)
    {
//         ossContext.fmtreeLocates += saRange.range.i2 - saRange.range.i1;
        //locate Interval with fm_tree
            auto uniIter = (ossContext.earlyLeaf && saRange.goRight) ? iter.revIter : iter.fwdIter;
            uniIter.vDesc.repLen = saRange.repLength;
            uniIter.vDesc.range = (ossContext.earlyLeaf && saRange.goRight) ? saRange.revRange : saRange.range;

            uint32_t counter = 0;
            //early leaf node calculation (corresponds to the last level of the fm_tree)
            if(ossContext.earlyLeaf){
                TContigsSum ssp = (saRange.goRight) ?
                    getRank(uniIter.index->sa.sparseString.indicators, saRange.parentRevRange.i1 - 1) :
                    getRank(uniIter.index->sa.sparseString.indicators, saRange.parentRange.i1 - 1);
                TContigsSum sep = (saRange.goRight) ?
                    getRank(uniIter.index->sa.sparseString.indicators, saRange.parentRevRange.i2 - 1) :
                    getRank(uniIter.index->sa.sparseString.indicators, saRange.parentRange.i2 - 1);
                auto lastChar = saRange.lastChar;

                if (saRange.goRight){
//                     std::cout << "\nrevIndex\n";

                    for(TContigsSum i = ssp; i < sep; ++i){
                        TContigsPos pos = getValue(uniIter.index->sa.sparseString.values, i);
                        setSeqOffset(pos, getSeqOffset(pos) - 1);
                        //Calculate position in forward text
                        setSeqOffset(pos, ossContext.sequenceLengths[getSeqNo(pos)] - getSeqOffset(pos) - 1);
                        if (genome[getSeqNo(pos)][getSeqOffset(pos)] == lastChar){
                            setSeqOffset(pos, getSeqOffset(pos) - (saRange.repLength - 1));
//                             std::cout << "FPos calc from reverseIndex: " << pos << "\n";
                            delegate(ossContext, needleId, saRange, pos);
                        }
                    }

                        //Debugging
                        for (TContigsSum i = saRange.range.i1; i < saRange.range.i2; ++i)
                        {
                            TContigsPos pos2 = iter.fwdIter.index->sa[i];
//                             std::cout << "All forward Pos: " << pos2 << "\t(" << i << ")\n";
                        }

                }else{
//                     std::cout << "\nfwdIndex\n";
                    for(TContigsSum i = ssp; i < sep; ++i){
                        TContigsPos pos = posAdd(getValue(uniIter.index->sa.sparseString.values, i), (-1));

                        if (genome[getSeqNo(pos)][getSeqOffset(pos)] == lastChar){
                            //Report position
//                             std::cout << "FPos: " << pos << "\n";
                            delegate(ossContext, needleId, saRange, pos);
                        }
                    }
                }
            }
         /*
            //Debugging
            for (TContigsSum i = uniIter.vDesc.range.i1; i < uniIter.vDesc.range.i2; ++i)
                        {
                TContigsPos pos2 = iter.fwdIter.index->sa[i];
                std::cout << "All forward Pos: " << pos2 << "\t(" << i << ")\n";
            }*/

            bool goRight2 = ossContext.earlyLeaf && saRange.goRight;
            fm_tree(ossContext, delegate, needleId, saRange, uniIter, goRight2, 0, static_cast<uint8_t>(ossContext.earlyLeaf), counter);
//             ossContext.fmtreeBacktrackings += 2 * counter;
    }
    else
    {
//         ossContext.defaultLocates += saRange.range.i2 - saRange.range.i1;
        for (uint32_t i = saRange.range.i1; i < saRange.range.i2; ++i)
        {
            TContigsPos pos = iter.fwdIter.index->sa[i];
//             ossContext.allPos += getSeqOffset(pos);
            delegate(ossContext, needleId, saRange, pos);
        }
    }
}

template<typename TContex,
         typename TContigSeqs,
         typename TDelegate,
         typename TIndex,
         typename TContigsSum>
inline void locateSARanges(TContex & ossContext,
                           TContigSeqs const & genome,
                           TDelegate & delegate,
                           uint32_t const needleId,
                           Iter<TIndex, VSTree<TopDown<> > > iter,
                           std::vector<seqan::SARange<TContigsSum> > & ranges)
{
    //Sort Ranges
    std::sort(ranges.begin(), ranges.end(), [](SARange<TContigsSum> & x, SARange<TContigsSum> & y){
        if(x.range.i1 == y.range.i1)
        {
            return(x.range.i2 > y.range.i2);
        }
        else
        {
            return(x.range.i1 < y.range.i1);
        }
    });


//     TContigsSum mymax = std::numeric_limits<TContigsSum>::max();


    //TODO choose parameter
//     int th_fm_tree = (pow(10, 3) * 4) ;   // th_fm_tree = pow(10, 3) * 4; //int expectedLF = pow(10, 4) * 2;
    uint64_t in = 0;
    while(in < ranges.size()){
        auto crange = ranges[in];
        ++in;
        if(in < ranges.size()){
            while(crange.range.i2 >= ranges[in].range.i1){
                if(crange.range.i2 < ranges[in].range.i2){
                    crange.range.i2 = ranges[in].range.i2;
                    if(crange.repLength < ranges[in].repLength)
                        //TODO update limOffsets of crange
                        crange.repLength = ranges[in].repLength;
                    if(crange.errors < ranges[in].errors)
                        crange.errors = ranges[in].errors;
                }
                ++in;
                //also stop if we scanned over all intervals
                if(in >= ranges.size())
                    break;
            } // merge intervals
        }

        ossContext.delegateFilteredOcc += crange.range.i2 - crange.range.i1;
        locate(ossContext, genome, delegate, needleId, crange, iter);
    }
}


template<typename TVector, typename TVSupport,
         typename TSALength>
inline TVector & getTVector(std::vector<std::pair<TVector, TVSupport> > & bitvectors,
           Pair<uint8_t, Pair<TSALength, TSALength>> const & brange)
{
    return bitvectors[brange.i1].first;
}

template<typename TVector, typename TVSupport,
         typename TSALength>
inline TVSupport & getTVSupport(std::vector<std::pair<TVector, TVSupport> > & bitvectors,
           Pair<uint8_t, Pair<TSALength, TSALength>> const & brange)
{
    return bitvectors[brange.i1].second;
}

template<typename TVector, typename TVSupport,
         typename TSALength>
inline TVector & getTVector(std::vector<std::pair<TVector, TVSupport>* > & bitvectors,
           Pair<uint8_t, Pair<TSALength, TSALength>> const & brange)
{
    return bitvectors[brange.i1]->first;
}

template<typename TVector, typename TVSupport,
         typename TSALength>
inline TVSupport & getTVSupport(std::vector<std::pair<TVector, TVSupport>* > & bitvectors,
           Pair<uint8_t, Pair<TSALength, TSALength>> const & brange)
{
    return bitvectors[brange.i1]->second;
}

template <typename TBitvectorPair,
          typename TSALength>
inline void getConsOnes(std::vector<TBitvectorPair > & bitvectors,
                Pair<uint8_t, Pair<TSALength, TSALength> > & inside_bit_interval,
                uint32_t const intervalsize,
                std::vector<std::pair<TSALength, TSALength>> & consOnesOutput)
{
    auto & b = getTVector(bitvectors, inside_bit_interval);
    TSALength k = inside_bit_interval.i2.i1;
    TSALength startOneInterval = inside_bit_interval.i2.i1;
    while(k < inside_bit_interval.i2.i2){
        TSALength interval = 0;
        //TODO delete second condition it should end with 1
        while(b[k + interval] == 0 && (k + interval) < inside_bit_interval.i2.i2){
            ++interval;
        }
        if(interval >= intervalsize){
            consOnesOutput.push_back(std::make_pair(startOneInterval, k));
            startOneInterval = k + interval;
        }
        k += interval;
        interval = 0;
        ++k;
    }
    consOnesOutput.push_back(std::make_pair(startOneInterval, k));
}

template <typename TSpec, typename TConfig,
          typename TDelegate, typename TDelegateD,
          typename TIndex,
          typename TNeedle,
          typename TBitvectorPair,
          size_t nbrBlocks,
          typename TSALength,
          typename TDir,
          typename TDistanceTag>
inline void filter_interval(OSSContext<TSpec, TConfig> & ossContext,
                            TDelegate & delegate,
                            TDelegateD & delegateDirect,
                            Iter<TIndex, VSTree<TopDown<> > > iter,
                            Pair <int8_t, int8_t> limOffsets,
                            TNeedle const & needle,
                            uint32_t needleId,
                            std::vector<TBitvectorPair > & bitvectors,
                            uint32_t const needleLeftPos,
                            uint32_t const needleRightPos,
                            uint8_t const errors,
                            OptimalSearch<nbrBlocks> const & s,
                            uint8_t const blockIndex,
                            Pair<uint8_t, Pair<TSALength, TSALength>> & inside_bit_interval,
                            TDir const & ,
                            TDistanceTag const &)
{
    typedef typename TConfig::TContigsSum       TContigsSum; //SAVALUE
    std::vector<std::pair<TContigsSum, TContigsSum> > consOnes;
    getConsOnes(bitvectors, inside_bit_interval, ossContext.normal.intervalsize, consOnes);
    uint32_t nos = ossContext.numberOfSequences;//countSequences(*iter.fwdIter.index);

    for(uint32_t i = 0; i < consOnes.size(); ++i){
        if (std::is_same<TDir, Rev>::value){
            iter.revIter.vDesc.range.i1 = consOnes[i].first + nos;
            iter.revIter.vDesc.range.i2 = consOnes[i].second + nos;
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter.revIter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, false, Rev(), TDistanceTag());
        }
        else
        {
            iter.fwdIter.vDesc.range.i1 = consOnes[i].first + nos;
            iter.fwdIter.vDesc.range.i2 = consOnes[i].second + nos;
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter.fwdIter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, false, Fwd(), TDistanceTag());
        }
    }
}

template <typename TSpec, typename TConfig,
          typename TDelegateD,
          typename TNeedle,
          typename TContigSeqs,
          typename TSAValue,
          size_t nbrBlocks>
inline void genomeSearch(OSSContext<TSpec, TConfig> & ossContext,
                         TDelegateD & delegateDirect,
                         TNeedle const & needle,
                         uint32_t needleId,
                         uint8_t errors,
                         OptimalSearch<nbrBlocks> const & s,
                         uint8_t const blockIndex,
                         TContigSeqs const & genome,
                         TSAValue const & sa_info,
                         std::array<uint32_t, nbrBlocks> & blockStarts,
                         std::array<uint32_t, nbrBlocks> & blockEnds)
{
    typedef typename TConfig::TContigsLen                    TContigsLen;
    TContigsLen seqOffset = getSeqOffset(sa_info);
    for(uint32_t j = 0; j < nbrBlocks - blockIndex; ++j){
        // compare bases to needle

        for(uint32_t k = blockStarts[j]; k <  blockEnds[j]; ++k){
            if(needle[k] != genome[getSeqNo(sa_info)][seqOffset + k]){
                ++errors;
            }
        }
        if(errors < s.l[blockIndex + j] || errors > s.u[blockIndex + j]){
            return;
        }
    }
    delegateDirect(ossContext, sa_info, posAdd(sa_info, length(needle)), needleId, errors);
}

template<typename TBitvectorPair,
         typename TSALength>
inline bool checkSinglePos(std::vector<TBitvectorPair > & bitvectors,
                           Pair<uint8_t, Pair<TSALength, TSALength> > const & brange,
                           uint32_t offset)
{
    if(bitvectors.empty()){
        return true;
    }
    else
    {
        auto & b = getTVector(bitvectors, brange);
        return (b[brange.i2.i1 + offset] == 1);
    }
}

template<typename TSAValue, typename TContigsLen>
inline void saPosOnFwd(TSAValue & sa_info,
                       TContigsLen const genomelength,
                       uint32_t const occLength)
{

    setSeqOffset(sa_info, genomelength - getSeqOffset(sa_info) - occLength);
//     sa_info.i2 = genomelength - sa_info.i2 - occLength;
}


//TODO TSAValue to TContigPos
template <typename TContex,
          typename TDelegateD,
          typename TString,
          typename TContigsLen,
          typename TSAValue,
          typename TNeedle>
inline void alignmentMyersBitvector(TContex & ossContext,
                                    TDelegateD & delegateDirect,
                                    TNeedle const & needle,
                                    uint32_t needleId,
                                    TString const & n_infix,
                                    TString const & ex_infix,
                                    TContigsLen const genomelength,
                                    TSAValue const & sa_info,
                                    uint8_t max_e,
                                    uint8_t overlap_l,
                                    uint8_t overlap_r,
                                    uint8_t intDel,
                                    bool usingReverseText)
{
/*
    std::cout << needleId << "\t" << sa_info << "Needle: " << "\n";
    std::cout << "   " << needle << "\n";
    std::cout << "   " << n_infix << "\n";
    std::cout << ex_infix << "\n\n";*/

    uint32_t needleL = length(needle);
    uint32_t ex_infixL = needleL + overlap_l + overlap_r;

    int32_t initialScore = globalAlignmentScore(ex_infix, needle, MyersBitVector());

 //assume more Insertions (in the read) than deletions
    int32_t ins_initialScore = globalAlignmentScore(n_infix, needle, MyersBitVector());

    if(ins_initialScore >= 0 - 2 * max_e || initialScore >= 0 - overlap_l - overlap_r - max_e + intDel) //MM creates one error D creates one error since now it also align to overlap
    {
        TSAValue sa_info_tmp = sa_info;
        //No Insertions or Deletions
//         cout << "E: " << (int)0 << endl;
        TString const infix_tmp0 = infix(ex_infix, overlap_l, ex_infixL - overlap_r); //use n_infix
        int32_t errors2 = 0 - globalAlignmentScore(infix_tmp0, needle, MyersBitVector()); //shold be uint32_t

        if(usingReverseText){
            saPosOnFwd(sa_info_tmp, genomelength, needleL);
        }

        //if the current best hit is not good it will not get reported
        int32_t bestHit = (errors2 <= max_e) ? errors2 : max_e + 1;
        TSAValue best_sa = sa_info_tmp;
        uint32_t best_length = needleL;
//         TString bestInfix = infix_tmp0;

        for(uint8_t e = 1; e <= max_e; ++e){
            if(bestHit <= e){
                //trimHit
//                 trimHit(needle, bestInfix, best_sa, best_length, bestHit);
                delegateDirect(ossContext, best_sa, posAdd(best_sa, best_length), needleId, bestHit);
                return;
            }

//             cout << "E: " << (int)e << endl;
            for(uint8_t del = 0; del <= e; ++del){
                //del is number of deletions
                uint8_t ins = e - del; //number of insertions
                sa_info_tmp = sa_info;
                int32_t occLength;

                if(del > 1 && ins == 0 || ins > 1 && del == 0){
                //only insertion or deletions
                    int16_t pos = (ins > del) ? 1 : (-1);
                    int16_t m = std::max(del,ins);
                    for(int16_t k = 0; k <= m; ++k)
                    {
                        if(!(0 <= overlap_l + (pos * k) && overlap_r >= 0 - (pos * (m - k))))
                            continue;
                        sa_info_tmp = sa_info;
//                         sa_info_tmp.i2 = sa_info_tmp.i2 + (pos * k);
                        setSeqOffset(sa_info_tmp, getSeqOffset(sa_info_tmp) + (pos * k));

                        TString const infix_tmp = infix(ex_infix, (pos * k) + overlap_l, ex_infixL - overlap_r - (pos * (m - k)));
                        errors2 = 0 - globalAlignmentScore(infix_tmp, needle, MyersBitVector());
                        if(errors2 < bestHit){
                            occLength = length(needle) - (pos * m);
                            if(usingReverseText){
                                saPosOnFwd(sa_info_tmp, genomelength, occLength);
                            }
//                             std::cout << "2  " << infix_tmp << "\t" << errors2 << "\t" << sa_info_tmp << "\t" << posAdd(sa_info_tmp, occLength) << "\n";

//                             delegateDirect(ossContext, sa_info_tmp, posAdd(sa_info_tmp, occLength), needleId, errors2);
                            bestHit = errors2;
                            best_sa = sa_info_tmp;
                            best_length = occLength;
//                             bestInfix = infix_tmp;
                        }
                    }
                }
                else
                {
                    //insertions left and deletion right
                    if(overlap_l >= del){
                        TString const infix_tmp = infix(ex_infix, overlap_l - del, ex_infixL - overlap_r - ins);
                        sa_info_tmp = sa_info;
//                         sa_info_tmp.i2 = sa_info_tmp.i2 - del;
                        setSeqOffset(sa_info_tmp, getSeqOffset(sa_info_tmp) - del);
                        errors2 = 0 - globalAlignmentScore(infix_tmp, needle, MyersBitVector());
                        if(errors2 < bestHit){
                            occLength = length(needle) - ins + del;
                            if(usingReverseText){
                                saPosOnFwd(sa_info_tmp, genomelength, occLength);
                            }
//                             std::cout << "3  " << infix_tmp << "\t" << errors2 << "\t" << sa_info_tmp << "\t" << posAdd(sa_info_tmp, occLength) << "\n";
//                             delegateDirect(ossContext, sa_info_tmp, posAdd(sa_info_tmp, occLength), needleId, errors2);
                            bestHit = errors2;
                            best_sa = sa_info_tmp;
                            best_length = occLength;
//                             bestInfix = infix_tmp;
                        }
                    }

                    //insertions right and deletion left
                    if(overlap_r >= del){
                        sa_info_tmp = sa_info; //just include del from before into the calculation and delete this
                        TString const infix_tmp = infix(ex_infix, overlap_l + ins, ex_infixL - overlap_r + del);
                        errors2 = 0 - globalAlignmentScore(infix_tmp, needle, MyersBitVector());
//                         sa_info_tmp.i2 = sa_info_tmp.i2 + ins;
                        setSeqOffset(sa_info_tmp, getSeqOffset(sa_info_tmp) + ins);
                        if(errors2 < bestHit){
                            occLength = length(needle) - ins + del;
                            if(usingReverseText){
                                saPosOnFwd(sa_info_tmp, genomelength, occLength);
                            }
//                             std::cout << "4  " << infix_tmp << "\t" << errors2 << "\t" << sa_info_tmp << "\t" << posAdd(sa_info_tmp, occLength) << "\n";
//                             delegateDirect(ossContext, sa_info_tmp, posAdd(sa_info_tmp, occLength), needleId, errors2);
                            bestHit = errors2;
                            best_sa = sa_info_tmp;
                            best_length = occLength;
//                             bestInfix = infix_tmp;
                        }
                    }
                }
            }
        }
        if(bestHit <= max_e){
//             trimHit(needle, bestInfix, best_sa, best_length, bestHit);
            delegateDirect(ossContext, best_sa, posAdd(best_sa, best_length), needleId, bestHit);
        }
    }
}

template <typename TContex,
          typename TDelegateD,
          typename TString,
          typename TContigsLen,
          typename TSAValue,
          typename TNeedle>
inline void inTextVerificationN(TContex & ossContext,
                                TDelegateD & delegateDirect,
                                TNeedle & needle,
                                uint32_t needleId,
                                TString & ex_infix,
                                TContigsLen const genomelength,
                                TSAValue const & sa_info,
                                uint8_t max_e,
                                uint8_t upper,
                                uint8_t lower,
                                bool usingReverseText)
{
    typedef ModifiedString<TNeedle, ModReverse>           TNeedleInfixRev;
    typedef ModifiedString<TString, ModReverse>           TStringInfixRev;
    typedef Finder<TString>                               TFinder;
    typedef Finder<TStringInfixRev>                       TFinder2;
    typedef AlignTextBanded<FindInfix,
                            NMatchesNone_,
                            NMatchesNone_>                    TMyersSpecInfix;
    typedef Myers<TMyersSpecInfix, True, void>                TAlgorithmInfix;
    typedef PatternState_<TNeedle, TAlgorithmInfix>  TPatternInfix;
    TPatternInfix patternInfix;

    typedef AlignTextBanded<FindPrefix,
                            NMatchesNone_,
                            NMatchesNone_>               TMyersSpec;
    typedef Myers<TMyersSpec, True, void>                TAlgorithm;
    typedef PatternState_<TNeedle, TAlgorithm>  TPattern;
    TPattern pattern;
    typedef PatternState_<TNeedleInfixRev, TAlgorithm>  TPatternRev;
    TPatternRev patternRev;


    ++ossContext.itvAttemps;

    //calc Score:
    uint8_t minErrors = max_e + 1;
    TFinder finderInfix(ex_infix);

//     std::cout << "Score from: \n" << ex_infix << "\n" << needle << "\n";
    while (find(finderInfix, needle, patternInfix, -static_cast<int>(length(needle)))) //TODO check
    {
        uint16_t currentErrors = -getScore(patternInfix);
        if(minErrors > currentErrors)
            minErrors = currentErrors;
    }

    if(minErrors > max_e)
        return;

//     if(minErrors > upper)
//         return;

    // check if match is outside of current search -> other search will find it no reason to report it now
    if(minErrors > upper || minErrors < lower)
        return;

    TFinder finder(ex_infix);
    uint8_t mErrors = max_e * 4;
    TContigsLen endPos = 0;
    while (find(finder, needle, pattern, -static_cast<int>(max_e * 4))) //TODO choose correct value
    {
        int currentEnd = position(finder) + 1;
        uint16_t currentErrors = -getScore(pattern);
//         std::cout << currentErrors << "\t" << currentEnd << "\n";
        if (currentErrors < mErrors)
        {
            mErrors = currentErrors;
            endPos = currentEnd;
        }
    }
    TString infixPrefix = infix(ex_infix, 0, endPos);

//     std::cout << "Cut one: " << "\n" << infixPrefix << "\n";


    TStringInfixRev infixRev(infixPrefix);
    TNeedleInfixRev needleRev(needle);
    TFinder2 finder2(infixRev);

    mErrors = max_e * 3;
    TContigsLen startPos = endPos;

    while (find(finder2, needleRev, patternRev, -static_cast<int>(max_e * 3))) //TODO choose correct value
    {
        int currentEnd = position(finder2) + 1;
        uint16_t currentErrors = -getScore(patternRev);
//         std::cout << currentErrors << "\t" << currentEnd << "\n";
        if (currentErrors < mErrors)
        {
            mErrors = currentErrors;
            startPos = currentEnd;
        }
    }

    TSAValue sa_info_tmp = sa_info;

//     std::cout << "final cut" << needleId << "\t" << posAdd(sa_info_tmp, endPos - startPos) << "\t" << posAdd(sa_info_tmp, endPos) << "\n";
//     std::cout << infix(ex_infix, endPos - startPos, endPos) << "\n\n";
    if(usingReverseText){
            saPosOnFwd(sa_info_tmp, genomelength, length(needle));
    }

    delegateDirect(ossContext, posAdd(sa_info_tmp, endPos - startPos), posAdd(sa_info_tmp, endPos), needleId, minErrors);
}

template <typename TSpec, typename TConfig,
          typename TDelegateD,
          typename TIndex,
          typename TNeedle,
          typename TBitvectorPair,
          size_t nbrBlocks,
          typename TSALength,
          typename TDir,
          typename TDistanceTag>
inline void directSearch(OSSContext<TSpec, TConfig> & ossContext,
                         TDelegateD & delegateDirect,
                         Iter<TIndex, VSTree<TopDown<> > > iter,
                         Pair <int8_t, int8_t> & limOffsets,
                         TNeedle const & needle,
                         uint32_t needleId,
                         std::vector<TBitvectorPair > & bitvectors,
                         uint32_t const needleLeftPos,
                         uint32_t const needleRightPos,
                         uint8_t const errors,
                         OptimalSearch<nbrBlocks> const & s,
                         uint8_t const blockIndex,
                         Pair<uint8_t, Pair<TSALength, TSALength>> const & brange,
                         TDir const & ,
                         TDistanceTag const &)
{
    typedef MapperTraits<TSpec, TConfig>                     TTraits;
    typedef typename TTraits::TSA                            TSA;
    typedef typename Size<TSA>::Type                         TSAPos;
    typedef typename Value<TSA>::Type                        TSAValue;
    typedef typename TTraits::TContigsPos                    TContigsPos;
    typedef typename TConfig::TContigsLen                    TContigsLen; //sa.i2
    typedef typename TConfig::TContigsSize                   TContigsSize; //sa.i1
    typedef typename TConfig::TContigsSum                    TContigsSum; //SAVALUE
//     typedef typename TConfig::TAlloc                         TAlloc;
//     typedef SeqStore<void, YaraContigsConfig<TAlloc> >       TContigs;
//     typedef typename TContigs::TSeqs                         TContigSeqs;
    typedef typename TTraits::TContigSeqs                    TContigSeqs;
    typedef typename InfixOnValue<TContigSeqs const>::Type   TContigSeqsInfix;
    typedef typename InfixOnValue<TNeedle const>::Type       TNeedleInfix;
//      typedef ModifiedString<TNeedleInfix, ModReverse>     TNeedleInfixRev;

    typedef typename TTraits::TReadSeqs                                                 TReadSeqs;
    typedef typename Size<TReadSeqs>::Type                                              TReadId;

    TContigSeqs const & genome = ossContext.contigSeqs;

//     auto const & genome = indexText(*iter.fwdIter.index);

//     bool test = needle;
    ModifiedString<TNeedle, ModReverse > revneedle(needle);

    if (std::is_same<TDistanceTag, EditDistance>::value){
        //TODO put this into a function
        //TODO if we are only interested in the best hit call return after delegate calls
        uint32_t needleL = length(needle);
        uint8_t max_e = ossContext.maxError;
        //there are different search with different upper and lower -> something goes wrong?
        uint8_t upper = s.u[s.u.size() - 1];
        uint8_t lower = s.l[s.l.size() - 1];

//         uint8_t overlap_l = max_e;
//         uint8_t overlap_r = max_e;

//         std::cout << "NLP: " << needleLeftPos << "\tNRP: " << needleRightPos << "\trepL: " << (int)repLength(iter) << "\trange: " <<  (int)needleRightPos - needleLeftPos - 1 << "\n";
//         std::cout << "Off: " << (int)limOffsets.i1 - max_e << ", " << (int)limOffsets.i2 - max_e << "\n";

        int8_t overlap_l;
        int8_t overlap_r;

        if(ossContext.delayITV)
        {
            overlap_l = limOffsets.i1;
            overlap_r = limOffsets.i2;
        }
        else
        {
            overlap_l = max_e;
            overlap_r = max_e;
        }

        for(TContigsSum r = 0; r < iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1; ++r)
        {
            if(checkSinglePos(bitvectors, brange, r)){
                TSAValue sa_info;
                TContigsLen chromlength;

                if(std::is_same<TDir, Rev>::value){
                    sa_info = iter.fwdIter.index->sa[iter.fwdIter.vDesc.range.i1 + r];
                    TContigsLen seqOffset = getSeqOffset(sa_info);
                    chromlength = length(genome[getSeqNo(sa_info)]);
//                     std::cout << "checking forward: " << sa_info << "\t" << chromlength << "\n";
                    if(!(needleLeftPos + overlap_l <= seqOffset && chromlength - 1 >= seqOffset - needleLeftPos + needleL - 1 + overlap_r))
                        continue;
//                     sa_info.i2 = sa_info.i2 - needleLeftPos;
                    setSeqOffset(sa_info, seqOffset - needleLeftPos);
                }
                else
                {
                    sa_info = iter.revIter.index->sa[iter.revIter.vDesc.range.i1 + r];
                    TContigsLen seqOffset = getSeqOffset(sa_info);
                    chromlength = length(genome[getSeqNo(sa_info)]);
//                     std::cout << "checking reverse: " << sa_info << "\t" << chromlength << "\n";
                    if(!(chromlength - 1 >= seqOffset + needleRightPos - 1 + overlap_r && seqOffset + needleRightPos - 1 - overlap_l >= length(needle) + 1))
                        continue;
//                     sa_info.i2 = chromlength - sa_info.i2 - needleRightPos + 1;
                    setSeqOffset(sa_info, chromlength - seqOffset - needleRightPos + 1);
                }
                //update seqOffset
                TContigsLen seqOffset = getSeqOffset(sa_info);
/*
                //types for globalAlignmentScore
                TContigSeqsInfix ex_infix = infix(genome[getSeqNo(sa_info)], seqOffset - overlap_l, seqOffset + needleL + overlap_r);
                TContigSeqsInfix n_infix = infix(genome[getSeqNo(sa_info)], seqOffset, seqOffset + needleL);
                alignmentMyersBitvector(ossContext, delegateDirect, needle, needleId, n_infix, ex_infix, chromlength, sa_info, max_e, overlap_l, overlap_r, intDel, false);
                */

                TReadId readId = getReadId(ossContext.readSeqs, needleId);
                int8_t overlap_l_tmp = (overlap_l <= seqOffset) ? overlap_l : 0;

//                 std::cout << ossContext.strata + s.l[s.l.size() - 1] << "\t" << ossContext.maxError << "\n";
                //TODO search <0, 2> need to search till 0 is found or 1
                if(ossContext.delayITV && (isMapped(ossContext.ctx, readId && getMinErrors(ossContext.ctx, readId) + ossContext.strata <= upper) || ossContext.strata + s.l[s.l.size() - 1] >= ossContext.maxError))
                {
                    delegateDirect(ossContext, posAdd(sa_info, 0 - overlap_l_tmp), posAdd(sa_info, needleL + overlap_r), needleId, 127);
                }
                else
                {
                    TContigSeqsInfix ex_infix = infix(genome[getSeqNo(sa_info)], seqOffset - overlap_l, seqOffset + needleL + overlap_r);
                    inTextVerificationN(ossContext, delegateDirect, needle, needleId, ex_infix, chromlength, posAdd(sa_info, 0 - overlap_l_tmp), max_e, upper, lower, false);
                }
//                 std::cout << posAdd(sa_info, -overlap_l) << "\t" << posAdd(sa_info, needleL + overlap_r) << "\n";
            }
        }
    }
    else
    {
        std::array<uint32_t, nbrBlocks> blockStarts;
        std::array<uint32_t, nbrBlocks> blockEnds;
        std::copy(std::begin(s.blockStarts) + blockIndex, std::end(s.blockStarts), std::begin(blockStarts));
        std::copy(std::begin(s.blockEnds) + blockIndex, std::end(s.blockEnds), std::begin(blockEnds));

        if(std::is_same<TDir, Rev>::value){
            //modify blockstart in case we are still inside a block
            if(needleRightPos - 1 > blockStarts[0] && needleRightPos - 1 < blockEnds[0])
                blockStarts[0] = needleRightPos - 1;

            for(TSAPos i = 0; i < iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1; ++i){
//                 if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
                if(checkSinglePos(bitvectors, brange, i)){
                    // mappability information is in reverse index order if we use the forward index
                    TSAValue sa_info = iter.fwdIter.index->sa[iter.fwdIter.vDesc.range.i1 + i];
                    TContigsLen seqOffset = getSeqOffset(sa_info);
                    TSAPos chromlength = length(genome[getSeqNo(sa_info)]);
                    //Info make sure we dont DS search something going over the chromosom edge
                    //check left chromosom boundry && check right chromosom boundry
                    if(!(needleLeftPos <= seqOffset && chromlength - 1 >= seqOffset - needleLeftPos + length(needle) - 1))
                        continue;

//                     sa_info.i2 = sa_info.i2 - needleLeftPos;
                    setSeqOffset(sa_info, seqOffset - needleLeftPos);

                    //search remaining blocks
                    genomeSearch(ossContext, delegateDirect, needle, needleId, errors, s, blockIndex, genome, sa_info, blockStarts, blockEnds);
                }
            }
        }
        else
        {
            //modify blockend in case we are still inside a block
            if(needleLeftPos > blockStarts[0] && needleLeftPos < blockEnds[0])
                blockEnds[0] = needleLeftPos;

            for(TSAPos i = 0; i < iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1; ++i){
//                 if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
                if(checkSinglePos(bitvectors, brange, i)){
                    TSAValue sa_info = iter.revIter.index->sa[iter.revIter.vDesc.range.i1 + i];
                    TContigsLen seqOffset = getSeqOffset(sa_info);
                    TSAPos chromlength = length(genome[getSeqNo(sa_info)]);
                    //check left chromosom boundry && check right chromosom boundry
                    if(!(chromlength - 1 >= seqOffset + needleRightPos - 1 && seqOffset + needleRightPos - 1 >= length(needle) + 1))
                        continue;
                    //calculate correct starting position of the needle  on the forward index
//                     sa_info.i2 = chromlength - sa_info.i2 - needleRightPos + 1;
                    setSeqOffset(sa_info, chromlength - seqOffset - needleRightPos + 1);

                    //search remaining blocks
                    genomeSearch(ossContext, delegateDirect, needle, needleId , errors, s, blockIndex, genome, sa_info, blockStarts, blockEnds);
                }
            }
        }
    }
}

template <typename TIndex,
          typename TSALength,
          typename TDir>
inline void request_bitvector_interval(Iter<TIndex, VSTree<TopDown<> > > iter,
                                       uint32_t numberOfSequences,
                                       uint8_t needed_bitvector,
                                       Pair<uint8_t, Pair<TSALength, TSALength>> & brangeOutput,
                                       TDir const & )
{
    Pair<TSALength, TSALength> dirrange = (std::is_same<TDir, Rev>::value) ? range(iter.fwdIter) : range(iter.revIter);
//     uint32_t nseq = countSequences(*iter.fwdIter.index);
    dirrange.i1 = dirrange.i1 - numberOfSequences;
    dirrange.i2 = dirrange.i2 - numberOfSequences;

    brangeOutput.i1 = needed_bitvector;
    brangeOutput.i2 = dirrange;

}

template <typename TIndex,
          typename TBitvectorPair,
          size_t nbrBlocks,
          typename TSALength>
inline void get_bitvector_interval_inside(Iter<TIndex, VSTree<TopDown<> > > iter,
                                          std::vector<TBitvectorPair > & bitvectors,
                                          uint32_t numberOfSequences,
                                          OptimalSearch<nbrBlocks> const & s,
                                          uint8_t const blockIndex,
                                          Pair<uint8_t, Pair<TSALength, TSALength>> & brangeOutput,
                                          bool const goToRight2)
{
    Pair<TSALength, TSALength> dirrange = (goToRight2) ? range(iter.revIter) : range(iter.fwdIter);
    uint8_t needed_bitvector;
    uint8_t bitvsize = bitvectors.size();
    if (goToRight2)
        needed_bitvector = bitvsize - s.max[blockIndex - 1];
    else
        needed_bitvector = s.min[blockIndex - 1] - 1;
    //TODO use request_bitvector_interval

//     uint32_t nseq = countSequences(*iter.fwdIter.index);
    dirrange.i1 = dirrange.i1 - numberOfSequences;
    dirrange.i2 = dirrange.i2 - numberOfSequences;

    brangeOutput.i1 = needed_bitvector;
    brangeOutput.i2 = dirrange;
}


template <typename TIndex,
          typename TBitvectorPair,
          size_t nbrBlocks,
          typename TSALength,
          typename TDir>
inline void get_bitvector_interval(Iter<TIndex, VSTree<TopDown<> > > iter,
                                   std::vector<TBitvectorPair > & bitvectors,
                                   uint32_t numberOfSequences,
                                   OptimalSearch<nbrBlocks> const & s,
                                   uint8_t const blockIndex,
                                   Pair<uint8_t, Pair<TSALength, TSALength>> & brangeOutput,
                                   TDir const & )
{
    uint8_t needed_bitvector;
    if (std::is_same<TDir, Rev>::value)
        needed_bitvector = s.min[blockIndex] - 1;
    else
        needed_bitvector = bitvectors.size() - s.max[blockIndex];// + 1 - 1

    request_bitvector_interval(iter, numberOfSequences, needed_bitvector, brangeOutput, TDir());
}

template<typename TContex,
         typename TIndex,
         typename TBitvectorPair,
         typename TSALength,
         size_t nbrBlocks>
inline bool testUnidirectionalFilter(TContex & ossContext,
                                     Iter<TIndex, VSTree<TopDown<> > > iter,
                                     std::vector<TBitvectorPair > & bitvectors,
                                     Pair<uint8_t, Pair<TSALength, TSALength>> & brange,
                                     OptimalSearch<nbrBlocks> const & s,
                                     uint8_t const blockIndex,
                                     bool const goToRight2)
{
    // need bitinterval from inside the pattern to filter according to the mappability form
    //therefore i also need to acces the block before because of that block i got mappability of both sides
    Pair<uint8_t, Pair<TSALength, TSALength>> bit_interval;
    uint32_t numberOfSequences = ossContext.numberOfSequences;
    get_bitvector_interval_inside(iter, bitvectors, numberOfSequences, s, blockIndex, bit_interval, goToRight2);
    auto & b2 = getTVector(bitvectors, bit_interval);

    //squash interval
    TSALength startPos = bit_interval.i2.i1, endPos = bit_interval.i2.i2;
    TSALength startPos2 = startPos;
    TSALength endPos2 = endPos;

    while(b2[startPos] == 0 && startPos < endPos)
        ++startPos;

    while(b2[endPos - 1] == 0 && endPos > startPos)
        --endPos;

    if(startPos > endPos){
        std::cerr << "Error bit vector has only zeroes this should have been checked by checkinterval" << "\n";
        exit(0);
    }

    TSALength ivalSize = brange.i2.i2 - brange.i2.i1;
    TSALength count = 0;

    if(ossContext.normal.testflipdensity){
        // order of bits
        bool last = b2[startPos];
        TSALength pos = startPos;
        while(pos < endPos){
            if(b2[pos] != last){
                ++count;
                last = !last;
            }
            ++pos;
        }
    }
    // if next condition is true then brange will be modified!!!!
    // it will contain mappability of the bitvector anchored at the other side of already searched needle

    // only interested in changes inside the supinterval (startPos - endPos)
    // allowed flips per intervalSize
    if(!ossContext.normal.testflipdensity || static_cast<float>(ivalSize) * ossContext.normal.invflipdensity - 1 > static_cast<float>(count)){
        brange.i1 = bit_interval.i1;
        brange.i2.i1 = startPos;
        brange.i2.i2 = endPos;
        return true;
    }
    return false;
}

template<typename TContex,
         typename TBitvectorPair,
         typename TSALength,
         size_t nbrBlocks>
inline ReturnCode checkInterval(TContex & ossContext,
                                std::vector<TBitvectorPair > & bitvectors,
                                Pair<uint8_t, Pair<TSALength, TSALength>> & brange,
                                OptimalSearch<nbrBlocks> const & s,
                                uint8_t const blockIndex)
{
    auto & b = getTVector(bitvectors, brange);
    auto & rb = getTVSupport(bitvectors, brange);
    rb.set_vector(&b);

    TSALength ivalOne = rb(brange.i2.i2) - rb(brange.i2.i1);
    TSALength ivalSize = brange.i2.i2 - brange.i2.i1;

    if(ossContext.normal.nomappability && ivalOne == 0)
        return ReturnCode::NOMAPPABILITY;

    if(ossContext.normal.directsearch && ossContext.itvCondition(s, blockIndex, ivalOne))
        return ReturnCode::DIRECTSEARCH;

    if(ossContext.normal.compmappable && ivalOne == (brange.i2.i2 - brange.i2.i1)) //TODO maybe allow some zeroes
        return ReturnCode::COMPMAPPABLE;
/*
    //equal or more than half zeroes
    if(ossContext.normal.suspectunidirectional && s.startUniDir <= blockIndex && static_cast<float>(ivalOne) / static_cast<float>(ivalSize) <= ossContext.normal.filter_th)
        return ReturnCode::SUSPECTUNIDIRECTIONAL;*/

    return ReturnCode::MAPPABLE;
}

template <typename TSpec, typename TConfig,
          typename TDelegate, typename TDelegateD,
          typename TIndex,
          typename TNeedle,
          typename TBitvectorPair,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline ReturnCode checkCurrentMappability(OSSContext<TSpec, TConfig> & ossContext,
                                          TDelegate & delegate,
                                          TDelegateD & delegateDirect,
                                          Iter<TIndex, VSTree<TopDown<> > > iter,
                                          Pair <int8_t, int8_t> limOffsets,
                                          TNeedle const & needle,
                                          uint32_t needleId,
                                          std::vector<TBitvectorPair > & bitvectors,
                                          uint32_t const needleLeftPos,
                                          uint32_t const needleRightPos,
                                          uint8_t const errors,
                                          OptimalSearch<nbrBlocks> const & s,
                                          uint8_t const blockIndex,
                                          uint8_t const minErrorsLeftInBlock,
                                          TDir const & ,
                                          TDistanceTag const &)
{
    typedef typename TConfig::TContigsSum   TContigsSum; //SAVALUE
    Pair<uint8_t, Pair<TContigsSum, TContigsSum>> bit_interval;
    uint32_t numberOfSequences = ossContext.numberOfSequences;
    get_bitvector_interval(iter, bitvectors, numberOfSequences, s, blockIndex, bit_interval, TDir());

    ReturnCode rcode = checkInterval(ossContext, bitvectors, bit_interval, s, blockIndex);

    switch(rcode){
        case ReturnCode::NOMAPPABILITY:
            return ReturnCode::FINISHED;

        case ReturnCode::DIRECTSEARCH:
        {
            directSearch(ossContext, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, bit_interval, TDir(), TDistanceTag());
            return ReturnCode::FINISHED;
        }

        case ReturnCode::COMPMAPPABLE:
        {
            std::vector<TBitvectorPair > empty_bitvectors;
            _optimalSearchSchemeChildren(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, empty_bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
            return ReturnCode::FINISHED;
        }

        default:
            return ReturnCode::MAPPABLE;
    }
}


template <typename TSpec, typename TConfig,
          typename TDelegate, typename TDelegateD,
          typename TIndex,
          typename TNeedle,
          typename TBitvectorPair,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline ReturnCode checkMappability(OSSContext<TSpec, TConfig> & ossContext,
                                   TDelegate & delegate,
                                   TDelegateD & delegateDirect,
                                   Iter<TIndex, VSTree<TopDown<> > > iter,
                                   Pair <int8_t, int8_t> limOffsets,
                                   TNeedle const & needle,
                                   uint32_t needleId,
                                   std::vector<TBitvectorPair > & bitvectors,
                                   uint32_t const current_needleLeftPos,
                                   uint32_t const current_needleRightPos,
                                   uint8_t const errors,
                                   OptimalSearch<nbrBlocks> const & s,
                                   uint8_t const blockIndex,
                                   bool const lastEdit,
                                   TDir const & ,
                                   TDistanceTag const &)
{
    typedef typename TConfig::TContigsSum       TContigsSum; //SAVALUE
    Pair<uint8_t, Pair<TContigsSum, TContigsSum>> bit_interval;
    uint32_t numberOfSequences = ossContext.numberOfSequences;

    get_bitvector_interval(iter, bitvectors, numberOfSequences, s, blockIndex, bit_interval, TDir());

    ReturnCode rcode = checkInterval(ossContext, bitvectors, bit_interval, s, blockIndex);
    switch(rcode)
    {
        case ReturnCode::NOMAPPABILITY:
            return ReturnCode::FINISHED;

        case ReturnCode::DIRECTSEARCH:
        {
            //search directly in Genome
            directSearch(ossContext, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, TDir(), TDistanceTag());
            return ReturnCode::FINISHED;
        }

        case ReturnCode::COMPMAPPABLE:
        {
            std::vector<TBitvectorPair > empty_bitvectors;
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, empty_bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, lastEdit, TDir(), TDistanceTag());
            return ReturnCode::FINISHED;
        }
/*
        case ReturnCode::SUSPECTUNIDIRECTIONAL:
        {
            //test unidirectional changes iter range if true
            //TODO modfy functions for TDIR
            bool goToRight2 = std::is_same<TDir, Rev>::value;
            if(testUnidirectionalFilter(ossContext, iter, bitvectors, bit_interval, s, blockIndex, goToRight2)){
                //range on iter was changed in function before
                filter_interval(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, TDir(), TDistanceTag());
                return ReturnCode::FINISHED;
            }
        }*/
        default:
            return ReturnCode::MAPPABLE;
    }
}

template <typename TContex,
          typename TDelegate,
          typename TDelegateDirect,
          typename TIndex,
          typename TNeedle,
          typename TBitvectorPair,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeDeletion(TContex & ossContext,
                                         TDelegate & delegate,
                                         TDelegateDirect & delegateDirect,
                                         Iter<TIndex, VSTree<TopDown<> > > iter,
                                         Pair <int8_t, int8_t> limOffsets,
                                         TNeedle const & needle,
                                         uint32_t needleId,
                                         std::vector<TBitvectorPair > & bitvectors,
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         bool const lastEdit,
                                         TDir const & /**/)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    if (minErrorsLeftInBlock == 0)
    {
        uint8_t const blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
        bool const goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

        if (goToRight2)
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex2, lastEdit, Rev(), EditDistance());
        else
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex2, lastEdit, Fwd(), EditDistance());
    }

    //    bool goToRight = std::is_same<TDir, Rev>::value;
    bool not_at_end = std::is_same<TDir, Rev>::value && needleRightPos != length(needle) + 1 || !std::is_same<TDir, Rev>::value && needleLeftPos != 0/* || true*/;

    if(std::is_same<TDir, Rev>::value)
        --limOffsets.i2;
    else
        --limOffsets.i1;

    if (not_at_end && maxErrorsLeftInBlock > 0 && goDown(iter, TDir()))
    {
        do
        {
            _optimalSearchSchemeDeletion(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors + 1, s, blockIndex, true, TDir());
        } while (goRight(iter, TDir()));
    }
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TIndex,
          typename TNeedle,
          typename TBitvectorPair,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeChildren(TContex & ossContext,
                                         TDelegate & delegate,
                                         TDelegateD & delegateDirect,
                                         Iter<TIndex, VSTree<TopDown<> > > iter,
                                         Pair <int8_t, int8_t> limOffsets,
                                         TNeedle const & needle,
                                         uint32_t needleId,
                                         std::vector<TBitvectorPair > & bitvectors,
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         uint8_t const minErrorsLeftInBlock,
                                         TDir const & ,
                                         TDistanceTag const &)
{
    bool goToRight = std::is_same<TDir, Rev>::value;
    auto newlimOffsets = limOffsets;
    if(goToRight)
        --newlimOffsets.i2;
    else
        --newlimOffsets.i1;

    if (goDown(iter, TDir()))
    {
       	uint32_t charsLeft = s.blocklength[blockIndex] - (needleRightPos - needleLeftPos - 1);
        do
        {
            bool delta = !ordEqual(parentEdgeLabel(iter, TDir()), needle[goToRight ? needleRightPos - 1 : needleLeftPos - 1]);
            if (!std::is_same<TDistanceTag, EditDistance>::value && minErrorsLeftInBlock > 0 && charsLeft + delta < minErrorsLeftInBlock + 1u)
            {
                continue;
            }
            int32_t needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t needleRightPos2 = needleRightPos + goToRight;
            //finished Block
            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                if (std::is_same<TDistanceTag, EditDistance>::value)
                {
                    //use delta instead of false if no mismatches are allowed
                    _optimalSearchSchemeDeletion(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, false, TDir());
                }
                else
                {
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
                    if (goToRight2)
                        _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, false, Rev(), TDistanceTag());
                    else
                        _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, false, Fwd(), TDistanceTag());
                }
            }
            else
            {

                //if want to disable mismatches at the start and end (!delta || not_at_end) && use delta instead of false
                _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, false, TDir(), TDistanceTag());
            }

            //Deletion
            bool not_at_end = std::is_same<TDir, Rev>::value && needleRightPos2 != length(needle) + 1 || !std::is_same<TDir, Rev>::value && needleLeftPos2 != 0/* || true*/;
            if (std::is_same<TDistanceTag, EditDistance>::value && not_at_end){
                _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, newlimOffsets, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors + 1, s, blockIndex, true, TDir(), TDistanceTag());
            }
        } while (goRight(iter, TDir()));
    }
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TIndex,
          typename TNeedle,
          typename TBitvectorPair,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeExact(TContex & ossContext,
                                      TDelegate & delegate,
                                      TDelegateD & delegateDirect,
                                      Iter<TIndex, VSTree<TopDown<> > > iter,
                                      Pair <int8_t, int8_t> & limOffsets,
                                      TNeedle const & needle,
                                      uint32_t needleId,
                                      std::vector<TBitvectorPair > & bitvectors,
                                      uint32_t const needleLeftPos,
                                      uint32_t const needleRightPos,
                                      uint8_t const errors,
                                      OptimalSearch<nbrBlocks> const & s,
                                      uint8_t const blockIndex,
                                      TDir const &,
                                      TDistanceTag const &)
{
    // not in last block and next Block is larger then current block
    bool goToRight2 = (blockIndex < s.pi.size() - 1) ? s.pi[blockIndex + 1] > s.pi[blockIndex] : s.pi[blockIndex] > s.pi[blockIndex - 1];
    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
    if (std::is_same<TDir, Rev>::value)
    {
        //search take rest of the block and search it reverse
        uint32_t infixPosLeft = needleRightPos - 1;
        uint32_t infixPosRight = needleLeftPos + s.blocklength[blockIndex] - 1;

	while(infixPosLeft < infixPosRight + 1){
          //  if (!goDown(iter, infix(needle, infixPosLeft, infixPosRight + 1), TDir()))
            if (!goDown(iter, needle[infixPosLeft], TDir()))
                return;
	    ++infixPosLeft;
	}
        if (goToRight2)
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, false, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, false, Fwd(), TDistanceTag());
    }
    else
    {
        // has to be signed, otherwise we run into troubles when checking for -1 >= 0u
        int32_t infixPosLeft = needleRightPos - s.blocklength[blockIndex] - 1;
        int32_t infixPosRight = needleLeftPos - 1;

        while (infixPosRight >= infixPosLeft)
        {
            if (!goDown(iter, needle[infixPosRight], TDir())){
                return;
            }
            --infixPosRight;
        }

        if (goToRight2)
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, false, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, false, Fwd(), TDistanceTag());
    }
}

template <typename TSpec, typename TConfig,
          typename TDelegate,
          typename TIndex,
          typename TBitvectorPair>
inline void filteredDelegate(OSSContext<TSpec, TConfig> & ossContext,
                             TDelegate & delegate,
                             Iter<TIndex, VSTree<TopDown<> > > iter,
                             Pair <int8_t, int8_t> limOffsets,
                             uint32_t needleId,
                             std::vector<TBitvectorPair > & bitvectors,
                             uint8_t const errors)
{
    typedef typename TConfig::TContigsSum       TContigsSum; //SAVALUE
    uint32_t numberOfSequences = ossContext.numberOfSequences;
    Pair<uint8_t, Pair<TContigsSum, TContigsSum>> left_bit_interval;
    request_bitvector_interval(iter, numberOfSequences, 0, left_bit_interval, Rev());

    TContigsSum rangeStart = iter.fwdIter.vDesc.range.i1;
    TContigsSum rangeEnd = iter.fwdIter.vDesc.range.i2;
    TContigsSum lastStart = 0;
    for(TContigsSum i = 0; i < rangeEnd - rangeStart; ++i)
    {
        auto & b = getTVector(bitvectors, left_bit_interval);
        if(b[left_bit_interval.i2.i1 + i] == 0 )
        {
            if(i != lastStart){
                iter.fwdIter.vDesc.range.i1 = rangeStart + lastStart;
                iter.fwdIter.vDesc.range.i2 = rangeStart + i - 1;
                //TODO add direction
//                 delegate(ossContext, iter, limOffsets, needleId, errors, DIR);
            }
            lastStart = i + 1;
        }
    }
    if(lastStart < rangeEnd - rangeStart){
        iter.fwdIter.vDesc.range.i1 = rangeStart + lastStart;
        iter.fwdIter.vDesc.range.i2 = rangeStart + rangeEnd - rangeStart;
        //TODO add direction
//         delegate(ossContext, iter, limOffsets, needleId, errors);
    }
}

template <typename TSpec, typename TConfig,
          typename TDelegate, typename TDelegateD,
          typename TIndex,
          typename TNeedle,
          typename TBitvectorPair,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchScheme(OSSContext<TSpec, TConfig> & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<TIndex, VSTree<TopDown<> > > iter,
                                 Pair <int8_t, int8_t> limOffsets,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 std::vector<TBitvectorPair > & bitvectors,
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 bool const lastEdit,
                                 TDir const & ,
                                 TDistanceTag const &)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;
    bool const done = minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1;
    bool const atBlockEnd = (blockIndex > 0) ? needleRightPos - needleLeftPos - 1 == s.blocklength[blockIndex - 1] : false;        //is not true if we finished needle
    bool const checkMappa = !bitvectors.empty();

    // Done. (Last step)
    if (done)
    {
        if(!lastEdit){
            bool goToRightT = (blockIndex < s.pi.size() - 1) ? s.pi[blockIndex + 1] > s.pi[blockIndex] : s.pi[blockIndex] > s.pi[blockIndex - 1];
            delegate(ossContext, iter, limOffsets, needleId, errors, goToRightT);
        }
        return;
    }

    if(atBlockEnd && checkMappa){
        ReturnCode rcode = checkMappability(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, lastEdit, TDir(), TDistanceTag());
        if(rcode == ReturnCode::FINISHED)
            return;
    }

    // Exact search in current block.
    if (maxErrorsLeftInBlock == 0)
    {
        _optimalSearchSchemeExact(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir(), TDistanceTag());
    }
    else if(!checkMappa && ossContext.itvConditionComp(iter, needleLeftPos, needleRightPos, errors, s, blockIndex))
    {
        typedef typename TConfig::TContigsSum   TContigsSum;
        //give emtpy bitvector and bitvector range sine we will not check mappability
        Pair<uint8_t, Pair<TContigsSum, TContigsSum>> dummy_bit_interval;
         directSearch(ossContext, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, dummy_bit_interval, TDir(), TDistanceTag());
    }

    // Approximate search in current block.
    else
    {
        bool not_at_end = std::is_same<TDir, Rev>::value && needleRightPos != length(needle) || !std::is_same<TDir, Rev>::value && needleLeftPos != 1/* || true*/;
        // Insertion
        if (std::is_same<TDistanceTag, EditDistance>::value && not_at_end)
        {
            bool const goToRight = std::is_same<TDir, Rev>::value;
            int32_t const needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t const needleRightPos2 = needleRightPos + goToRight;
            auto newlimOffsets = limOffsets;
            if(goToRight)
                ++newlimOffsets.i2;
            else
                ++newlimOffsets.i1;

            //if we are at the end of block we need to add possible deletions because _optimalSearchScheme does not check it
            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                uint8_t const minErrorsLeftInBlock2 = (s.l[blockIndex] > (errors + 1)) ? (s.l[blockIndex] - (errors + 1)) : 0;
                if (minErrorsLeftInBlock2 == 0)
                {
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

                    if (goToRight2)
                        _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, newlimOffsets, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex2, true, Rev(), TDistanceTag());
                    else
                        _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, newlimOffsets, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex2, true, Fwd(), TDistanceTag());
                }
            }
            else
            {
                _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, newlimOffsets, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex, true, TDir(), TDistanceTag());
            }
        }

        //checkCurrentMappability (inside a Block)
        if(!atBlockEnd && checkMappa && ossContext.inBlockCheckMappabilityCondition(needleLeftPos, needleRightPos, s, blockIndex))
        {
            ReturnCode rcode = checkCurrentMappability(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
            if(rcode == ReturnCode::FINISHED)
                return;
        }

        _optimalSearchSchemeChildren(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
    }
}


template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TIndex,
          typename TBitvectorPair,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDistanceTag>
inline void _optimalSearchScheme(TContex & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<TIndex, VSTree<TopDown<> > > it,
                                 std::vector<TBitvectorPair > & bitvectors,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 OptimalSearch<nbrBlocks> const & s,
                                 TDistanceTag const &)
{
    bool initialDirection = s.pi[1] > s.pi[0];
    uint32_t max_e = ossContext.maxError; //unify overlap between all used schemes
    if(initialDirection)
        _optimalSearchScheme(ossContext, delegate, delegateDirect, it, Pair<int8_t, int8_t>(max_e, max_e), needle, needleId, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, false, Rev(), TDistanceTag());
    else
        _optimalSearchScheme(ossContext, delegate, delegateDirect, it, Pair<int8_t, int8_t>(max_e, max_e), needle, needleId, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, false, Fwd(), TDistanceTag());
}

template <typename TSpec, typename TConfig,
          typename TDelegate,
          typename TDelegateD,
          typename TIndex,// typename TIndexSpec,
          typename TBitvectorPair,
          typename TNeedle,
          size_t nbrBlocks, size_t N,
          typename TDistanceTag>
inline void _optimalSearchScheme(OSSContext<TSpec, TConfig> & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<TIndex, VSTree<TopDown<> > > it,
                                 std::vector<TBitvectorPair > & bitvectors,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss,
                                 TDistanceTag const &)
{
    typedef MapperTraits<TSpec, TConfig>                     TTraits;
    typedef typename TTraits::TContigsPos                    TContigsPos;
    typedef typename TConfig::TContigsSum                    TContigsSum;
    typedef typename TTraits::TContigSeqs                    TContigSeqs;

    TContigSeqs const & genome = ossContext.contigSeqs;
    if(!ossContext.noSAfilter){
        std::vector<seqan::SARange<TContigsSum> > ranges;

        auto delegateRange = [&ranges](OSSContext<TSpec, TConfig> & ossContext, auto const & iter, auto & limOffsets, auto const needleId, uint8_t const errors, bool const goRight)
        {
            //only required for stats
            ossContext.delegateOcc += iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1;
        //         std::cout << "Using new delegate: \n";
            uint32_t readId = getReadId(ossContext.readSeqs, needleId);
            SARange<TContigsSum> range;
            range.errors = errors;
            range.range = Pair<TContigsSum, TContigsSum> (iter.fwdIter.vDesc.range.i1, iter.fwdIter.vDesc.range.i2);
            range.revRange = Pair<TContigsSum, TContigsSum> (iter.revIter.vDesc.range.i1, iter.revIter.vDesc.range.i2);
            range.parentRange = Pair<TContigsSum, TContigsSum> (iter.fwdIter._parentDesc.range.i1, iter.fwdIter._parentDesc.range.i2);
            range.parentRevRange = Pair<TContigsSum, TContigsSum> (iter.revIter._parentDesc.range.i1, iter.revIter._parentDesc.range.i2);
            range.repLength = repLength(iter);
            range.limOffsets = limOffsets;
            range.lastChar = (goRight) ? iter.revIter.vDesc.lastChar : iter.fwdIter.vDesc.lastChar;
            range.goRight = goRight;
            ranges.push_back(range);

            setMapped(ossContext.ctx, readId);
            setMinErrors(ossContext.ctx, readId, errors);

        };

        for (auto & s : ss){
            _optimalSearchScheme(ossContext, delegateRange, delegateDirect, it, bitvectors, needle, needleId, s, TDistanceTag());
        }

        if(!ranges.empty()){
            ++ossContext.filteredOccsOfRead;
            locateSARanges(ossContext, genome, delegate, needleId, it, ranges);
        }
    }
    else
    {
//         std::cout << "no SA filter\n";
        auto delegateFMTree = [&delegate](OSSContext<TSpec, TConfig> & ossContext, auto const & iter, auto & limOffsets, auto const needleId, uint8_t const errors, bool const goRight)
        {
            TContigSeqs const & genome = ossContext.contigSeqs;

            //only required for stats
//             ossContext.delegateOcc += iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1;

            SARange<TContigsSum> saRange;
            saRange.repLength = repLength(iter);
            saRange.limOffsets = limOffsets;
            //TODO check errors
            saRange.errors = errors;

            uint32_t readId = getReadId(ossContext.readSeqs, needleId);
            setMapped(ossContext.ctx, readId);
            setMinErrors(ossContext.ctx, readId, errors);


            if(ossContext.fmTreeThreshold < iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1 && ossContext.samplingRate > 1)
            {
                ossContext.fmtreeLocates += iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1;
                auto uniIter = (ossContext.earlyLeaf && goRight) ? iter.revIter : iter.fwdIter;
                uint32_t counter = 0;

                //early leaf node calculation (corresponds to the last level of the fm_tree)
                if(ossContext.earlyLeaf){
                    TContigsSum ssp = getRank(uniIter.index->sa.sparseString.indicators, uniIter._parentDesc.range.i1 - 1);
                    TContigsSum sep = getRank(uniIter.index->sa.sparseString.indicators, uniIter._parentDesc.range.i2 - 1);
                    auto lastChar = uniIter.vDesc.lastChar;

                    if (goRight){
//                         std::cout << "\nrevIndex\n";
                        for(TContigsSum i = ssp; i < sep; ++i){
//                        for(TContigsSum i = iter.revIter._parentDesc.range.i1; i < iter.revIter._parentDesc.range.i2; ++i){
//                            if(getValue(iter.revIter.index->sa.sparseString.indicators, i)){
                             TContigsPos pos = getValue(uniIter.index->sa.sparseString.values, i);
//                             TContigsPos pos = iter.revIter.index->sa[i];
                             setSeqOffset(pos, getSeqOffset(pos) - 1);
                        //     std::cout << "Rev Pos: " << pos << "\t" << i << "\n";
                             //Calculate last position of backtracking in reverse Index on the forward Text;
                            setSeqOffset(pos, ossContext.sequenceLengths[getSeqNo(pos)] - getSeqOffset(pos) - 1);
                      //      std::cout << "End Pos: " << pos << "\n";
                            //Calculate start position

                         //   std::cout << "lastChar " << lastChar << "\n";
                            if (genome[getSeqNo(pos)][getSeqOffset(pos)] == lastChar){
                                //searched one less char the we use range from parent
                                setSeqOffset(pos, getSeqOffset(pos) - (saRange.repLength - 1));
        //                        std::cout << "Genome: " << infix(genome, pos, posAdd(pos, 140)) << "\n";
                                delegate(ossContext, needleId, saRange, pos);
                            }
                            else
                            {
                                setSeqOffset(pos, getSeqOffset(pos) - (saRange.repLength - 1));
                            }
                        }

                    }
                    else
                    {
                        for(TContigsSum i = ssp; i < sep; ++i){
                            TContigsPos pos = posAdd(getValue(uniIter.index->sa.sparseString.values, i), (-1));
                            if (genome[getSeqNo(pos)][getSeqOffset(pos)] == lastChar){
                                //Report position
                                delegate(ossContext, needleId, saRange, pos);
                            }
                        }
                    }
                }

                bool goRight2 = ossContext.earlyLeaf && goRight;
                fm_tree(ossContext, delegate, needleId, saRange, uniIter, goRight2, 0, static_cast<uint8_t>(ossContext.earlyLeaf), counter);
            }
            else
            {
//                 ossContext.defaultLocates += saRange.range.i2 - saRange.range.i1;
                for (TContigsSum i = iter.fwdIter.vDesc.range.i1; i < iter.fwdIter.vDesc.range.i2; ++i)
                {
                    TContigsPos pos = iter.fwdIter.index->sa[i];
                    delegate(ossContext, needleId, saRange, pos);
                }
            }
        };

        for (auto & s : ss)
            _optimalSearchScheme(ossContext, delegateFMTree, delegateDirect, it, bitvectors, needle, needleId, s, TDistanceTag());
    }
}

/*
template <typename TSpec, typename TConfig,
          typename TDelegate,
          typename TDelegateD,
          typename TIndex,// typename TIndexSpec,
          typename TBitvectorPair,
          typename TNeedle,
          size_t nbrBlocks, size_t N,
          typename TDistanceTag>
inline void _optimalSearchScheme(OSSContext<TSpec, TConfig> & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<TIndex, VSTree<TopDown<> > > it,
                                 std::vector<TBitvectorPair > & bitvectors,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss,
                                 TDistanceTag const &)
{
    for (auto & s : ss)
        _optimalSearchScheme(ossContext, delegate, delegateDirect, it, bitvectors, needle, needleId, s, TDistanceTag());
}*/

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TIndex,
//           typename TChar, typename TStringSpec,
          typename TBitvectorPair,
          typename TNeedle,
          size_t nbrBlocks, size_t N,
          typename TDistanceTag>
inline void
find(TContex & ossContext,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     TIndex & index,
     std::vector<TBitvectorPair > & bitvectors,
     TNeedle/*String<TChar, TStringSpec>*/ const & needle,
     uint32_t needleId,
     std::array<OptimalSearch<nbrBlocks>, N> const & ss,
     TDistanceTag const & )
{
    auto scheme = ss;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(needle));
    _optimalSearchSchemeComputeChronBlocklength(scheme);
//     Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);
    Iter<TIndex, VSTree<TopDown<> > > it(index);
    _optimalSearchScheme(ossContext, delegate, delegateDirect, it, bitvectors, needle, needleId, scheme, TDistanceTag());

}


template <size_t minErrors, size_t maxErrors,
          typename TSpec, typename TConfig,
          typename TDelegate, typename TDelegateD,
          typename TIndex,
          typename TBitvectorPair,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void
find(OSSContext<TSpec, TConfig> & ossContext,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     TIndex & index,
     std::vector<TBitvectorPair > & bitvectors,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & )
{
    typedef MapperTraits<TSpec, TConfig>                                                TTraits;
    typedef typename Iterator<StringSet<TNeedle, TStringSetSpec> const, Rooted>::Type   TReadIt;
    typedef typename Reference<TReadIt>::Type                                           TReadRef;
    typedef typename TTraits::TReadSeqs                                                 TReadSeqs;
    typedef typename Size<TReadSeqs>::Type                                              TReadId;

    //TODO use readLength to only calculate scheme for bitvectors use real read length to search for it
    auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
    calcConstParameters(scheme);

    std::cout << "Scheme: " << minErrors << "\t" << maxErrors << ">\n";
    std::vector<TBitvectorPair * > lbitvectors;
    //load Bitvectors needed for scheme (Blocklength and chronblockLengths have to be calculated therefore I need to assume needle length)
    //use either specified maxumum readLength given as input or length of first needle
    linkBitvectors(ossContext, scheme, bitvectors, lbitvectors);

    // Iterate over all reads.
    iterate(needles, [&](TReadIt const & readIt)
    {
        bool skip = false;
        TReadRef it = value(readIt);
        TReadId readId = getReadId(ossContext.readSeqs, position(readIt));
//             std::cout << "ReadId: " << readId << "\n";

        if(isMapped(ossContext.ctx, readId)){
//                 std::cout << "MinErrors: " << (int)getMinErrors(ossContext.ctx, readId) << "\n";
            if(static_cast<uint8_t>(getMinErrors(ossContext.ctx, readId)) + ossContext.strata < minErrors){
//                     std::cout << "Skip" << "\n";
                skip = true;
            }
        }

        if(!skip){
//             std::cout << "Search\n";
            find(ossContext, delegate, delegateDirect, index, bitvectors, it, position(readIt), scheme, TDistanceTag());
        }
    }, Rooted(), typename TTraits::TThreading());
}


// Index<Void, BidirectionalIndex<TIndexSpec> > & index, ??
template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TBitvectorPair,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 std::vector<TBitvectorPair > & bitvectors, // cant be const since TVSupport.set_vector(&TVector)
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 TDistanceTag const & )
{
    switch (maxErrors)
    {
        case 0: find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                break;
        case 1: find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                break;
        case 2: find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                break;
        case 3: find<0, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                break;
        case 4: find<0, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                break;
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}

//no strata (needed for one Scheme Best X)
template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 TDistanceTag const & )
{
    std::vector<std::pair<TBitvector, TSupport>> empty_bitvectors;
    find(minErrors, maxErrors, ossContext, delegate, delegateDirect, index, empty_bitvectors, needles, TDistanceTag());
}


// for find2_index_approx.h find function
template <typename TDelegate,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
     const int maxErrors,
     TDelegate & delegate,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & )
{
    switch (maxErrors)
    {
        case 0: find<0, 0>(delegate, index, needles, TDistanceTag());
                break;
        case 1: find<0, 1>(delegate, index, needles, TDistanceTag());
                break;
        case 2: find<0, 2>(delegate, index, needles, TDistanceTag());
                break;
        case 3: find<0, 3>(delegate, index, needles, TDistanceTag());
                break;
        case 4: find<0, 4>(delegate, index, needles, TDistanceTag());
                break;
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}

}



#endif
