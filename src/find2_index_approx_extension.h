#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_

#include <iostream>
#include <sdsl/bit_vectors.hpp>
// #include "find2_index_approx_unidirectional.h"





namespace seqan{

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
    delegateDirect(ossContext, sa_info, posAdd(sa_info, length(needle)), errors, needleId);
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
/*
//skip mismatches in the start and end adjust sa_value and occlength accordingly
template <typename TNeedle, typename TSAValue, typename TInfix>
inline void trimHit(TNeedle const & needle, TInfix const & infix, TSAValue & sa_info, uint32_t & occlength, uint8_t errors)
{
    // trimming does not work everytime ex: TCCCCCCCCGCTA
    //                                      CCCCCCCCGCTA -> only after comparing C to G it is clear we need and insertion
    //                                                      in the beginning
    std::cout << "Trimming: \n";
    std::cout << "Needle" << sa_info << "\t" << occlength << "\n" << needle << "\n";
    std::cout << infix << "\n";
    uint32_t k = 0;
    while(needle[k] != infix[k] && k < errors + 2)
    {
        ++k;
    }

    //check if edit operation was required
    if(k == errors + 1)
        k = 0;

    uint32_t len = length(needle);
    uint32_t l = 1;

    uint32_t left_errors = errors + 2 - k;
    while(needle[len - l] != infix[occlength - l] && l < errors + 2 - k)
    {
        ++l;
    }
    // check if edit operation was required
    if(l == left_errors - 1)
        l = 1;


    auto seqOffset = getSeqOffset(sa_info);
    setSeqOffset(sa_info, seqOffset + k);
    occlength = occlength - k - l + 1;

    std::cout << "trimmed:" << sa_info << "\t" << occlength << "\n" << infix << "\n" << k << "\t" << l - 1 << "\n";
}*/


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

    TContigSeqs const & genome = ossContext.contigSeqs;

//     auto const & genome = indexText(*iter.fwdIter.index);

//     bool test = needle;
    ModifiedString<TNeedle, ModReverse > revneedle(needle);

    if (std::is_same<TDistanceTag, EditDistance>::value){
        //TODO put this into a function
        //TODO if we are only interested in the best hit call return after delegate calls
        uint32_t needleL = length(needle);
        uint32_t max_e = s.u[s.u.size() - 1];
        uint8_t intIns = 0;
        uint8_t intDel = 0;
        //calculate net sum of internal Insertions - Deletions
/*
        if(repLength(iter) < needleRightPos - needleLeftPos - 1)
            intIns = needleRightPos - needleLeftPos - 1 - repLength(iter);
        else
            intDel = repLength(iter) - (needleRightPos - needleLeftPos - 1);*/

//         uint8_t overlap_l = max_e;
//         uint8_t overlap_r = max_e;

//         std::cout << "NLP: " << needleLeftPos << "\tNRP: " << needleRightPos << "\trepL: " << (int)repLength(iter) << "\trange: " <<  (int)needleRightPos - needleLeftPos - 1 << "\n";
//         std::cout << "Off: " << (int)limOffsets.i1 - max_e << ", " << (int)limOffsets.i2 - max_e << "\n";

        uint8_t overlap_l = limOffsets.i1;
        uint8_t overlap_r = limOffsets.i2;

/*
        uint8_t overlap_l;
        uint8_t overlap_r;
        if(s.pi[0] == 1){overlap_l = 0;}
        else if(needleLeftPos == 0){overlap_l = intDel;}
        else{overlap_l = max_e;}

        if(s.pi[s.pi.back()] == s.pi.size()){overlap_r = 0;}
        else if(needleRightPos == needleL + 1){overlap_r = intDel;}
        else{overlap_r = max_e;}*/

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
                    setSeqOffset(sa_info, seqOffset - needleLeftPos); //maybe -1
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
                delegateDirect(ossContext, posAdd(sa_info, -overlap_l), posAdd(sa_info, needleL + overlap_r), needleId, 127);
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

//     ivalOne < (s.pi.size() - blockIndex - 1 + ossContext.normal.directsearchblockoffset) * ossContext.normal.directsearch_th
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

        if (!goDown(iter, infix(needle, infixPosLeft, infixPosRight + 1), TDir()))
            return;

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
          typename TNeedle,
          typename TBitvectorPair>
inline void filteredDelegate(OSSContext<TSpec, TConfig> & ossContext,
                             TDelegate & delegate,
                             Iter<TIndex, VSTree<TopDown<> > > iter,
                             Pair <int8_t, int8_t> limOffsets,
                             TNeedle const & needle,
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
                delegate(ossContext, iter, limOffsets, needle, needleId, errors, false);
            }
            lastStart = i + 1;
        }
    }
    if(lastStart < rangeEnd - rangeStart){
        iter.fwdIter.vDesc.range.i1 = rangeStart + lastStart;
        iter.fwdIter.vDesc.range.i2 = rangeStart + rangeEnd - rangeStart;
        delegate(ossContext, iter, limOffsets, needle, needleId, errors, false);
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
        auto newlimOffsets = limOffsets;
//         std::cout << "Done" << "\n";
        //last input only matters for unidirectional searches (has to be false in this case)
        if(!lastEdit/*true*/){
            if(checkMappa){
                filteredDelegate(ossContext, delegate, iter, limOffsets, needle, needleId, bitvectors, errors);
            }
            else
            {
//                 std::cout << "Calling delegate" << "\n";
                delegate(ossContext, iter, limOffsets, needle, needleId, errors, false);
            }
        }
        return;
    }
/*
    if(atBlockEnd && checkMappa){
        ReturnCode rcode = checkMappability(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, lastEdit, TDir(), TDistanceTag());
        if(rcode == ReturnCode::FINISHED)
            return;
    }*/

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
/*
        //checkCurrentMappability (inside a Block)
        if(!atBlockEnd && checkMappa && ossContext.inBlockCheckMappabilityCondition(needleLeftPos, needleRightPos, s, blockIndex))
        {
            ReturnCode rcode = checkCurrentMappability(ossContext, delegate, delegateDirect, iter, limOffsets, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
            if(rcode == ReturnCode::FINISHED)
                return;
        }*/

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
    uint32_t max_e = s.u[s.u.size() - 1];
    if(initialDirection)
        _optimalSearchScheme(ossContext, delegate, delegateDirect, it, Pair<int8_t, int8_t>(max_e, max_e), needle, needleId, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, false, Rev(), TDistanceTag());
    else
        _optimalSearchScheme(ossContext, delegate, delegateDirect, it, Pair<int8_t, int8_t>(max_e, max_e), needle, needleId, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, false, Fwd(), TDistanceTag());
}

template <typename TSpec, typename TConfig,
          typename TDelegate, typename TDelegateD,
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
    typedef MapperTraits<TSpec, TConfig>                      TTraits;
    typedef typename TTraits::TBiIter                         TIter;
    typedef typename TTraits::TSparseIter                     TSparseIter;

    for (auto & s : ss)
        _optimalSearchScheme(ossContext, delegate, delegateDirect, it, bitvectors, needle, needleId, s, TDistanceTag());
}


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
//     uint32_t k = 0;
    iterate(needles, [&](TReadIt const & readIt)
    {
        bool skip = false;
        TReadRef it = value(readIt);
        if(ossContext.bestXMapper){
            TReadId readId = getReadId(ossContext.readSeqs, position(readIt));
//             std::cout << "ReadId: " << readId << "\n";

            //advance condition for multiple SearchSchemes
            /*
            bool m = isMapped(ossContext.ctx, readId);
            uint8_t minE = getMinErrors(ossContext.ctx, readId);
            bool search = (maxErrors == ossContext.maxError ||  maxErrors - minErrors == ossContext.strata) && !m ||
                             m && std::min(minE + ossContext.strata, ossContext.maxError) == maxErrors;*/
            if(isMapped(ossContext.ctx, readId)){
//                 std::cout << "MinErrors: " << (int)getMinErrors(ossContext.ctx, readId) << "\n";
                if(static_cast<uint8_t>(getMinErrors(ossContext.ctx, readId)) + ossContext.strata < minErrors){
//                     std::cout << "Skip" << "\n";
                    skip = true;
                }
            }
        }

        if(!skip){
//             std::cout << "Search\n";
            find(ossContext, delegate, delegateDirect, index, bitvectors, it, position(readIt), scheme, TDistanceTag());
        }
//         k++;
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
