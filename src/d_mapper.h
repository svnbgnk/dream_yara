// ==========================================================================
//                                 d_mapper.h
// ==========================================================================
// Copyright (c) 2017-2022, Temesgen H. Dadi, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Temesgen H. Dadi or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL TEMESGEN H. DADI OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
// ==========================================================================


#ifndef APP_YARA_DIS_MAPPER_H_
#define APP_YARA_DIS_MAPPER_H_

#include "oss_context.h"

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class DisOptions
// --------------------------------------------------------------------------

struct DisOptions : public Options
{
public:
    CharString              IndicesDirectory;
    CharString              filterFile;
    CharString              superOutputFile;

    CharString              MappabilityDirectory;
    bool                    ossOff = false;
    bool                    noITV = false;
    bool                    noSAfilter = false;
    bool                    noDelayITV = false;
    bool                    determineExactSecondaryPos = true;
    bool                    noMappability = false;
    bool                    earlyLeaf = false;
    bool                    compare = false;
    double                  hammingDpercentage = 0.0;
    uint32_t                threshold = 11;
    uint32_t                itvOccThreshold = 10;
    uint32_t                fmTreeThreshold = 1000;
    uint32_t                fmTreeBreak = 10;
    uint32_t                startBin = 0;
    uint32_t                readLength = 0;

    double                  loadFilter      = 0.0;
    double                  filterReads     = 0.0;
    double                  copyReads       = 0.0;
    double                  copyAlignments  = 0.0;
    double                  moveCigars      = 0.0;

    bool                    skipSamHeader = false;

    uint32_t                kmerSize = 20;
    uint32_t                numberOfBins = 0;

    uint32_t                currentBinNo = 0;
    uint64_t                filteredReads = 0;
    std::vector<uint32_t>   contigOffsets;

    std::vector<std::vector<uint32_t>>          origReadIdMap;
    std::map<uint32_t, String<CigarElement<>>>  collectedCigars;

    FilterType      filterType = BLOOM;

    std::vector<std::string> filterTypeList = {"bloom", "kmer_direct", "none"};

    uint32_t getContigOffsets()
    {
        return contigOffsets[currentBinNo];
    }

    uint16_t getThreshold(uint16_t readLen)
    {
        uint16_t maxError = errorRate * readLen;

        // same as readLen - kmerSize + 1 - (maxError * kmerSize);
        if (kmerSize * (1 + maxError) > readLen)
            return 0;

        return readLen - kmerSize * (1 + maxError) + 1;
    }

};

template <typename TSpec, typename TConfig>
struct DelegateDirect
{
    typedef MapperTraits<TSpec, TConfig>                   Traits;
    typedef typename Traits::TMatch                        TMatch;
    typedef typename Traits::TMatchesAppender              TMatches;
    typedef typename Traits::TContigsPos                  TContigsPos;

    typedef typename Traits::TReadSeqs                     TReadSeqs;
    typedef typename Size<TReadSeqs>::Type                  TReadId;

    TMatches &          matches;

    DelegateDirect(TMatches & matches) :
        matches(matches)
    {}

    template <typename TNeedleId, typename TMatchErrors>
    void operator() (OSSContext<TSpec, TConfig> & ossContext, TContigsPos const & start, TContigsPos const & end, TNeedleId const needleId, TMatchErrors errors)
    {
        TReadId readId = getReadId(ossContext.readSeqs, needleId);
        TMatch hit;
        setContigPosition(hit, start, end);
//         SEQAN_ASSERT_LEQ(5, getSeqOffset(end) - getSeqOffset(start));
        hit.errors = errors;
        setReadId(hit, ossContext.readSeqs, needleId); // needleId is used to determine if read is reverse complement

        appendValue(matches, hit, Generous(), typename Traits::TThreading());

        if(errors <= ossContext.maxError){
            ++ossContext.itvOccs;
            setMapped(ossContext.ctx, readId);
            setMinErrors(ossContext.ctx, readId, errors);
        }
        else
        {
            ++ossContext.itvJobs;
        }
    }
};


template <typename TSpec, typename TConfig>
struct Delegate
{
    typedef MapperTraits<TSpec, TConfig>                    Traits;
    typedef typename Traits::TMatch                        TMatch;
    typedef typename Traits::TMatchesAppender              TMatches;
    typedef typename Traits::TContigsPos                   TContigsPos;

//     typedef typename TConfig::TContigsSum                       TContigsSum;


    typedef typename Traits::TReadSeqs                     TReadSeqs;
    typedef typename Size<TReadSeqs>::Type                  TReadId;


    TMatches &          matches;
    bool const          noOverlap;
    uint8_t const       maxError;

    Delegate(TMatches & matches,
             bool noOverlap,
             uint8_t maxError) :
        matches(matches),
        noOverlap(noOverlap),
        maxError(maxError)
    {}

    template <typename TNeedleId, typename TSARange, typename TPos>
    void operator() (OSSContext<TSpec, TConfig> & ossContext, TNeedleId const needleId, TSARange const & rangeInfo, TPos & pos)
    {
        //TODO move into function calls
        int16_t overlap_l = maxError;
        int16_t overlap_r = maxError;
        uint32_t occLength = rangeInfo.repLength;
        TMatch hit;
        if(noOverlap){
            setContigPosition(hit, pos, posAdd(pos, occLength));
        }
        else
        {
            overlap_l = (overlap_l <=  getSeqOffset(pos)) ? overlap_l : 0;
            setContigPosition(hit, posAdd(pos, -overlap_l), posAdd(pos, occLength + overlap_r));
        }

        hit.errors = rangeInfo.errors;
        setReadId(hit, ossContext.readSeqs, needleId); // needleId is used to determine if read is reverse complement

/*
        if (!noOverlap && (getMember(hit, ContigEnd()) - getMember(hit, ContigBegin())) <= ossContext.readLength)
        {
            std::cout << "After adding overlap Match\n";
            write(std::cout, hit);
            std::cout << "maxError: " << maxError << "\n";
            std::cout << "pos " << pos << "\t" << occLength << "\n";
        }*/


        appendValue(matches, hit, Generous(), typename Traits::TThreading()); //TODO does this make any sense (always single occ)

        if(!ossContext.noSAfilter && getSeqOffset(pos) == 0)
        {
//             std::cout << "append additional Match\n";
            TMatch hit2 = hit;
            hit2.errors = 127;
            setContigPosition(hit2, pos, posAdd(pos, occLength));
            appendValue(matches, hit2, Generous(), typename Traits::TThreading());
        }

        //read Context is done as soon as a range is reported
//             TReadId readId = getReadId(ossContext.readSeqs, needleId);
//             setMapped(ossContext.ctx, readId);
//             setMinErrors(ossContext.ctx, readId, errors);
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function appendStats()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void appendStats(Mapper<TSpec, TMainConfig> & mainMapper, Mapper<TSpec, TConfig> & childMapper)
{
    mainMapper.stats.loadContigs    += childMapper.stats.loadContigs;
    mainMapper.stats.loadReads      += childMapper.stats.loadReads;
    mainMapper.stats.collectSeeds   += childMapper.stats.collectSeeds;
    mainMapper.stats.findSeeds      += childMapper.stats.findSeeds;
    mainMapper.stats.classifyReads  += childMapper.stats.classifyReads;
    mainMapper.stats.rankSeeds      += childMapper.stats.rankSeeds;
    mainMapper.stats.extendHits     += childMapper.stats.extendHits;
    mainMapper.stats.loadingBitvectors += childMapper.stats.loadingBitvectors;
    mainMapper.stats.optimumSearch +=  childMapper.stats.optimumSearch;
    mainMapper.stats.inTextVerification +=  childMapper.stats.inTextVerification;
    mainMapper.stats.sortMatches    += childMapper.stats.sortMatches;
    mainMapper.stats.compactMatches += childMapper.stats.compactMatches;
    mainMapper.stats.selectPairs    += childMapper.stats.selectPairs;
    mainMapper.stats.verifyMatches  += childMapper.stats.verifyMatches;
    mainMapper.stats.alignMatches   += childMapper.stats.alignMatches;
    mainMapper.stats.writeMatches   += childMapper.stats.writeMatches;
    mainMapper.stats.rescuedReads   += childMapper.stats.rescuedReads;
}
/*

// ----------------------------------------------------------------------------
// Function copyMatches()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void copyMatches(Mapper<TSpec, TMainConfig> & mainMapper, Mapper<TSpec, TConfig> & childMapper, DisOptions & disOptions)
{
    start(mainMapper.timer);
    typedef typename MapperTraits<TSpec, TMainConfig>::TMatch             TMatch;
    typedef typename MapperTraits<TSpec, TMainConfig>::TThreading         TThreading;
    typedef typename MapperTraits<TSpec, TMainConfig>::TMatchesAppender   TMatchesAppender;

    TMatchesAppender appender(mainMapper.matchesByCoord);

    uint32_t matchCount = length(childMapper.matchesByCoord);
    for (uint32_t i = 0; i < matchCount; ++i)
    {
        TMatch currentMatch;
        uint32_t readId         = childMapper.matchesByCoord[i].readId;
        uint32_t origReadId     = readId;
        if(disOptions.filterType != NONE)
            origReadId = disOptions.origReadIdMap[disOptions.currentBinNo][readId];

        currentMatch.readId        = origReadId;
        currentMatch.contigId      = childMapper.matchesByCoord[i].contigId + disOptions.getContigOffsets();
        currentMatch.isRev         = childMapper.matchesByCoord[i].isRev;
        currentMatch.contigBegin   = childMapper.matchesByCoord[i].contigBegin;
        currentMatch.contigEnd     = childMapper.matchesByCoord[i].contigEnd;
        currentMatch.errors        = childMapper.matchesByCoord[i].errors;
        appendValue(appender, currentMatch, Generous(), TThreading());
    }
    stop(mainMapper.timer);
    disOptions.copyAlignments += getValue(mainMapper.timer);
}*/

template<typename TGMatch, typename TMatch>
inline void copyMatch(TGMatch & globalMatch, TMatch const & match, DisOptions & disOptions)
{
    uint32_t readId = match.readId;
    uint32_t origReadId = readId;
    if(disOptions.filterType != NONE)
        origReadId = disOptions.origReadIdMap[disOptions.currentBinNo][readId];
    globalMatch.readId = origReadId;
    globalMatch.contigId      = match.contigId + disOptions.getContigOffsets();
    globalMatch.isRev         = match.isRev;
    globalMatch.contigBegin   = match.contigBegin;
    globalMatch.contigEnd     = match.contigEnd;
    globalMatch.errors        = match.errors;
}

// ----------------------------------------------------------------------------
// Function copyMatches()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void copyMatches(Mapper<TSpec, TMainConfig> & mainMapper, Mapper<TSpec, TConfig> & childMapper, DisOptions & disOptions)
{
     start(mainMapper.timer);
    typedef typename MapperTraits<TSpec, TMainConfig>::TMatch             TMatch;
    typedef typename MapperTraits<TSpec, TMainConfig>::TMatches           TMatches;
//     typedef StringSet<TMatches, Segment<TMatches> >                       TMatchesSet;
//     typedef String<TMatch, Segment<TMatches> >                            TMatchesSegment; //still wrong?
    typedef typename Iterator<TMatches, Standard>::Type                   TMatchesIterator;
    typedef typename MapperTraits<TSpec, TMainConfig>::TThreading         TThreading;
    typedef typename MapperTraits<TSpec, TMainConfig>::TMatchesAppender   TMatchesAppender;

    TMatchesAppender appender(mainMapper.matchesByCoord);


    uint32_t matchCount = length(childMapper.matchesSetByCoord);
    for (uint32_t i = 0; i < matchCount; ++i)
    {
        auto const & matches = childMapper.matchesSetByCoord[i];
        auto matchIt = begin(matches, Standard());
        auto matchEnd = end(matches, Standard());
        while(matchIt != matchEnd){
             TMatch currentMatch;
             if(isValid(*matchIt)){
                copyMatch(currentMatch, *matchIt, disOptions);
                appendValue(appender, currentMatch, Generous(), TThreading());
             }
             ++matchIt;
        }
    }
    stop(mainMapper.timer);
    disOptions.copyAlignments += getValue(mainMapper.timer);
}

// ----------------------------------------------------------------------------
// Function display Cigars()
// ----------------------------------------------------------------------------
inline CharString cigars2String(String<CigarElement<> > const & cigStr)
{
    CharString cs;
    for (CigarElement<> const & cigElement : cigStr)
    {
        appendNumber(cs, cigElement.count);
        appendValue(cs, cigElement.operation);
    }
    return cs;
}

inline StringSet<CharString> cigarsSet2String(StringSet<String<CigarElement<> >, Segment< String<CigarElement<> >> > & cigStrSet)
{
    StringSet<CharString> css;
    for (String<CigarElement<>> const & cigStr : cigStrSet)
    {
        appendValue(css, cigars2String(cigStr));
    }
    return css;
}


// ----------------------------------------------------------------------------
// Function copyCigars()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void copyCigars(Mapper<TSpec, TMainConfig> & mainMapper, Mapper<TSpec, TConfig> & childMapper, DisOptions & disOptions)
{
    start(mainMapper.timer);
    typedef typename MapperTraits<TSpec, TConfig>::TMatch             TMatch;
    uint32_t matchCount = length(childMapper.primaryMatches);
    for (uint32_t i = 0; i < matchCount; ++i)
    {
        TMatch currentMatch = childMapper.primaryMatches[i];
        uint32_t readId         = currentMatch.readId;
        uint32_t origReadId     = readId;
        if(disOptions.filterType != NONE)
            origReadId = disOptions.origReadIdMap[disOptions.currentBinNo][readId];

        setMapped(mainMapper.ctx, origReadId);
        setMinErrors(mainMapper.ctx, origReadId, currentMatch.errors);

        if (!isPaired(mainMapper.ctx, origReadId) && isPaired(childMapper.ctx, readId)) // First paired match
        {
            setPaired(mainMapper.ctx, origReadId);
            mainMapper.primaryMatchesProbs[origReadId] = childMapper.primaryMatchesProbs[readId];
        }

        if(getMinErrors(mainMapper.ctx, origReadId) == currentMatch.errors)
        {
            disOptions.collectedCigars[origReadId] = childMapper.cigars[readId];
        }
    }
    stop(mainMapper.timer);
    disOptions.moveCigars += getValue(mainMapper.timer);
}

// ----------------------------------------------------------------------------
// Function transferCigars()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TMainConfig>
inline void transferCigars(Mapper<TSpec, TMainConfig> & mainMapper, DisOptions & disOptions)
{
    start(mainMapper.timer);

    typedef MapperTraits<TSpec, TMainConfig>         TTraits;
    typedef typename TTraits::TCigarsPos             TCigarsPos;
    typedef typename TTraits::TThreading             TThreading;

    resize(mainMapper.cigars.limits, getReadsCount(mainMapper.reads.seqs)+1, 0);
    for(auto iter = disOptions.collectedCigars.begin(); iter != disOptions.collectedCigars.end(); ++iter)
    {
        mainMapper.cigars.limits[iter->first + 1] = length(iter->second);
        append(mainMapper.cigarString, iter->second);
    }
    partialSum(mainMapper.cigars.limits, TThreading());
    assign(mainMapper.cigars.positions, prefix(mainMapper.cigars.limits, length(mainMapper.cigars.limits) - 1));

    // If only the primary matches were aligned, we use the identity modifier
    setHost(mainMapper.primaryCigars, mainMapper.cigars);
    setCargo(mainMapper.primaryCigars, mainMapper.primaryCigarPositions);
    assign(mainMapper.primaryCigarPositions, seqan::Range<TCigarsPos>(0, length(mainMapper.primaryMatches)), Exact());

    stop(mainMapper.timer);
    disOptions.moveCigars += getValue(mainMapper.timer);
}



template<size_t minErrors, size_t maxErrors,
         typename TIndex, typename TContigSeqs, typename TMatch>
int testReadOcc(TIndex & index, TContigSeqs & text, TMatch & match, uint32_t len, int threshold, bool const edit, bool const editMappa)
{
    int mErrors = maxErrors;
    uint64_t dist = 3 * mErrors;

//     std::vector<SHit> hits;
    std::set<uint64_t> hits;
    uint64_t occSize = 0;

    auto delegate = [&hits](auto & it, DnaString const & needle, uint8_t errors)
    {
        for(uint64_t i = it.fwdIter.vDesc.range.i1; i < it.fwdIter.vDesc.range.i2; ++i){
            hits.insert(i);
        }
//         for (auto occ : getOccurrences(iter)){
//             SHit me(occ);
//             hits.push_back(me);
//         }
    };

    StringSet<Dna5String> readOcc;

//      std::cout << getMember(*newIt, Errors()) << "\t < 5 , " << getMember(*newIt, ContigBegin()) << " >\t" << getMember(*newIt, ContigEnd()) << "\n";

    int64_t seqNo = getMember(match, ContigId());
    int64_t seqOffset = getMember(match, ContigBegin());
    int64_t seqOffsetEnd = seqOffset + len;//getMember(match, ContigEnd());
    bool rC = onReverseStrand(match);
    if(edit) //reverse this //TODO revert this !!!!
        seqOffsetEnd += maxErrors;
    if(!edit){ //edit
        Dna5String part = infix(text[seqNo], seqOffset, seqOffsetEnd);
        appendValue(readOcc, part);
    }
    else
    {
        if(seqOffset < mErrors || length(text[seqNo]) < seqOffsetEnd)
            return(666);
        for(int64_t off = -mErrors; off <= mErrors; ++off){
            Dna5String part = infix(text[seqNo], seqOffset + off, seqOffsetEnd + off);
            appendValue(readOcc, part);
        }
    }
/*
    Dna5String part = infix(text[seqNo], seqOffset, seqOffsetEnd);
    std::cout << "Search occ <" <<  seqNo << ", " << seqOffset << ">" << "\tRC: " << rC << "\n" << part << "\n";*/
    if(!edit){
        std::cout << readOcc[0] << "\n";
    }else{
//         for(int i = 0; i < length(readOcc); ++i)
//             std::cout << readOcc[i] << "\n";
        std::cout << readOcc[mErrors] << "\n";
    }

    typedef Pair <uint64_t, uint64_t>       TOcc;
    std::vector<TOcc > occs;

    if(!edit){//edit
        find<minErrors, maxErrors>(delegate, index, readOcc[0], HammingDistance());
        std::cout << hits.size() << " hits!!!!!!!!!!" << "\n";
    }
    else
    {
        if(!editMappa)
            find<minErrors, maxErrors>(delegate, index, readOcc[mErrors], HammingDistance());
        else
            find<minErrors, maxErrors>(delegate, index, readOcc[mErrors], EditDistance());

        std::cout << "Hits before del: " << hits.size() << "\n";
        occSize = 1;
//         typedef Pair <uint64_t, uint64_t>       TOcc;
//         std::vector<TOcc > occs;
        occs.reserve(hits.size());
        for(std::set<uint64_t>::iterator hitsIt = hits.begin(); hitsIt != hits.end(); ++hitsIt)
                        occs.push_back(index.fwd.sa[*hitsIt]);
        std::sort(occs.begin(), occs.end(),
                  [](TOcc & x, TOcc & y){
                        if(getSeqNo(x) == getSeqNo(y))
                            return getSeqOffset(x) < getSeqOffset(y);
                        else
                        return getSeqNo(x) < getSeqNo(y);
                    });

        uint64_t prev = 0;
        for(uint64_t i = 1; i < occs.size(); ++i){
            if(!(getSeqNo(occs[i]) == getSeqNo(occs[prev]) &&
                getSeqOffset(occs[i]) + dist >= getSeqOffset(occs[prev]) &&
                getSeqOffset(occs[i]) <= getSeqOffset(occs[prev]) + dist))
            {
                prev = i;
                ++occSize;
            }
        }
        std::cout << occSize << " hits!!!!!!!!!!" << "\n";


//         std::sort(hits.begin(), hits.end(), sHit_smaller);
//         hits.erase(std::unique(hits.begin(), hits.end(), sHit_similar<2 * maxErrors>), hits.end());
//         std::cout << hits.size() << " hits!!!!!!!!!!" << "\n";

        int k = 0;
        while(/*hits.size()*/occSize < threshold && k < length(readOcc)){
            hits.clear();
            if(!editMappa)
                find<minErrors, maxErrors>(delegate, index, readOcc[k], HammingDistance());
            else
                find<minErrors, maxErrors>(delegate, index, readOcc[k], EditDistance());
            ++k;

            //use sort and erase

            std::cout << "Hits before del: " << hits.size() << "\n";
            occSize = 1;
            occs.clear();
            occs.reserve(hits.size());
            for(std::set<uint64_t>::iterator hitsIt = hits.begin(); hitsIt != hits.end(); ++hitsIt)
                            occs.push_back(index.fwd.sa[*hitsIt]);
            std::sort(occs.begin(), occs.end(),
                    [](TOcc & x, TOcc & y){
                            if(getSeqNo(x) == getSeqNo(y))
                                return getSeqOffset(x) < getSeqOffset(y);
                            else
                            return getSeqNo(x) < getSeqNo(y);
                        });

            prev = 0;
            for(uint64_t i = 1; i < occs.size(); ++i){
                if(!(getSeqNo(occs[i]) == getSeqNo(occs[prev]) &&
                    getSeqOffset(occs[i]) + dist >= getSeqOffset(occs[prev]) &&
                    getSeqOffset(occs[i]) <= getSeqOffset(occs[prev]) + dist))
                {
                    prev = i;
                    ++occSize;
                }
            }
            std::cout << occSize << " hits!!!!!!!!!!e" << "\n";

//             std::cout << "Hits before del: " << hits.size() << "\n";
//             std::sort(hits.begin(), hits.end(), sHit_smaller);
//             hits.erase(std::unique(hits.begin(), hits.end(), sHit_similar<2 * maxErrors>), hits.end());
//             std::cout << hits.size() << " hits!!!!!!!!!!e" << "\n";
        }
    }

//     int nhits = (int)occSize;//hits.size();
//     if(occSize < 10){
//         std::cout << "Something went wrong2" << "\n";
//         std::cout << "Print occs\n";
//         for(int i = 0; i < occs.size(); ++i)
//             std::cout << occs[i] << "\n";
//     }
    return occSize;
}





template<typename TIndex, typename TContigSeqs, typename TMatch>
int testReadOcc(TIndex & index, TContigSeqs & text, TMatch & match, uint8_t maxErrors, uint32_t len, int threshold, bool const edit, bool const editMappa)
{
    int hits = 0;

    switch (maxErrors)
    {
        case 0: hits = testReadOcc<0, 0>(index, text, match, len, threshold, edit, editMappa);
                break;
        case 1: hits = testReadOcc<0, 1>(index, text, match, len, threshold, edit, editMappa);
                break;
        case 2: hits = testReadOcc<0, 2>(index, text, match, len, threshold, edit, editMappa);
                break;
        case 3: hits = testReadOcc<0, 3>(index, text, match, len, threshold, edit, editMappa);
                break;
        case 4: hits = testReadOcc<0, 4>(index, text, match, len, threshold, edit, editMappa);
                break;
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                std::exit(1);
    }
    return hits;
}

template<typename TContigSeqs, typename TMatch>
void printOcc(TContigSeqs & text, TMatch & match, int64_t ex = 0)
{
//     typedef typename InfixOnValue<TContigSeqs const>::Type         TInfix;

    StringSet<Dna5String> readOcc;
    int64_t seqNo = getMember(match, ContigId());
    int64_t seqOffset = getMember(match, ContigBegin());
    int64_t seqOffsetEnd = getMember(match, ContigEnd()); // seqOffset + len;

    Dna5String part = infix(text[seqNo], seqOffset - ex, seqOffsetEnd + ex);

    if(onReverseStrand(match))
    {
        Dna5StringReverseComplement revc(part);
        appendValue(readOcc, revc);
    }
    else
    {
        appendValue(readOcc, part);
    }
    std::cout << readOcc[0] << "\n";
}

template<typename TMatch>
inline bool matchSmaller(TMatch const & a, TMatch const & b){
    return getSortKey(a, ContigBegin()) < getSortKey(b, ContigBegin());
}

template<typename TMatch>
inline uint32_t getReadIdOSS(TMatch const & a)
{
    return(a.readId);
}

template<typename TMatch>
inline uint8_t getErrorsOSS(TMatch const & a)
{
    return(a.errors);
}

template<typename TSpec>
inline bool isSimilar(Match<TSpec> const & a, Match<TSpec> const & b, uint8_t errors)
{
    int64_t offSetA = getMember(a, ContigBegin());
    int64_t offSetB = getMember(b, ContigBegin());
/*
    if(errors > 3)
    {
        std::cout << "large Distance" << "\n";
        std::cout << "errors " << (int)errors << "\n";
        std::cout << "offSetA " << offSetA << "\n";
        std::cout << "offSetB " << offSetB << "\n";
        std::cout << "des " << (offSetA + errors >= offSetB && offSetB + errors >= offSetA) << "\n";
        std::cout << "oth " << (contigEqual(a, b) && strandEqual(a, b)) << "\n";
    }*/

    return contigEqual(a, b) && strandEqual(a, b) && offSetA + errors >= offSetB && offSetB + errors >= offSetA;

}


template <typename TSpec, typename TConfig, typename TSpec2, typename TConfig2>
void compareHits(Mapper<TSpec, TConfig> & me,
                 Mapper<TSpec2, TConfig2> & me2,
                 uint8_t maxError,
                 uint8_t strata,
                 DisOptions & disOptions)
{
            for(int i = 0; i < length(me.matchesSetByCoord); ++i){
                auto const & matches = me.matchesSetByCoord[i];
                auto matchIt = begin(matches, Standard());
                auto matchEnd = end(matches, Standard());
                while(matchIt != matchEnd){
                    write(std::cout, *matchIt);
                    ++matchIt;
                }
            }

            std::cout << "Print OSS Matches" << "\n\n\n";
            for(int i = 0; i < length(me2.matchesSetByCoord); ++i){
                auto const & matches = me2.matchesSetByCoord[i];
                auto matchIt = begin(matches, Standard());
                auto matchEnd = end(matches, Standard());
                while(matchIt != matchEnd){
                    write(std::cout, *matchIt);
                    ++matchIt;
                }
            }


            std::cout << "\n\n";
//             bool wrong = false;

            std::cout << "compare lengths: " << length(me.matchesSetByCoord) << "\t" << length(me2.matchesSetByCoord) << "\n";

            for(int i = 0; i < length(me.matchesSetByCoord); ++i){
                auto const matches = me.matchesSetByCoord[i];
                auto matchIt = begin(matches, Standard());
                auto matchEnd = end(matches, Standard());

                //get min Error
                uint32_t readId = getReadIdOSS(*matchIt);
                uint8_t minErrors = getMinErrors(me.ctx, readId);

                auto const matchesOSS = me2.matchesSetByCoord[i];
                auto matchItOSS = begin(matchesOSS, Standard());
                auto matchEndOSS = end(matchesOSS, Standard());


                while(matchIt != matchEnd){
                    bool same = isDuplicate(*matchIt, *matchItOSS, ContigBegin());
                    bool similar = isSimilar(*matchIt, *matchItOSS, maxError);
                    if(!isSimilar(*matchIt, *matchItOSS, 2 * maxError)){

                        //Yara Strata Error
                        if(getErrorsOSS(*matchIt) > strata + minErrors){
//                             std::cout << "Yara outside of strata\n";
//                             write(std::cout, *matchIt);
                            ++matchIt;
                            continue;
                        }

                        //test if yara missed the hit
                        auto matchIt_temp = matchItOSS;
                        bool smaller = true;
                        bool same2 = false;
                        while(matchIt_temp != matchEndOSS && !same2 && smaller){
                            ++matchIt_temp;
                            same2 = isSimilar(*matchIt, *matchIt_temp, 2 * maxError); //isDuplicate(*matchIt, *matchIt_temp, ContigBegin());
                            if(same2){
                                break;
                            }
                            smaller = matchSmaller(*matchIt_temp, *matchIt); //TODO move into while
                        }
                        //Skip missed hits
                        if(same2){
                            while(matchItOSS != matchIt_temp){
                                std::cout << "Yara missed\n";
                                write(std::cout, *matchItOSS);
                                std::cout << "needle: \n   " << me.reads.seqs[readId] << "\n";
                                printOcc(me.contigs.seqs, *matchItOSS, maxError);
                                ++matchItOSS;
                            }
/*
                            std::cout << "found Same:\n";
                            write(std::cout, *matchIt);
                            write(std::cout, *matchItOSS);*/
                        }
                        else
                        //verify if hit was missed because of mappability
                        {
                            std::cout << "Need to verify this: " << "\n";
                            write(std::cout, *matchIt);
                            std::cout << "Compared to this one: " << "\n";
                            write(std::cout, *matchItOSS);
//
                            Dna5String tneedle = me.reads.seqs[readId];
                            Dna5StringReverseComplement revN(tneedle);
                            std::cout << "needle: \n  " << tneedle << "\n  " << revN << "\n";

                            //TODO use true for editMappability
                            int nhits = testReadOcc(me.biIndex, me.contigs.seqs, *matchIt, maxError, disOptions.readLength, disOptions.threshold, true, false); //TODO add Threshold as input option for me
                            bool wrong = false; //revert this
                            if(nhits < disOptions.threshold)
                                wrong = true;

                            if(wrong)
                                std::cout << "Something went wrong\n";
                            ++matchIt;
                        }
                    }
                    else
                    {
                        if(!same)
                        {
                            if(similar){
                                std::cout << "Only Similar\n";
                            }
                            else
                            {
                                std::cout << "Only large distance Similar\n";
                            }
                            write(std::cout, *matchIt);
                            write(std::cout, *matchItOSS);
                            std::cout << "needle: \n   " << me.reads.seqs[readId] << "\n";
                            printOcc(me.contigs.seqs, *matchIt, maxError);
                            printOcc(me.contigs.seqs, *matchItOSS, maxError);

                        }
                        ++matchIt;
                        ++matchItOSS;

                        if(matchItOSS == matchEndOSS)
                            --matchItOSS;
                    }
                }
            }

}

/*
template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void trimHit(Mapper<TSpec, TConfig> & me, TReadSeqs & readSeqs)
{
    // trimming does not work everytime ex: GTCCCCCCCCGCTA
    //                                      ACCCCCCCCGCTA -> only after comparing C to G it is clear we need and insertion
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


// ----------------------------------------------------------------------------
// Function inTextVerification()
// ----------------------------------------------------------------------------

template<typename TSpec, typename TConfig, typename TMatch, typename TNeedle, typename TError>
inline bool inTextVerification(Mapper<TSpec, TConfig> & me, TMatch & match, TNeedle const & needle, TError maxErrors)
{
    typedef MapperTraits<TSpec, TConfig>                        TTraits;
    typedef typename TConfig::TContigsLen                       TContigsLen;
    typedef typename TConfig::TContigsSize                      TContigsSize;

    typedef typename TTraits::TContigSeqs                       TContigSeqs;
    typedef typename InfixOnValue<TContigSeqs const>::Type      TContigSeqsInfix;
    typedef typename InfixOnValue<TNeedle const>::Type          TNeedleInfix;

    TContigSeqs & contigSeqs = me.contigs.seqs;

//     std::cout << "Needle: \n" << needle << "\n";
    TContigsSize seqNo = getMember(match, ContigId());
    TContigsLen seqOffset = getMember(match, ContigBegin());
    TContigsLen seqOffsetEnd = getMember(match, ContigEnd());;
    TContigSeqsInfix text = infix(contigSeqs[seqNo], seqOffset, seqOffsetEnd);
//     std::cout << text << "\n";

    typedef Finder<TContigSeqsInfix>                        TFinder;
    typedef AlignTextBanded</*FindPrefix*/FindInfix,
                            NMatchesNone_,
                            NMatchesNone_>                  TMyersSpec;
    typedef Myers<TMyersSpec, True, void>                   TAlgorithm;
    typedef PatternState_<TContigSeqsInfix, TAlgorithm>     TPattern;
    TPattern pattern;

    int minErrors = maxErrors + 1;
    TFinder finder(text);

    while (find(finder, needle, pattern, -static_cast<int>(maxErrors + 1)))
    {
        int currentErrors = -getScore(pattern);
        if(minErrors > currentErrors)
            minErrors = currentErrors;

//         std::cout << "cScore: " << currentErrors << "\t";
    }

//     std::cout << "\nScore: " << (int)minErrors << "\n";
    if(minErrors <= maxErrors)
        match.errors = minErrors;

    return(minErrors <= maxErrors);
}

// ----------------------------------------------------------------------------
// Function inTextVerification()
// ----------------------------------------------------------------------------

template<typename TSpec, typename TConfig, typename TMatch, typename TNeedle, typename TError>
inline bool inTextVerificationE(Mapper<TSpec, TConfig> & me, TMatch & match, TNeedle const & needle, TError maxErrors, bool const ossMatch = false)
{
    typedef MapperTraits<TSpec, TConfig>                        TTraits;
    typedef typename TConfig::TContigsLen                       TContigsLen;
    typedef typename TConfig::TContigsSize                      TContigsSize;


    typedef typename TTraits::TContigSeqs                       TContigSeqs;
    typedef typename Value<TContigSeqs const>::Type              TContigSeq;
//     typedef typename Value<TContigSeqs>::Type                       TContig;
    typedef typename StringSetPosition<TContigSeqs const>::Type TContigPos;
    typedef typename InfixOnValue<TContigSeq const>::Type      TContigSeqInfix;
    typedef typename InfixOnValue<TNeedle const>::Type          TNeedleInfix;
    typedef ModifiedString<TNeedle, ModReverse>                 TNeedleInfixRev;
    typedef ModifiedString<TContigSeqInfix, ModReverse>        TStringInfixRev;

    typedef String<GapAnchor<int> >            TGapAnchors;
    typedef AnchorGaps<TGapAnchors>            TAnchorGaps;
    typedef Gaps<TContigSeqInfix, TAnchorGaps>     TContigGaps;
    typedef Gaps<TNeedle, TAnchorGaps>             TReadGaps;

    TContigSeqs & contigSeqs = me.contigs.seqs;

    TContigsSize seqNo = getMember(match, ContigId());
    TContigsLen seqOffset = getMember(match, ContigBegin());
    TContigsLen seqOffsetEnd = getMember(match, ContigEnd());


    TContigSeqInfix text = infix(contigSeqs[seqNo], seqOffset, seqOffsetEnd);


    // Fix rare edge Case
    if (seqOffset == 0 && length(text) <= length(needle))
    {
        TGapAnchors contigAnchors;
        TGapAnchors readAnchors;
        TContigGaps contigGaps(text, contigAnchors);
        TReadGaps needleGaps(needle, readAnchors);

//         TContigSeqInfix text_c = infix(contigSeqs[seqNo], seqOffset, seqOffsetEnd);
//         TContigGaps contigGaps(text_c);
//         Gaps<TNeedle> needleGaps(needle);

        Score<int, Simple> scoringScheme(0, -999, -1000);
        int score = globalAlignment(contigGaps, needleGaps, scoringScheme, AlignConfig<true, false, false, true>()) / -999;
        clipSemiGlobal(contigGaps, needleGaps);
/*
        std::cout << contigGaps << "\n" << needleGaps << "\n";
        std::cout << beginPosition(contigGaps) << "\n" << endPosition(contigGaps) << "\n";
        std::cout << beginPosition(needleGaps) << "\n" << endPosition(needleGaps) << "\n";
        std::cout << "Score " << (int)score << "\n";
        std::cout << "Before:\n";
        write(std::cout, match);*/

        setErrors(match, score);

        TContigPos contigBegin(getMember(match, ContigId()), getMember(match, ContigBegin()));
        TContigPos contigEnd = contigBegin;
        contigEnd = posAdd(contigEnd, endPosition(contigGaps));
        contigBegin = posAdd(contigBegin, beginPosition(contigGaps));
        setContigPosition(match, contigBegin, contigEnd);
//         write(std::cout, match);
        return score <= maxErrors;
    }

    typedef Finder<TContigSeqInfix>                        TFinder;
    typedef Finder<TStringInfixRev>                       TFinder2;
    typedef AlignTextBanded</*FindPrefix*/FindInfix,
                            NMatchesNone_,
                            NMatchesNone_>                     TMyersSpecInfix;
    typedef Myers<TMyersSpecInfix, True, void>                 TAlgorithmInfix;
    typedef PatternState_<TContigSeqInfix, TAlgorithmInfix>   TPatternInfix;
    TPatternInfix patternInfix;

    typedef AlignTextBanded<FindPrefix,
                            NMatchesNone_,
                            NMatchesNone_>               TMyersSpec;
    typedef Myers<TMyersSpec, True, void>                TAlgorithm;
    typedef PatternState_<TNeedle, TAlgorithm>  TPattern;
    TPattern pattern;
    typedef PatternState_<TNeedleInfixRev, TAlgorithm>  TPatternRev;
    TPatternRev patternRev;

    int minErrors;
    // for ossMatch we know match is valid and only need to determine start and end position of the match (we added overlap so suffix filter works)
    if(!ossMatch){
        minErrors = maxErrors + 1;
        TFinder finderInfix(text);

        while (find(finderInfix, needle, patternInfix, -static_cast<int>(maxErrors + 1)))
        {
            int currentErrors = -getScore(patternInfix);
            if(minErrors > currentErrors)
                minErrors = currentErrors;
        }

        if(minErrors <= maxErrors)
            setErrors(match, minErrors);
//             match.errors = minErrors;
        else
            return false;
    }
    else
    {
        minErrors = 0;
    }


    TFinder finder(text);
    int mErrors = maxErrors * 4;
    TContigsLen endPos = 0;
    int tmpErrors = 99999;
    while (find(finder, needle, pattern, -static_cast<int>(maxErrors * 4)))
    {
        int currentEnd = position(finder) + 1;
        int currentErrors = -getScore(pattern);
        if (currentErrors <= tmpErrors)
            tmpErrors = currentErrors;
//         if (getValue(text, currentEnd) != back(needle))
//             ++currentErrors;
//         std::cout << currentErrors << "\t" << currentEnd << "\n";
        if (currentErrors <= mErrors)
        {
            mErrors = currentErrors;
            endPos = currentEnd;
        }
    }

    SEQAN_ASSERT_NOT(endPos == 0);
/*
    if(endPos == 0){
        std::cout << "EndPos failed\n" << "Current: " << tmpErrors << "\n";
        write(std::cout, match);
        std::cout << text << "\n" << needle << "\n";
    }*/

//     TContigSeqInfix infixPrefix = infix(text, 0, endPos);
    TStringInfixRev infixRev(text);
    TNeedleInfixRev needleRev(needle);
    TFinder2 finder2(infixRev);

    mErrors = maxErrors * 4;
    TContigsLen startPos = endPos;

    while (find(finder2, needleRev, patternRev, -static_cast<int>(maxErrors * 4)))
    {
        int currentEnd = position(finder2) + 1;
        int currentErrors = -getScore(patternRev);
//         std::cout << currentErrors << "\t" << currentEnd << "\n";
        if (currentErrors <= mErrors)
        {
            mErrors = currentErrors;
            startPos = currentEnd;
        }
    }


    //there is no need to report the minimum error or that a read mapped since this case in already in-text-Verification inside the optimal search schemes.
    TContigPos contigBegin(getMember(match, ContigId()), getMember(match, ContigBegin()));
    TContigPos contigEnd = contigBegin;
    contigEnd = posAdd(contigEnd, endPos);
    contigBegin = posAdd(contigBegin, length(text) - startPos);
    setContigPosition(match, contigBegin, contigEnd);

    return(minErrors <= maxErrors);
}


// ----------------------------------------------------------------------------
// Function _mapReadsImpl()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig, typename TReadSeqs, typename TSeqsSpec>
inline void _mapReadsImpl(Mapper<TSpec, TConfig> & me,
                          Mapper<TSpec, TMainConfig>  & mainMapper,
                          auto & mybiIndex,
                          StringSet<TReadSeqs, TSeqsSpec> & readSeqs,
                          DisOptions & disOptions)
{
    uint32_t len;
    if(disOptions.readLength != 0)
        len = disOptions.readLength;
    else
        len = length(readSeqs[0]);

    disOptions.readLength = len;

    me.maxError = std::ceil(me.options.errorRate * len);
    me.strata = std::ceil(disOptions.strataRate * len);
    Mapper<void, TConfig> me2(disOptions);

    //tracking
    uint32_t itvJobsDone = 0;

    if(!disOptions.ossOff){

    initReadsContext(me, readSeqs);

    typedef MapperTraits<TSpec, TConfig>                        TTraits;
    typedef typename TTraits::TMatchesAppender                  TMatchesAppender;
    typedef Delegate<TSpec, TConfig>                            Delegate;
    typedef DelegateDirect<TSpec, TConfig>                      DelegateDirect;
//     typedef DelegateUnfiltered<TSpec, TConfig>                  DelegateUnfiltered;

    typedef typename TTraits::TContigSeqs                       TContigSeqs;

/*
    typedef typename TTraits::TIndexConfig                      TIndexConfig;
    TIndexConfig::SAMPLING = me.samplingRate;

    typedef FMIndex<void, TIndexConfig>                             TIndexSpec;
    typedef BidirectionalIndex<TIndexSpec>                          TBiIndexSpec;
    typedef Index<typename TIndexConfig::Text, TIndexSpec>          TIndex;
    typedef Index<typename TIndexConfig::Text, TBiIndexSpec>        TBiIndex;

    TBiIndex & mybiIndex = me.biIndex;

    */

    TMatchesAppender appender(me.matchesByCoord);
    bool noOverlap = disOptions.noSAfilter && (disOptions.determineExactSecondaryPos || disOptions.noDelayITV || disOptions.hammingDistance);
    Delegate delegate(appender, noOverlap, me.maxError);
    DelegateDirect delegateDirect(appender);
//     DelegateUnfiltered delegateUnfiltered(appender, noOverlap);
    TContigSeqs & contigSeqs = me.contigs.seqs;

    OSSContext<TSpec, TConfig> ossContext(me.ctx, appender, readSeqs, contigSeqs);

    if(disOptions.verbose > 1){
        std::cout << "Mapping " << length(readSeqs) << " reads\n";
        std::cout << "maxError: " << (int)me.maxError << "\t" << (int)me.strata << "\n";
    }

    if(!disOptions.noMappability){
        start(me.timer);
        CharString bPath = me.options.mappabilitySubDirectory;
        bPath += "/";
        if(disOptions.verbose > 1)
            std::cout << "\nLoading Bitvectors: " << bPath << "\n";
        loadAllBitvectors(bPath, me.bitvectors, me.bitvectorsMeta, len, (disOptions.verbose > 1));

        if(!me.bitvectors.empty() && disOptions.verbose > 1){
            std::cout << "Bit vectors loaded. Number: " << me.bitvectors.size() << "\n";
            std::cout << "Index Size: " << length(me.biIndex.fwd.sa) << "\n";
            std::cout << "Number of Sequences: " << length(me.contigs.seqs) << "\n";

            std::cout << "Length of Bitvectors: " << me.bitvectors[0].first.size() << "\n";
            ossContext.bitvectorsMeta = me.bitvectorsMeta;
        }
        stop(me.timer);
        me.stats.loadingBitvectors += getValue(me.timer);
        if(disOptions.verbose > 1)
            std::cout << "Loading Bitvectors time:\t\t\t" << me.timer << std::endl;
    }

// copy parameters to ossContext

//     YaraFMConfig<uint16_t, uint32_t, uint32_t> myConfig{}; //myConfig.SAMPLING
    if(disOptions.verbose > 1)
        std::cout << "SAMPLING RATE: " << me.samplingRate << "\n";

    ossContext.loadInputParameters(me.maxError, me.strata, len, length(me.contigs.seqs), me.samplingRate, disOptions.fmTreeThreshold, disOptions.fmTreeBreak);
    ossContext.itv = !disOptions.noITV;
    ossContext.normal.suspectunidirectional = false;
//     ossContext.saFilter = !disOptions.noSAfilter;
    ossContext.delayITV = !disOptions.noDelayITV;
    ossContext.earlyLeaf = disOptions.earlyLeaf;
    ossContext.itvOccThreshold = disOptions.itvOccThreshold;
    ossContext.noSAfilter = disOptions.noSAfilter;
    ossContext.hammingDpercentage = disOptions.hammingDpercentage;

    start(me.timer);

    if(disOptions.hammingDistance)
        find(0, me.maxError, me.strata, ossContext, delegate, delegateDirect, mybiIndex, me.bitvectors, readSeqs, HammingDistance());
    else
        find(0, me.maxError, me.strata, ossContext, delegate, delegateDirect, mybiIndex, me.bitvectors, readSeqs, EditDistance());


    if (me.options.verbose > 0)
    {
        typedef MapperTraits<TSpec, TMainConfig>                    TTraits2;
        unsigned long mappedReads = count(me.ctx.mapped, true, typename TTraits2::TThreading());
        me.stats.mappedReads += mappedReads;

        if (me.options.verbose > 0)
            std::cout << "1Mapped reads:\t\t\t" << mappedReads << std::endl;
    }

    stop(me.timer);
    me.stats.optimumSearch += getValue(me.timer);
    if(disOptions.verbose > 0)
    {
        std::cout << "\nOptimum Search time:\t\t" << me.timer << std::endl;
        std::cout << "Matches count:\t\t\t" << lengthSum(me.matchesByCoord) << std::endl;
        std::cout << "Default Locates:\t\t" << ossContext.defaultLocates << std::endl;
        std::cout << "FM-Tree Locates:\t\t" << ossContext.fmtreeLocates << std::endl;
        std::cout << "FM-Tree FM-LF-Mappings:\t\t" << ossContext.fmtreeBacktrackings << std::endl;
        std::cout << "FM-Tree Thres LF-Mappings:\t\t" << ossContext.fmtreeBreakLocates << std::endl;

    }


    //save OSS hits
    if(disOptions.compare && !ossContext.delayITV){
        me2.matchesByCoord = me.matchesByCoord;
        me2.matchesSetByCoord = me.matchesSetByCoord;
        clear(me.matchesSetByCoord);
        clear(me.matchesByCoord);
        aggregateMatchesOSS(me2, readSeqs);
    }
    else if(!ossContext.delayITV)
    {
        std::cout << "ITV during Search\n";
        aggregateMatchesOSS(me, readSeqs);

    }
    else
    {/*
        aggregateMatchesITV_debug(me, readSeqs);
        std::cout << "Print duplicate Matches\n";
        for(int i = 0; i < 20; ++i){
//             std::cout << "Is Read mapped: " << isMapped(ossContext.ctx, i) << "\n";
//             std::cout << "ReadmapperCont: " << isMapped(me.ctx, i) << "\n";
            auto const & matches = me.matchesSetByCoord[i];
            auto matchIt = begin(matches, Standard());
            auto matchEnd = end(matches, Standard());
            while(matchIt != matchEnd){
                write(std::cout, *matchIt);
                ++matchIt;
            }
        }*/

         //ITV after filtering
        //filtering after cordinates + Intervalsize
        aggregateMatchesITV(me, readSeqs);

/*

        std::cout << "Print Matches\n";
        for(int i = 0; i < 20; ++i){
//             std::cout << "Is Read mapped: " << isMapped(ossContext.ctx, i) << "\n";
//             std::cout << "ReadmapperCont: " << isMapped(me.ctx, i) << "\n";
            auto const & matches = me.matchesSetByCoord[i];
            auto matchIt = begin(matches, Standard());
            auto matchEnd = end(matches, Standard());
            while(matchIt != matchEnd){
                write(std::cout, *matchIt);
                ++matchIt;
            }
        }*/

        std::cout << "Hits after filtering: " << lengthSum(me.matchesSetByCoord) << "\n\n";

        typedef typename TConfig::TContigsLen                       TContigsLen;
        typedef typename TConfig::TContigsSize                      TContigsSize;


        start(me.timer);
        uint32_t valids = 0;
        uint32_t oss = 0;
        uint32_t merges = 0;
        uint32_t dups = 0;

        #pragma omp parallel for schedule(dynamic) num_threads(disOptions.threadsCount)
        for(int i = 0; i < length(me.matchesSetByCoord); ++i){
//                 std::cout << "New read" << i << "\n";
            auto const & matches = me.matchesSetByCoord[i];
            auto matchIt = begin(matches, Standard());
            auto matchEnd = end(matches, Standard());
            auto lValid = matchIt;
            int extendedLength = 0;

            while(matchIt != matchEnd){
                bool ossMatch = getMember(*matchIt, Errors()) <= me.maxError;

                oss += ossMatch;

                //check Neighborhood
                if(matchIt != lValid && extendedLength <= (me.maxError * 10) &&
                   !disOptions.determineExactSecondaryPos &&
                   getMember(*matchIt, ContigId()) == getMember(*lValid, ContigId()) &&
                   !(onForwardStrand(*matchIt) ^ onForwardStrand(*lValid)) &&
                   getMember(*lValid, ContigBegin()) + 2 * me.maxError >= getMember(*matchIt, ContigBegin()))
                {
                    bool valid = true;
                    if(!ossMatch)
                    {
                        uint32_t readSeqId = getReadSeqId(*matchIt, readSeqs);
                        uint32_t readId = getReadId(readSeqs, readSeqId);
                        valid = inTextVerification(me, *matchIt, readSeqs[readSeqId], me.maxError);
                        ++itvJobsDone;
                        if(valid){
                            setMapped(me.ctx, readId);
                            setMinErrors(me.ctx, readId, getMember(*matchIt, Errors()));
                        }
                    }
                    // extend the last valid match if the current match is valid
                    if(valid){
                        int32_t shift = static_cast<int32_t>(getMember(*matchIt, ContigEnd())) - getMember(*lValid, ContigEnd());
                        if(shift > 0)
                        {
                            extendedLength += shift;
                            //merge oss and ITV if extendedLength is under threshold into the lastValid match
                            shiftEnd(*lValid, shift);
                        }
                        ++merges;
                    }
                    setInvalid(*matchIt);
                    ++dups;
                }
                else
                {
                    // match is not in Neighborhood of another match (to the left) do ITV or we want to determine the exact position of all (repeating) matches
                    if(!ossMatch || !disOptions.noSAfilter)
                    {
                        uint32_t readSeqId = getReadSeqId(*matchIt, readSeqs);
                        uint32_t readId = getReadId(readSeqs, readSeqId);
                        bool valid;
                        if(disOptions.determineExactSecondaryPos || !disOptions.noSAfilter)
                        {
//                             write(std::cout, *matchIt);

                            if(!disOptions.alignSecondary)
                                valid = inTextVerificationE(me, *matchIt, readSeqs[readSeqId], me.maxError, ossMatch);
                            else
                                valid = true;

//                                 write(std::cout, *matchIt);
                        }
                        else
                        {
                            valid = inTextVerification(me, *matchIt, readSeqs[readSeqId], me.maxError);
                        }
                        if(!ossMatch)
                            ++itvJobsDone;
                        if(valid){
                            setMapped(me.ctx, readId);
                            setMinErrors(me.ctx, readId, getMember(*matchIt, Errors()));
                            ++valids;
                            lValid = matchIt;
                            extendedLength = 0;
                        }
                        else
                        {
                            setInvalid(*matchIt);
                        }
                    }
                    else
                    {
                        lValid = matchIt;
                        extendedLength = 0;
                    }
                }
                ++matchIt;
            }
        }
/*
        std::cout << "Print aligned Matches\n";
        for(int i = 0; i < 1; ++i){
//             std::cout << "Is Read mapped: " << isMapped(ossContext.ctx, i) << "\n";
//             std::cout << "ReadmapperCont: " << isMapped(me.ctx, i) << "\n";
            auto const & matches = me.matchesSetByCoord[i];
            auto matchIt = begin(matches, Standard());
            auto matchEnd = end(matches, Standard());
            while(matchIt != matchEnd){
                write(std::cout, *matchIt);
                ++matchIt;
            }
        }*/


        stop(me.timer);
        me.stats.inTextVerification += getValue(me.timer);

        if(disOptions.verbose > 1){
            std::cout << "Hits accepted: " << valids << " OSS: " << oss << " merges: " << merges << " dups: " << dups << " from: " << lengthSum(me.matchesSetByCoord) << "\n";
        }
    }
    if(disOptions.verbose > 1){
        std::cout << "Unique Matches count:\t\t\t" << lengthSum(me.matchesSetByCoord) << std::endl;

        std::cout << "filtered Occs of read on suffix array: \t" << ossContext.filteredOccsOfRead << "\n";
        std::cout << "OSS Matches: \t\t\t\t" << ossContext.delegateOcc << "\n";
        std::cout << "filtered OSS Matches: \t\t\t" << ossContext.delegateFilteredOcc << "\n";
        std::cout << "itv Attemps: \t\t\t\t" << ossContext.itvAttemps << "\n";
        std::cout << "itv Matches: \t\t\t\t" << ossContext.itvOccs << "\n";
        std::cout << "itv Jobs: \t\t\t\t" << ossContext.itvJobs << "\n";
        std::cout << "itv Jobs done: \t\t\t\t" << itvJobsDone << "\n\n";
    }

    }



    if (me.options.verbose > 1)
    {
        typedef MapperTraits<TSpec, TMainConfig>                    TTraits2;
        unsigned long mappedReads = count(me.ctx.mapped, true, typename TTraits2::TThreading());
        me.stats.mappedReads += mappedReads;

        std::cout << "2Mapped reads:\t\t\t" << mappedReads << "\n";
    }

    if(disOptions.ossOff || disOptions.compare)
    {
        std::cout << "Using Seed and Extension: \n";
        initReadsContext(me, readSeqs);
        initSeeds(me, readSeqs);

        collectSeeds<0>(me, readSeqs);
        findSeeds<0>(me, 0);
        classifyReads(me);
        collectSeeds<1>(me, readSeqs);
        collectSeeds<2>(me, readSeqs);
        findSeeds<0>(me, 1);
        findSeeds<0>(me, 2);
        rankSeeds(me);
        reserveMatches(me);
        extendHits<0>(me, 0);
        extendHits<0>(me, 1);
        extendHits<0>(me, 2);
        clearSeeds(me);
        clearHits(me);

        initSeeds(me, readSeqs);
        collectSeeds<1>(me, readSeqs);
        findSeeds<1>(me, 1);
        collectSeeds<2>(me, readSeqs);
        findSeeds<1>(me, 2);
        rankSeeds(me);
        // TODO(esiragusa): filter out hits with distance < 1.
        extendHits<1>(me, 1);
        extendHits<1>(me, 2);
        clearSeeds(me);
        clearHits(me);

        if (me.options.sensitivity > LOW)
        {
            initSeeds(me, readSeqs);
            collectSeeds<2>(me, readSeqs);
            findSeeds<2>(me, 2);
            rankSeeds(me);
            // TODO(esiragusa): filter out hits with distance < 2.
            extendHits<2>(me, 2);
            clearHits(me);
            clearSeeds(me);
        }
        if(disOptions.verbose > 0){
            std::cerr << "\n Major search time:\t\t" << me.stats.extendHits +  me.stats.findSeeds << std::endl;
            std::cerr << "Matches count:\t\t\t" << lengthSum(me.matchesByCoord) << std::endl;
        }

        aggregateMatches(me, readSeqs);

        if(disOptions.verbose > 0)
        {
            std::cerr << "Unique Matches count:\t\t\t" << lengthSum(me.matchesSetByCoord)/*length(me.matchesByCoord)*/ << std::endl;
        }

        //TODO use optimal Search Schemes again to get "correct Positions"

        if(disOptions.compare){
//             compareHits(me, me2, me.maxError, me.strata, disOptions);
        }
    }

    rankMatches(me, me.reads.seqs);
    if (me.options.verifyMatches)
        verifyMatches(me);
//     alignMatches(me);
    copyMatches(mainMapper, me, disOptions);
//     copyCigars(mainMapper, me, disOptions);
    appendStats(mainMapper, me);
}


template <typename TSpec, typename TConfig, typename TMainConfig, typename TIndex, typename TReadSeqs, typename TSeqsSpec>
inline void _mapReadsImplOSS(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig>  & mainMapper, TIndex & index,
                             StringSet<TReadSeqs, TSeqsSpec> & readSeqs, DisOptions & disOptions)
{
	//disOptions.origReadIdMap[disOptions.currentBinNo][i] //get orginal id from ith-read is done in copy matches so dont need that
    initReadsContext(me, readSeqs);


    typedef MapperTraits<TSpec, TConfig>                        TTraits;
    typedef typename TTraits::TReads::TSeq                      TRead;

    typedef typename TTraits::TSA                               TSA;
    typedef typename TTraits::TContigsPos                       TContigsPos;
    typedef typename Size<TSA>::Type                            TSAPos;
    typedef typename Value<TSA>::Type                           TSAValue;


    typedef typename TTraits::TMatchesAppender                  TMatchesAppender;

    typedef typename TTraits::TMatch                             TMatch;
    typedef typename TTraits::TReadSeq                            TReadSeq;
    typedef typename Size<TReadSeq>::Type                       TReadSeqSize;
    typedef typename Size<TReadSeqs>::Type                      TReadId;

    typedef TReadSeqSize                                         THitErrors;

    typedef typename TTraits::TContigSeqs                        TContigSeqs;

//     typedef typename TTraits::TBiIndex                          TBiIndex;
//     typedef typename TTraits::TIndex                            TIndex;




    typedef typename Iterator<StringSet<TReadSeqs, TSeqsSpec> const, Rooted>::Type TReadIt;
    typedef typename Reference<TReadIt>::Type                                     TReadRef;

    // Iterate over all reads.
    int k = 0;
    iterate(readSeqs, [&](TReadIt const & readIt)
    {
        k++;
        TReadRef it = value(readIt);
//         std::cout << "Test Iterate: " << k << "\n";
//         std::cout << it << "\n";
    }, Rooted(), typename TTraits::TThreading());



    TRead read = readSeqs[0];

    std::cout << "Test cout" << "\n";
    std::cout << read << "\n";

// //     Iterator<TIndex, TopDown<> >::Type it(index);
//     Iter<Index<>, VSTree<TopDown<> > > it(index);
//     Iter<TBiIndex, VSTree<TopDown<> > > iter(me.biIndex);


    Iter<TIndex, VSTree<TopDown<> > > iter(index);


    std::cout << iter.fwdIter.index->sa[0] << "\n";




    for(TSAPos r = iter.fwdIter.vDesc.range.i1; r < iter.fwdIter.vDesc.range.i1 + 5; ++r){
        TSAValue saPost = iter.fwdIter.index->sa[r];
        std::cout << saPost << "\n";
    }

    iter.fwdIter.vDesc.range.i1 = 20;
    iter.fwdIter.vDesc.range.i1 = 30;

    iter.revIter.vDesc.range.i1 = 30;
    iter.revIter.vDesc.range.i1 = 40;

    TMatchesAppender appender(me.matchesByCoord);
    TReadId readID = 7;
    uint8_t cerrors = 1;
    for (auto occ : getOccurrences(iter)){
        TMatch prototype;
        setReadId(prototype, readID);
        TContigsPos start = occ;
        TContigsPos end = posAdd(start, repLength(iter));
        setContigPosition(prototype, start, end); //fix this
        THitErrors errors = cerrors;
        prototype.errors = errors;
        setMapped(me.ctx, readID);
	//appendValue(appender, prototype, Generous(), typename TTraits::TThreading());
    }

//     for(TSAPos r = getValueI1(range(iter.fwdIter)); r < getValueI1(range(iter.fwdIter)) + 5; ++r)
//         std::cout << r << "\n";


    TSAValue saPos = iter.fwdIter.index->sa[10];
    TContigsPos start = saPos;
//     setSeqOffset(saValue, suffixLength(saValue, me.contigSeqs) - seedLength); //for direct modification
    TContigsPos end = posAdd(saPos, 100);

    std::cout << "Value1: " << start << "\n";
    std::cout << "Value: " << end << "\n";


//     TMatchesAppender appender(me.matchesByCoord);
    TMatch prototype;

    //check if read is reverse!!
//     TReadId readID = 7;
    setReadId(prototype, readID);
    setMapped(me.ctx, readID);
    setContigPosition(prototype, start, end);
    THitErrors errors = 2;
    prototype.errors = errors;


    std::cout << "Number of Matches: " << length(me.matchesByCoord) << "\n";
    appendValue(appender, prototype, Generous(), typename TTraits::TThreading());

    std::cout << length(me.matchesByCoord) << "\n";

    TMatch myMatch = me.matchesByCoord[0];
    std::cout << "Start Match: " << "\n";
    std::cout << myMatch.contigBegin << "\n";
    std::cout << myMatch.contigEnd << "\n";
    std::cout << myMatch.contigId << "\n";
    std::cout << myMatch.errors << "\n";
    std::cout << myMatch.readId << "\n";


    TContigSeqs & contigSeqs = me.contigs.seqs;
    OSSContext<TSpec, TConfig> ossContext(me.ctx, appender, readSeqs, contigSeqs);

/*
    typedef Delegate<TTraits>          TDelegate;
    TDelegate delegate(appender);

    delegate(ossContext, start, end, errors, readID);
    std::cout << "experimental Delegate call: "<< length(me.matchesByCoord) << "\n";*/

    //does not work since i cant use templates needed for ossContext


//     template<typename TContext>
//     auto delegate = [](TContext & ossContext, TContigsPos const & start, TContigsPos const & end, TMatchErrors errors, TReadId const needleId)
//     {
//         TMatch hit;
//         setContigPosition(hit, start, end);
//         hit.errors = errors;
//         setReadId(hit, ossContext.readSeqs, needleId);
// //         hit.readId = needleId;
//         appendValue(ossContext.matches, hit, Generous(), typename Traits::TThreading());
//     };

//     delegate(ossContext, start, end, errors, readID);




    typedef typename Id<TContigSeqs>::Type                    TContigSeqsId;
    typedef typename Value<TContigSeqs>::Type                 TContigSeqsString;
    typedef typename Size<TContigSeqsString>::Type            TContigSeqsSize;
    typedef typename TTraits::TContigsLen                     TContigsLen;
    typedef typename InfixOnValue<TContigSeqs const>::Type    TContigSeqsInfix;
//     typedef ModifiedString<THaystackInfix, ModReverse> THaystackInfixRev;




    TContigSeqsId contigSeqsId = getSeqNo(start);
    TContigsLen saLoc = getSeqOffset(start);
    TContigSeqsSize contigSeqsLength = length(contigSeqs[contigSeqsId]);
//     setSeqOffset(saValue, suffixLength(saValue, me.contigSeqs) - seedLength);
//     std::cout << "Test access to contigSeqs" << "\n";
//     std::cout << "Contig ID: " << contigSeqsId << " saLoc: " << saLoc << "\n";
//     std::cout << contigSeqs[contigSeqsId][saLoc] << "\n";
//     std::cout << contigSeqs[contigSeqsId][saLoc + 1] << "\n";
//     std::cout << contigSeqs[contigSeqsId][saLoc + 2] << "\n";
//
//     std::cout << "Test access infices of contigSeqs" << "\n";
    TContigSeqsInfix myinfix = infix(contigSeqs, start, end);
//     std::cout << myinfix << "\n";

//     IsSameType<TDistance, HammingDistance>::VALUE



/*
    goDown(iter, Fwd());
    bool delta = !ordEqual(parentEdgeLabel(iter, Fwd()), readSeqs[0][0]);
    goRight(iter, Fwd());
    delta = !ordEqual(parentEdgeLabel(iter, Fwd()), readSeqs[0][0]);

    aggregateMatches(me, readSeqs);
    rankMatches(me, me.reads.seqs);
    if (me.options.verifyMatches)
        verifyMatches(me);
    alignMatches(me);
    copyMatches(mainMapper, me, disOptions);
    copyCigars(mainMapper, me, disOptions);
    appendStats(mainMapper, me);*/
}




// ----------------------------------------------------------------------------
// Function clasifyLoadedReads()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TMainConfig, typename TFilter>
inline void clasifyLoadedReads(Mapper<TSpec, TMainConfig>  & mainMapper, TFilter const & filter, DisOptions & disOptions)
{
    if (disOptions.verbose > 1)
        std::cerr << "Clasify loaded Reads in " << disOptions.numberOfBins << " bin\n";

    start(mainMapper.timer);

    uint32_t numReads = getReadsCount( mainMapper.reads.seqs);
    uint16_t avgReadLen = lengthSum(mainMapper.reads.seqs) / (numReads * 2);
    uint16_t threshold = disOptions.getThreshold(avgReadLen);

    disOptions.origReadIdMap.clear();
    disOptions.origReadIdMap.resize(disOptions.numberOfBins);

    if (threshold == 0)
    {
        std::cerr <<"[WARNING!] 0 k-mer is required to filter a read!\n";
        std::cerr <<"All reads will pass filteration and be mapped everywhere.\n ";
        std::cerr <<"This will be extremly slow.\n ";
        std::cerr <<"Choose an approprate error rate based on kmer size and read length\n";
        for (uint32_t readID = 0; readID < numReads; ++readID)
        {
            for (uint32_t binNo = 0; binNo < disOptions.numberOfBins; ++binNo)
            {
               disOptions.origReadIdMap[binNo].push_back(readID);
            }
        }
    }
    // if paired classify only one pair
    if (IsSameType<typename TMainConfig::TSequencing, PairedEnd>::VALUE)
        numReads = getPairsCount( mainMapper.reads.seqs);

    uint32_t numThr = disOptions.threadsCount;
    uint32_t batchSize = numReads/numThr;
    if(batchSize * numThr < numReads) ++batchSize;

    std::vector<std::future<void>> tasks;


    for (uint32_t taskNo = 0; taskNo < numThr; ++taskNo)
    {
        tasks.emplace_back(std::async([=, &mainMapper, &disOptions, &filter] {
            for (uint32_t readID = taskNo*batchSize; readID < numReads && readID < (taskNo +1) * batchSize; ++readID)
            {
                std::vector<bool> selectedBins(disOptions.numberOfBins, false);
                filter.whichBins(selectedBins, mainMapper.reads.seqs[readID], threshold);
                filter.whichBins(selectedBins, mainMapper.reads.seqs[readID + numReads], threshold);

                if (IsSameType<typename TMainConfig::TSequencing, PairedEnd>::VALUE)
                {
                    filter.whichBins(selectedBins, mainMapper.reads.seqs[readID + 2*numReads], threshold);
                    filter.whichBins(selectedBins, mainMapper.reads.seqs[readID + 3*numReads], threshold);
                }

                for (uint32_t binNo = 0; binNo < disOptions.numberOfBins; ++binNo)
                {
                    if(selectedBins[binNo])
                    {
                        mtx.lock();
                        disOptions.origReadIdMap[binNo].push_back(readID);
                        mtx.unlock();
                    }
                }
            }
        }));
    }
    for (auto &&task : tasks)
    {
        task.get();
    }

    stop(mainMapper.timer);
    disOptions.filterReads += getValue(mainMapper.timer);

    if (disOptions.verbose > 1)
    {
        for (uint32_t binNo = 0; binNo < disOptions.numberOfBins; ++binNo)
        {
            std::cerr << "bin " << binNo << "\t" << disOptions.origReadIdMap[binNo].size() << std::endl;
        }
    }
}


// ----------------------------------------------------------------------------
// Function loadFilteredReads()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void loadFilteredReads(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig>  & mainMapper, DisOptions & disOptions)
{

    start(mainMapper.timer);
    if(disOptions.filterType == NONE)
    {
        for (uint32_t i = 0; i< getReadSeqsCount(mainMapper.reads.seqs); ++i)
        {
            appendValue(me.reads.seqs, mainMapper.reads.seqs[i]);
        }
    }
    else
    {
        uint32_t numReads = getReadsCount( mainMapper.reads.seqs);
        uint32_t numFilteredReads = disOptions.origReadIdMap[disOptions.currentBinNo].size();

        //load forward reads
        for (uint32_t i = 0; i< numFilteredReads; ++i)
        {
            uint32_t orgId = disOptions.origReadIdMap[disOptions.currentBinNo][i];
            appendValue(me.reads.seqs, mainMapper.reads.seqs[orgId]);
        }

        // if paired classify only one pair
        if (IsSameType<typename TMainConfig::TSequencing, PairedEnd>::VALUE)
        {
            uint32_t numPairs = getPairsCount( mainMapper.reads.seqs);
            //load mates
            for (uint32_t i = 0; i< numFilteredReads; ++i)
            {
                uint32_t orgId = disOptions.origReadIdMap[disOptions.currentBinNo][i];
                appendValue(me.reads.seqs, mainMapper.reads.seqs[orgId + numPairs]);
                disOptions.origReadIdMap[disOptions.currentBinNo].push_back(orgId + numPairs);
            }
            numFilteredReads *= 2; //now we have twice the reads
        }

        //load reverse reads
        for (uint32_t i = 0; i< numFilteredReads; ++i)
        {
            uint32_t orgId = disOptions.origReadIdMap[disOptions.currentBinNo][i];
            appendValue(me.reads.seqs, mainMapper.reads.seqs[orgId + numReads]);
            disOptions.origReadIdMap[disOptions.currentBinNo].push_back(orgId + numReads);
        }
    }
    stop(mainMapper.timer);
    disOptions.copyReads += getValue(mainMapper.timer);
    disOptions.filteredReads += getReadsCount(me.reads.seqs);
}

// ----------------------------------------------------------------------------
// Function mapReads()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void mapReads(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig>  & mainMapper, auto & mybiIndex, DisOptions & disOptions)
{
    _mapReadsImpl(me, mainMapper, mybiIndex, me.reads.seqs, disOptions);

//     _mapReadsImplOSS(me, mainMapper, me.biIndex, me.reads.seqs, disOptions);
}




template <typename TSpec, typename TConfig, typename TMainConfig>
inline void runMapper(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig> & mainMapper, DisOptions & disOptions)
{
//     typedef String<TChar, TAllocConfig> TString;
//     typedef StringSet<TString, Owner<ConcatDirect<SizeSpec_<TSeqNo, TSeqPos> > > > TStringSet;
//
//     using TFMIndexConfig = TGemMapFastFMIndexConfig<TBWTLen>;
//     TFMIndexConfig::SAMPLING = opt.sampling;



    typedef MapperTraits<TSpec, TConfig>                        TTraits;
//     typedef typename TTraits::TIndexConfig                      TIndexConfig;
    using TIndexConfig = YaraFMConfig<typename TConfig::TContigsSize, typename TConfig::TContigsLen, typename TConfig::TContigsSum, typename TConfig::TAlloc>;

//     using TTraits::TIndexConfig                      TIndexConfig;
    TIndexConfig::SAMPLING = me.samplingRate;

    using TIndexSpec = FMIndex<void, TIndexConfig>;
//     typedef BidirectionalIndex<TIndexSpec>                          TBiIndexSpec;
//     typedef Index<typename TIndexConfig::Text, TIndexSpec>          TIndex;
    using TBiIndex = Index<typename TIndexConfig::Text, BidirectionalIndex<TIndexSpec> >;

    TBiIndex mybiIndex;

    loadFilteredReads(me, mainMapper, disOptions);
    std::cout << "loaded Reads"<< "\n" << length(me.reads.seqs) << "\n";
    if (empty(me.reads.seqs)) return;
    loadContigs(me);
    std::cout << "loaded Contigs" << "\n" << length(me.contigs.seqs)<< "\n";
//     loadContigsIndex(me);
    loadContigsBiIndex(me, mybiIndex);
    std::cout << "loaded Index" << "\nIndexSize: ";
//     std::cout << length(me.biIndex.fwd.sa) << "\n";
    mapReads(me, mainMapper, mybiIndex, disOptions);
    std::cout << "Mapped Reads" << "\n";
}


// ----------------------------------------------------------------------------
// Function spawnMapper()
// ----------------------------------------------------------------------------
template <typename TContigsSize,
          typename TContigsLen,
          typename TContigsSum,
          typename TSpec,
          typename TMainConfig,
          typename TThreading,
          typename TSequencing,
          typename TSeedsDistance>
void spawnMapper(Options const & options,
                 Mapper<TSpec, TMainConfig> & mainMapper,
                 DisOptions & disOptions,
                 TThreading const & /*threading*/,
                 TSequencing const & /*sequencing*/,
                 TSeedsDistance const & /*distance*/)
{

    uint32_t numFilteredReads = (disOptions.filterType == NONE) ? getReadSeqsCount(mainMapper.reads.seqs) : disOptions.origReadIdMap[disOptions.currentBinNo].size();
    if(disOptions.mmap || numFilteredReads < disOptions.allocThreshold)
    {
        if(disOptions.verbose > 1)
            std::cout << "Memory Mapping\n";
        typedef ReadMapperConfig<TThreading, TSequencing, TSeedsDistance, TContigsSize, TContigsLen, TContigsSum>  TConfig;
        Mapper<void, TConfig> mapper(options);
        runMapper(mapper, mainMapper, disOptions);
    }
    else
    {
        if(disOptions.verbose > 1)
            std::cout << "Memory Allocation\n";
        typedef ReadMapperConfig<TThreading, TSequencing, TSeedsDistance, TContigsSize, TContigsLen, TContigsSum, Alloc<>>  TConfig;

        Mapper<void, TConfig> mapper(options);
        mapper.samplingRate = options.samplingRate;
        runMapper(mapper, mainMapper, disOptions);
    }
}
// ----------------------------------------------------------------------------
// Function configureMapper()
// ----------------------------------------------------------------------------
template <typename TContigsSize,
          typename TContigsLen,
          typename TSpec,
          typename TMainConfig,
          typename TThreading,
          typename TSequencing,
          typename TSeedsDistance>
void configureMapper(Options const & options,
                     Mapper<TSpec, TMainConfig> & mainMapper,
                     DisOptions & disOptions,
                     TThreading const & threading,
                     TSequencing const & sequencing,
                     TSeedsDistance const & distance)
{
    if (options.contigsSum <= MaxValue<uint32_t>::VALUE)
    {
        spawnMapper<TContigsSize, TContigsLen, uint32_t>(options, mainMapper, disOptions, threading, sequencing, distance);
    }
    else
    {
        throw RuntimeError("Maximum number of contigs exceeded. Recompile with -DYARA_LARGE_CONTIGS=ON dMapper.h.");
//         spawnMapper<TContigsSize, TContigsLen, uint64_t>(options, mainMapper, disOptions, threading, sequencing, distance);
    }
}

template <typename TContigsSize,
          typename TSpec,
          typename TMainConfig,
          typename TThreading,
          typename TSequencing,
          typename TSeedsDistance>
void configureMapper(Options const & options,
                     Mapper<TSpec, TMainConfig> & mainMapper,
                     DisOptions & disOptions,
                     TThreading const & threading,
                     TSequencing const & sequencing,
                     TSeedsDistance const & distance)
{
    if (options.contigsMaxLength <= MaxValue<uint32_t>::VALUE)
    {
        configureMapper<TContigsSize, uint32_t>(options, mainMapper, disOptions, threading, sequencing, distance);
    }
    else
    {
#ifdef DR_YARA_LARGE_CONTIGS
        configureMapper<TContigsSize, uint64_t>(options, mainMapper, disOptions, threading, sequencing, distance);
#else
        throw RuntimeError("Maximum contig length exceeded. Recompile with -DDR_YARA_LARGE_CONTIGS=ON.");
#endif
    }
}

template <typename TSpec, typename TMainConfig>
void configureMapper(Options const & options,
                     Mapper<TSpec, TMainConfig> & mainMapper,
                     DisOptions & disOptions)
{
    typedef typename MapperTraits<TSpec, TMainConfig>::TThreading       TThreading;
    typedef typename MapperTraits<TSpec, TMainConfig>::TSequencing      TSequencing;
    typedef typename MapperTraits<TSpec, TMainConfig>::TSeedsDistance   TSeedsDistance;

    if (options.contigsSize <= MaxValue<uint8_t>::VALUE)
    {
        configureMapper<uint8_t>(options, mainMapper, disOptions, TThreading(), TSequencing(), TSeedsDistance());
    }
    else if (options.contigsSize <= MaxValue<uint16_t>::VALUE)
    {
        configureMapper<uint16_t>(options, mainMapper, disOptions, TThreading(), TSequencing(), TSeedsDistance());
    }
    else
    {
#ifdef DR_YARA_LARGE_CONTIGS
        configureMapper<uint32_t>(options, mainMapper, disOptions, TThreading(), TSequencing(), TSeedsDistance());
#else
        throw RuntimeError("Maximum number of contigs exceeded. Recompile with -DYARA_LARGE_CONTIGS=ON.");
#endif
    }
}


template <typename TSpec = void>
struct AllContigs
{
    typedef SeqStore<TSpec, YaraContigsConfig<> >   TContigs;

    TContigs            contigs;
};


// ----------------------------------------------------------------------------
// Function loadAllContigs()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void loadAllContigs(Mapper<TSpec, TConfig> & mainMapper, DisOptions & disOptions)
{
    typedef typename MapperTraits<TSpec, TConfig>::TContigs          TContigs;

    CharString allContigsFile = disOptions.IndicesDirectory;
    allContigsFile += "allContigs";
//     std::cout << "AllContigsFile name " << allContigsFile << "\n";
    TContigs tmpContigs;
    if (!open(mainMapper.contigs, toCString(allContigsFile), OPEN_RDONLY)){
        start(mainMapper.timer);
//         std::cout << "Load with Alloc\n";
        try
        {
            AllContigs<> allContigs;
            for (uint32_t i=0; i < disOptions.numberOfBins; ++i)
            {

                if(disOptions.verbose > 1)
                    std::cout << "Load contig from " << i << "bin.\n";
                TContigs tmpContigs;
                CharString fileName;
                appendFileName(fileName, disOptions.IndicesDirectory, i);

                if (!open(tmpContigs, toCString(fileName), OPEN_RDONLY))
                    throw RuntimeError("Error while opening reference file.");
//                 append(mainMapper.contigs.seqs, tmpContigs.seqs);
//                 append(mainMapper.contigs.names, tmpContigs.names);
                append(allContigs.contigs.seqs, tmpContigs.seqs);
                append(allContigs.contigs.names, tmpContigs.names);
            }

            if (!save(/*mainMapper.contigs*/allContigs.contigs, toCString(allContigsFile)))
                throw RuntimeError("Error while saving all Contig references.");
            append(mainMapper.contigs.seqs, std::move(allContigs.contigs.seqs));
            append(mainMapper.contigs.names, std::move(allContigs.contigs.names));
        }
        catch (BadAlloc const & /* e */)
        {
            throw RuntimeError("Insufficient memory to load the reference.");
        }
    }
    else
    {
//         std::cout << "From Object\n";
    }
    stop(mainMapper.timer);
    mainMapper.stats.loadContigs += getValue(mainMapper.timer);

    if (mainMapper.options.verbose > 1)
        std::cout << "Loading reference:\t\t\t" << mainMapper.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function rankMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void rankMatches2(Mapper<TSpec, TConfig> & me, TReadSeqs const & readSeqs)
{
    typedef MapperTraits<TSpec, TConfig>                    TTraits;
    typedef typename TTraits::TMatch                        TMatch;
    typedef typename TTraits::TMatchesSet                   TMatchesSet;
    typedef typename TTraits::TMatchesViewSet               TMatchesViewSet;
    typedef typename Value<TMatchesSet const>::Type         TMatchesSetValue;
    typedef typename Value<TMatchesViewSet const>::Type     TMatchesViewSetValue;
    typedef typename Iterator<TMatchesViewSet const, Standard>::Type TMatchesViewSetIt;
    typedef typename Size<TReadSeqs>::Type                  TReadId;
    typedef typename Size<TMatchesSetValue>::Type           TMatchesSize;
    typedef std::uniform_int_distribution<TMatchesSize>     TMatchesRnd;

    start(me.timer);
    // Create a position modifier of the matches from the identity permutation.
    assign(me.matchesPositions, seqan::Range<TMatchesSize>(0, length(me.matchesByCoord)), Exact());
    setHost(me.matchesByErrors, me.matchesByCoord);
    setCargo(me.matchesByErrors, me.matchesPositions);

    // Bucket matches in the position modifier.
    setHost(me.matchesSetByErrors, me.matchesByErrors);
    assign(stringSetLimits(me.matchesSetByErrors), stringSetLimits(me.matchesSetByCoord), Exact());
    assign(stringSetPositions(me.matchesSetByErrors), stringSetPositions(me.matchesSetByCoord), Exact());

    // Sort matches by pairing info. iff possible

    // Sort matches by errors.
    forEach(me.matchesSetByErrors, sortMatches<TMatchesViewSetValue, Errors>, typename TTraits::TThreading());

    // Select all co-optimal matches.
    assign(me.optimalMatchesSet, me.matchesSetByErrors);
    clipMatches(me.optimalMatchesSet, countMatchesInBestStratum<TMatchesViewSetValue>, typename TTraits::TThreading());

    // Select all sub-optimal matches.
    assign(me.matchesSet, me.matchesSetByErrors);
    clipMatches(me.matchesSet, [&](TMatchesViewSetValue const & matches)
                {
                    if (empty(matches)) return TMatchesSize(0);

                    TReadId readId = getMember(front(matches), ReadId());

                    return countMatchesInStrata(matches, getReadStrata<TMatch>(me.options, length(readSeqs[readId])));
                },
                typename TTraits::TThreading());

    // Append an invalid match to matches by coord.
    resize(me.matchesByCoord, length(me.matchesByCoord) + 1, Exact());
    setInvalid(back(me.matchesByCoord));
    // Update matches by errors.
    resize(me.matchesPositions, length(me.matchesPositions) + 1, Exact());
    setPosition(me.matchesByErrors, length(me.matchesByErrors) - 1, length(me.matchesByCoord) - 1);

    // Initialize primary matches.
    setHost(me.primaryMatches, me.matchesByErrors);
    assign(me.primaryMatchesPositions, stringSetPositions(me.matchesSetByErrors), Exact());
    setCargo(me.primaryMatches, me.primaryMatchesPositions);

    // Choose primary matches among best matches.
    iterate(me.optimalMatchesSet, [&](TMatchesViewSetIt const & matchesIt)
            {
                // Use one generator per thread.
                std::default_random_engine generator;

                TReadId readId = position(matchesIt, me.optimalMatchesSet);
                TMatchesViewSetValue const & matches = value(matchesIt);

                // Set unmapped reads as invalid.
                if (empty(matches))
                {
                    setPosition(me.primaryMatches, readId, length(me.matchesByErrors) - 1);
                }
                // Choose match at random.
                else
                {
                    TMatchesRnd rnd(0, length(matches) - 1);
                    setPosition(me.primaryMatches, readId, position(me.primaryMatches, readId) + rnd(generator));
                }
            },
            Standard(), typename TTraits::TThreading());

    stop(me.timer);
    me.stats.sortMatches += getValue(me.timer);
    if (me.options.verbose > 1)
        std::cout << "Sorting time:\t\t\t" << me.timer << std::endl;

    // Update mapped reads.
    transform(me.ctx.mapped, me.primaryMatches, isValid<typename TTraits::TMatchSpec>, typename TTraits::TThreading());

    if (me.options.verbose > 0)
    {
        unsigned long mappedReads = count(me.ctx.mapped, true, typename TTraits::TThreading());
        me.stats.mappedReads += mappedReads;

        if (me.options.verbose > 1)
            std::cout << "Mapped reads:\t\t\t" << mappedReads << std::endl;
    }

    if (IsSameType<typename TConfig::TSequencing, SingleEnd>::VALUE) return;

    // Update paired reads.
    if (me.options.verbose > 0)
    {
        unsigned long pairedReads = count(me.ctx.paired, true, typename TTraits::TThreading());
        me.stats.pairedReads += pairedReads;

        if (me.options.verbose > 1)
        {
            std::cerr << "Pairing time:\t\t\t" << me.timer << std::endl;
            std::cerr << "Paired reads:\t\t\t" << pairedReads << std::endl;
        }
    }
}

// ----------------------------------------------------------------------------
// Function openOutputFile()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig>
inline void openOutputFile(Mapper<TSpec, TConfig> & mainMapper, DisOptions & disOptions)
{
    typedef MapperTraits<TSpec, TConfig>            TTraits;
    typedef typename TTraits::TContigs              TContigs;
    typedef typename TTraits::TContigSeqs           TContigSeqs;
    typedef typename Value<TContigSeqs>::Type       TContigSeq;

    String<uint32_t> allContigLengths;

    start(mainMapper.timer);
    try
    {
        for (uint32_t i=0; i < disOptions.numberOfBins; ++i)
        {
            TContigs tmpContigs;
            String<uint32_t> tmpContigLengths;
            CharString fileName;
            appendFileName(fileName, disOptions.IndicesDirectory, i);

            if (!open(tmpContigs, toCString(fileName), OPEN_RDONLY))
                throw RuntimeError("Error while opening reference file.");

            resize(tmpContigLengths, length(tmpContigs.seqs));
            transform(tmpContigLengths, tmpContigs.seqs, [](TContigSeq const & seq) { return length(seq); });
            append(allContigLengths, tmpContigLengths);

            append(mainMapper.contigs.names, tmpContigs.names);
        }
    }
    catch (BadAlloc const & /* e */)
    {
        throw RuntimeError("Insufficient memory to load the reference.");
    }
    stop(mainMapper.timer);
    mainMapper.stats.loadContigs += getValue(mainMapper.timer);

    if (mainMapper.options.verbose > 1)
        std::cout << "Loading reference:\t\t\t" << mainMapper.timer << std::endl;

    bool opened = false;

    if (empty(mainMapper.options.outputFile))
    {
        // Output to cout.
        if (mainMapper.options.uncompressedBam)
        {
            // Turn off BAM compression.
            setFormat(mainMapper.outputFile, mainMapper.options.outputFormat);
            opened = _open(mainMapper.outputFile, std::cout, Nothing(), False());
        }
        else
        {
            opened = open(mainMapper.outputFile, std::cout, mainMapper.options.outputFormat);
        }
    }
    else
    {
        // Output to file.
        opened = open(mainMapper.outputFile, toCString(mainMapper.options.outputFile), OPEN_WRONLY | OPEN_CREATE);
    }

    if (!opened) throw RuntimeError("Error while opening output file.");

    setContigNames(context(mainMapper.outputFile), mainMapper.contigs.names);

    // Fill contig lengths.
    resize(contigLengths(context(mainMapper.outputFile)), length(allContigLengths));
    assign(contigLengths(context(mainMapper.outputFile)), allContigLengths);

    typedef FileFormat<BamFileOut>::Type    TOutputFormat;
    TOutputFormat of;
    assign(of, Bam());

    if(mainMapper.outputFile.format.tagId == of.tagId || !disOptions.skipSamHeader)
    {
        // Write header.
        BamHeader header;
        fillHeader(header, mainMapper.options);
        writeHeader(mainMapper.outputFile, header);
    }
}

// ----------------------------------------------------------------------------
// Function prepairMainMapper()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TMainConfig, typename TFilter>
inline void prepairMainMapper(Mapper<TSpec, TMainConfig> & mainMapper, TFilter const & filter, DisOptions & disOptions)
{
    initReadsContext(mainMapper, mainMapper.reads.seqs);
    setHost(mainMapper.cigars, mainMapper.cigarString);
    disOptions.collectedCigars.clear();
    if (IsSameType<typename TMainConfig::TSequencing, PairedEnd>::VALUE)
        resize(mainMapper.primaryMatchesProbs, getReadsCount(mainMapper.reads.seqs), 0.0, Exact());
    if(disOptions.filterType != NONE)
        clasifyLoadedReads(mainMapper, filter, disOptions);
}

// ----------------------------------------------------------------------------
// Function finalizeMainMapper()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TMainConfig>
inline void finalizeMainMapper(Mapper<TSpec, TMainConfig> & mainMapper, DisOptions & disOptions)
{
    /*
    std::cout << "Print Matches\n";
    for(int i = 0; i < 1; ++i){
            auto matchIt = begin(matches, Standard());
            auto matchEnd = end(matches, Standard());
            while(matchIt != matchEnd){
                write(std::cout, *matchIt);
                ++matchIt;
            }
    }*/
    aggregateMatches(mainMapper, mainMapper.reads.seqs);
/*
//     std::cout << "Print Filtered Matches\n";
    for(int i = 0; i < length(mainMapper.matchesSetByCoord); ++i){
            auto const & matches = mainMapper.matchesSetByCoord[i];
            auto matchIt = begin(matches, Standard());
            auto matchEnd = end(matches, Standard());
            while(matchIt != matchEnd){
                  if (getMember(*matchIt, ContigEnd()) - getMember(*matchIt, ContigBegin()) >= 200){
                    std::cout << "After sort Match\n";
                    write(std::cout, *matchIt);
                }
//                 write(std::cout, *matchIt);
                ++matchIt;
            }
    }*/




    rankMatches2(mainMapper, mainMapper.reads.seqs);

    /*
    for(int i = 0; i < length(mainMapper.matchesSetByCoord); ++i){
        auto const & matches = mainMapper.matchesSetByCoord[i];
            auto matchIt = begin(matches, Standard());
            auto matchEnd = end(matches, Standard());
            while(matchIt != matchEnd){
                  if (getMember(*matchIt, ContigEnd()) - getMember(*matchIt, ContigBegin()) >= 200){
                    std::cout << "After rank Match\n";
                    write(std::cout, *matchIt);
                }
//                 write(std::cout, *matchIt);
                ++matchIt;
            }
    }*/

//     transferCigars(mainMapper, disOptions);
    std::cout << "Load All Contigs\n";
    loadAllContigs(mainMapper, disOptions);
    std::cout << "Align all matches\n";
    alignMatches(mainMapper);

    writeMatches(mainMapper);
    clearMatches(mainMapper);
    clearAlignments(mainMapper);
    clearReads(mainMapper);
}

// ----------------------------------------------------------------------------
// Function sortedBins()
// ----------------------------------------------------------------------------
std::vector<uint32_t> sortedBins(DisOptions const & disOptions)
{
    std::vector<uint32_t> sortedBinIndex(disOptions.numberOfBins);
    iota(sortedBinIndex.begin(), sortedBinIndex.end(), 0);

    // sort indexes based on comparing values in v
    std::sort(sortedBinIndex.begin(), sortedBinIndex.end(),
         [&disOptions](size_t i1, size_t i2) {return disOptions.origReadIdMap[i1].size() > disOptions.origReadIdMap[i2].size();});

    return sortedBinIndex;
}

// ----------------------------------------------------------------------------
// Function runDisMapper()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TMainConfig, typename TFilter>
inline void runDisMapper(Mapper<TSpec, TMainConfig> & mainMapper, TFilter const & filter, DisOptions & disOptions)
{
    configureThreads(mainMapper);

    // Open output file and write header.
    openOutputFile(mainMapper, disOptions);
    openReads(mainMapper);

    while (true)
    {
        if (mainMapper.options.verbose > 1) printRuler(std::cout);
        loadReads(mainMapper);
        if (empty(mainMapper.reads.seqs)) break;
        prepairMainMapper(mainMapper, filter, disOptions);

        if(disOptions.filterType == NONE){
            for(uint32_t i = disOptions.startBin; i < disOptions.numberOfBins; ++i){
                //TODO iterate over reads batches using one Index
                std::cout << "\n\nIn bin Number: " << i << "\n";
                disOptions.currentBinNo = i;
                Options options = mainMapper.options;
                appendFileName(options.contigsIndexFile, disOptions.IndicesDirectory, i);
                appendFileName(options.mappabilitySubDirectory, disOptions.MappabilityDirectory, i);
                // also load SamplingRate;
                if (!openContigsLimits(options))
                    throw RuntimeError("Error while opening reference file.");
                configureMapper<TSpec, TMainConfig>(options, mainMapper, disOptions);
            }
        }
        else
        {
            for (auto i: sortedBins(disOptions))
            {
                disOptions.currentBinNo = i;
                Options options = mainMapper.options;
                appendFileName(options.contigsIndexFile, disOptions.IndicesDirectory, i);
                appendFileName(options.mappabilitySubDirectory, disOptions.MappabilityDirectory, i);
                if (!openContigsLimits(options))
                    throw RuntimeError("Error while opening reference file.");
                configureMapper<TSpec, TMainConfig>(options, mainMapper, disOptions);
            }
        }

        finalizeMainMapper(mainMapper, disOptions);
    }
    closeReads(mainMapper);
    closeOutputFile(mainMapper);
}

// ----------------------------------------------------------------------------
// Function spawnDisMapper()
// ----------------------------------------------------------------------------
template <typename TContigsSize, typename TContigsLen, typename TContigsSum,
typename TThreading, typename TSequencing, typename TSeedsDistance>
inline void spawnDisMapper(DisOptions & disOptions,
                           TThreading const & /* tag */,
                           TSequencing const & /* tag */,
                           TSeedsDistance const & /* tag */)
{
    disOptions.outputFile = disOptions.superOutputFile;
    typedef ReadMapperConfig<TThreading, TSequencing, TSeedsDistance, TContigsSize, TContigsLen, TContigsSum, MMap<>>  TMainConfig;
    Mapper<void, TMainConfig> disMapper(disOptions);

    Timer<double> timer;

    start(timer);
    start(disMapper.timer);

    if (disOptions.filterType == BLOOM)
    {
        SeqAnBloomFilter<> filter  (toCString(disOptions.filterFile));

        disOptions.kmerSize = filter.getKmerSize();

//        if(filter.getNumberOfBins() != disOptions.numberOfBins)
//            std::cerr << "[WARNING] Provided number of bins (" << disOptions.numberOfBins << ")differs from that of the bloom filter (" << filter.getNumberOfBins() << ")";

        stop(disMapper.timer);
        disOptions.loadFilter += getValue(disMapper.timer);
        runDisMapper(disMapper, filter, disOptions);
    }
    else if (disOptions.filterType == KMER_DIRECT)
    {
        SeqAnKDXFilter<> filter (toCString(disOptions.filterFile));

        disOptions.kmerSize = filter.getKmerSize();

//        if(filter.getNumberOfBins() != disOptions.numberOfBins)
//            std::cerr << "[WARNING] Provided number of bins (" << disOptions.numberOfBins << ")differs from that of the bloom filter (" << filter.getNumberOfBins() << ")\n";

        stop(disMapper.timer);
        disOptions.loadFilter += getValue(disMapper.timer);
        runDisMapper(disMapper, filter, disOptions);
    }
    else
    {
        // dummy filter in case of nofilter option
        SeqAnBloomFilter<> filter(64, 3, 20, 1);

        stop(disMapper.timer);
        disOptions.loadFilter += getValue(disMapper.timer);
        runDisMapper(disMapper, filter, disOptions);
    }
    stop(timer);
    if (disMapper.options.verbose > 0)
    {
        double total = getValue(timer) / 100.0;

        std::cout << "\nFilter loading time:\t\t" << disOptions.loadFilter << " sec" << "\t\t" << disOptions.loadFilter / total << " %" << std::endl;
        std::cout << "Reads filtering time:\t\t" << disOptions.filterReads << " sec" << "\t\t" << disOptions.filterReads / total << " %" << std::endl;
        std::cout << "Reads copying time:\t\t" << disOptions.copyReads << " sec" << "\t\t" << disOptions.copyReads / total << " %" << std::endl;
        std::cout << "Alignments copying time:\t" << disOptions.copyAlignments << " sec" << "\t\t" << disOptions.copyAlignments / total << " %" << std::endl;
        std::cout << "Cigars moving time:\t\t" << disOptions.moveCigars << " sec" << "\t\t" << disOptions.moveCigars / total << " %" << std::endl;

        printStats(disMapper, timer);
		std::cout << "Avg reads per bin:\t\t" << (double)disOptions.filteredReads / disOptions.numberOfBins << std::endl;
	}
}

#endif  // #ifndef APP_YARA_MAPPER_H_
