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
    bool                    noMappability = false;
    uint32_t                startBin = 0;
    uint32_t                readLength = 0;

    double                  loadFilter      = 0.0;
    double                  filterReads     = 0.0;
    double                  copyReads       = 0.0;
    double                  copyAlignments  = 0.0;
    double                  moveCigars      = 0.0;

    bool                    skipSamHeader = false;

    uint32_t                kmerSize = 20;
    uint32_t                numberOfBins = 64;

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
    typedef MapperTraits<TSpec, TConfig>                    Traits;
    typedef typename Traits::TMatch                        TMatch;
    typedef typename Traits::TMatchesAppender              TMatches;


    TMatches &          matches;

    DelegateDirect(TMatches & matches) :
        matches(matches)
    {}

    template <typename TContext, typename TContigsPos, typename TReadId, typename TMatchErrors>
    void operator() (TContext & ossContext, TContigsPos const & start, TContigsPos const & end, TReadId const needleId, TMatchErrors errors)
    {
        TMatch hit;
        setContigPosition(hit, start, end);
        hit.errors = errors;
        setReadId(hit, ossContext.readSeqs, needleId);

        TReadId readId = getReadId(ossContext.readSeqs, needleId);
        setMapped(ossContext.ctx, readId);
        setMinErrors(ossContext.ctx, readId, errors);


//         std::cout << "Direct hit" << start << "end" << end << " (" << (int)getSeqOffset(end) - getSeqOffset(start) << ")" << "\t" << needleId << "\terrors" << (int) errors << "\n";
//         setMinErrors(ossContext.ctx, needleId, errors);
//         hit.readId = needleId;

//         THit hit = { range(indexIt), (TSeedId)position(seedsIt), errors };

        appendValue(matches, hit, Generous(), typename Traits::TThreading());
    }
};

template <typename TTraits>
struct Delegate
{
//     typedef MapperTraits<TSpec, TConfig>                   Traits;
    typedef typename TTraits::TMatch                        TMatch;
    typedef typename TTraits::TMatchesAppender              TMatches;
    typedef typename TTraits::TContigsPos                   TContigsPos;
//     typedef typename TTraits::TSA                           TSA;
//     typedef typename Size<TSA>::Type                        TSAPos;
//     typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename TTraits::TReadSeqs                     TReadSeqs;
    typedef typename Size<TReadSeqs>::Type                  TReadId;


    TMatches &          matches;

    Delegate(TMatches & matches) :
        matches(matches)
    {}

    template <typename TContext, typename TNeedleId, typename TMatchErrors>
    void operator() (TContext & ossContext, auto const & iter, TNeedleId const needleId, TMatchErrors const errors, bool const rev)
    {
        TReadId readId = getReadId(ossContext.readSeqs, needleId);
        uint32_t occLength = repLength(iter);
        for (TContigsPos occ : getOccurrences(iter)){
//         for (TSAPos i = iter.fwdIter.vDesc.range.i1; i < iter.fwdIter.vDesc.range.i2; ++i){
//             TSAValue saPos = iter.fwdIter.index->sa[i];
//              std::cout << occ << "\n";
            TMatch hit;
            setContigPosition(hit, occ, posAdd(occ, occLength));
            hit.errors = errors;
            setReadId(hit, ossContext.readSeqs, needleId); // needleId is used to determine if read is reverse complement

            setMapped(ossContext.ctx, readId);
            setMinErrors(ossContext.ctx, readId, errors); //TODO move out of loop

//             std::cout << "hit" << occ << "end" << posAdd(occ, occLength) << " (" << occLength << ")" << "\t" << needleId << "\terrors" << (int) errors << "\n";
//             std::cout << "isMapped: "<< isMapped(ossContext.ctx, readId) << "\n";
//             std::cout << "MinErrors: "<< (int)getMinErrors(ossContext.ctx, readId) << "\n";

    //         hit.readId = readId;
//TODO maybe use this form
    //         THit hit = { range(indexIt), (TSeedId)position(seedsIt), errors };

            appendValue(matches, hit, Generous(), typename TTraits::TThreading()); //does this make any sense (always single occ)
        }
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
    mainMapper.stats.sortMatches    += childMapper.stats.sortMatches;
    mainMapper.stats.compactMatches += childMapper.stats.compactMatches;
    mainMapper.stats.selectPairs    += childMapper.stats.selectPairs;
    mainMapper.stats.verifyMatches  += childMapper.stats.verifyMatches;
    mainMapper.stats.alignMatches   += childMapper.stats.alignMatches;
    mainMapper.stats.writeMatches   += childMapper.stats.writeMatches;
    mainMapper.stats.rescuedReads   += childMapper.stats.rescuedReads;
}


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



struct SHit{
    Pair <uint64_t, uint64_t> occ;
    template <typename TOcc>
    SHit(TOcc & inOcc):
        occ(inOcc)
    {}
};

template<size_t minErrors, size_t maxErrors,
         typename TIndex, typename TContigSeqs, typename TMatch>
int testReadOcc(TIndex & index, TContigSeqs & text, TMatch & match, uint8_t maErrors, int threshold) //add Distance //Distance Mappability
{
//     int mErrors = maErrors;

    std::vector<SHit> hits;

    auto delegate = [&hits](auto & iter, DnaString const & needle, uint8_t errors)
    {
        for (auto occ : getOccurrences(iter)){
            SHit me(occ);
            hits.push_back(me);
        }
    };

    StringSet<Dna5String> readOcc;

//      std::cout << getMember(*newIt, Errors()) << "\t < 5 , " << getMember(*newIt, ContigBegin()) << " >\t" << getMember(*newIt, ContigEnd()) << "\n";
/*
    if(false){ //edit
        Dna5String part = infix(text[getSeqNo(getMember(match, ContigId()))], getSeqOffset(getMember(match, ContigBegin())), getSeqOffset(getMember(match, ContigEnd())));
        appendValue(readOcc, part);
    }else{
//         bool lim = testhit.occ.i2 - mErrors > 0 && testhit.occ.i2 + seqan::length(testhit.read) + 2 * mErrors < length(genome[testhit.occ.i1]);
        if(getSeqOffset(getMember(match, ContigBegin()) < mErrors ||
            length(text[getSeqNo(getMember(match, ContigId()))]) < getSeqOffset(getMember(match, ContigEnd())) + mErrors)
            return(666);
        for(int off = -maxErrors; off <= mErrors; ++off){
            Dna5String part = infix(text[getSeqNo(getMember(match, ContigId()))],
                                    static_cast<int64_t>(getSeqOffset(getMember(match, ContigBegin()))) + off,
                                    static_cast<int64_t>(getSeqOffset(getMember(match, ContigEnd()))) + off);
            appendValue(readOcc, part);
        }
    }

    if(false){//edit
        find<minErrors, maxErrors>(delegate, index, readOcc[0], HammingDistance());
        std::cout << hits.size() << " hits!!!!!!!!!!" << "\n";
    }
    else
    {
        if(false)
            find<minErrors, maxErrors>(delegate, index, readOcc[mErrors], HammingDistance());
        else
            find<minErrors, maxErrors>(delegate, index, readOcc[mErrors], EditDistance());

        std::cout << hits.size() << " hits!!!!!!!!!!" << "\n";

        int k = 0;
        while(hits.size() < threshold && k < length(readOcc)){
            hits.clear();
            if(false)
                find<minErrors, maxErrors>(delegate, index, readOcc[k], HammingDistance());
            else
                find<minErrors, maxErrors>(delegate, index, readOcc[k], EditDistance());
            ++k;
            //use sort and erase
            std::cout << hits.size() << " hits!!!!!!!!!!e" << "\n";
        }
    }*/

    int nhits = hits.size();
    return nhits;
}





template<typename TIndex, typename TContigSeqs, typename TMatch>
int testReadOcc(TIndex & index, TContigSeqs & text, TMatch & match, uint8_t maxErrors, int threshold) //add Distance //Distance Mappability
{
    int hits = 0;

    switch (maxErrors)
    {
        case 1: hits = testReadOcc<0, 1>(index, text, match, maxErrors, threshold);//add Distance //Distance Mappability
                break;
        case 2: hits = testReadOcc<0, 2>(index, text, match, maxErrors, threshold);//add Distance //Distance Mappability
                break;
        case 3: hits = testReadOcc<0, 3>(index, text, match, maxErrors, threshold);//add Distance //Distance Mappability
                break;
        case 4: hits = testReadOcc<0, 4>(index, text, match, maxErrors, threshold);//add Distance //Distance Mappability
                break;
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                std::exit(1);
    }
    return hits;
}


// ----------------------------------------------------------------------------
// Function _mapReadsImpl()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig, typename TIndex, typename TBiIndex, typename TReadSeqs, typename TSeqsSpec>
inline void _mapReadsImpl(Mapper<TSpec, TConfig> & me,
                          Mapper<TSpec, TMainConfig>  & mainMapper,
                          TIndex & index,
                          TBiIndex & biindex,
                          StringSet<TReadSeqs, TSeqsSpec> & readSeqs,
                          DisOptions & disOptions)
{

    uint32_t len;
    if(disOptions.readLength != 0)
        len = disOptions.readLength;
    else
        len = length(readSeqs[0]);

    uint16_t maxError = me.options.errorRate * len;
    uint16_t strata = disOptions.strataRate * len;

    if(!disOptions.ossOff){

    initReadsContext(me, readSeqs);

    typedef MapperTraits<TSpec, TConfig>                        TTraits;
    typedef typename TTraits::TMatchesAppender                  TMatchesAppender;
    typedef Delegate<TTraits>                                   Delegate;
    typedef DelegateDirect<TSpec, TConfig>                      DelegateDirect;

    typedef typename TTraits::TContigSeqs                       TContigSeqs;


//     reserve(me.matchesByCoord, countHits(me) / 3); //some formula includding distance errors and mappability

    TMatchesAppender appender(me.matchesByCoord);
    Delegate delegate(appender);
    DelegateDirect delegateDirect(appender);
    TContigSeqs & contigSeqs = me.contigs.seqs;

    OSSContext<TSpec, TConfig> ossContext(me.ctx, appender, readSeqs, contigSeqs);

    if(disOptions.verbose > 1){
        std::cout << "maxError: " << maxError << "\t" << strata << "\n";
    }

    if(!disOptions.noMappability){
        start(me.timer);
        CharString bPath = me.options.mappabilitySubDirectory;
        bPath += "/";
        if(disOptions.verbose > 1)
            std::cout << "\nLoading Bitvectors: " << bPath << "\n";
        loadAllBitvectors(bPath, me.bitvectors, me.bitvectorsMeta, len);
        std::cout << "Bit vectors loaded. Number: " << me.bitvectors.size() << "\n";

        if(!me.bitvectors.empty()){
            std::cout << "Index Size: " << length(me.biIndex.fwd.sa) << "\n";
            std::cout << "Number of Sequences: " << length(me.contigs.seqs) << "\n";

            std::cout << "Length of Bitvectors: " << me.bitvectors[0].first.size() << "\n";
            ossContext.bitvectorsMeta = me.bitvectorsMeta;
        }
        stop(me.timer);
        me.stats.loadingBitvectors += getValue(me.timer);
        if(disOptions.verbose > 1)
            std::cerr << "Loading Bitvectors time:\t\t\t" << me.timer << std::endl;
    }

    // copy parameters to ossContext
    bool mscheme = true;
    ossContext.setReadContextOSS(maxError, strata, mscheme);
    ossContext.readLength = len;
    ossContext.numberOfSequences = length(me.contigs.seqs);

    start(me.timer);
    if(mscheme){
        find(0, maxError, strata, ossContext, delegate, delegateDirect, me.biIndex, me.bitvectors, readSeqs, EditDistance());
    }else{
        find(0, maxError, ossContext, delegate, delegateDirect, me.biIndex, me.bitvectors, readSeqs, EditDistance());
    }
    stop(me.timer);
    me.stats.optimumSearch += getValue(me.timer);
    if(disOptions.verbose > 0)
    {
        std::cerr << "\nOptimum Search time:\t\t" << me.timer << std::endl;
        std::cerr << "Matches count:\t\t\t" << lengthSum(me.matchesByCoord) << std::endl;
    }

//     if(addMyOwnOption) //IsSameType<typename TConfig::TSeedsDistance, EditDistance>::VALUE
//         find(0, maxError, ossContext, delegate, delegateDirect, me.biIndex, readSeqs, EditDistance());
//     else
//         find(0, maxError, ossContext, delegate, delegateDirect, me.biIndex, readSeqs, HammingDistance());


//     std::cout << "Finished Searching: " << "\n";
//     typedef typename TTraits::TMatch                             TMatch;
//     for(int i = 0; i < length(me.matchesByCoord); ++i){
//         std::cout << "Match: " << "\n";
//         TMatch myMatch = me.matchesByCoord[i];
//         std::cout << myMatch.contigBegin << "\n";
//         std::cout << myMatch.contigEnd << "\n";
//         std::cout << "ContigId: " << myMatch.contigId << "\n";
//         std::cout << "Errors: " << myMatch.errors << "\n";
//         std::cout << "ReadId: "<< myMatch.readId << "\n" << "\n";
//     }
//     std::cout << "End of Matches ____________________________________________________________________________________________________ " << "\n";

    aggregateMatchesOSS(me, readSeqs);


    if(disOptions.verbose > 0)
    {
        std::cerr << "Unique Matches count:\t\t\t" << lengthSum(me.matchesSetByCoord)/*length(me.matchesByCoord)*/ << std::endl;
    }

    }
/*
    auto sets = me.matchesSetByCoord;
    auto matches = me.matchesByCoord;

    clear(me.matchesSetByCoord);
    clear(me.matchesByCoord);*/


    if(disOptions.ossOff/* || true*/)
    {
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

    }



/*
    for(int i = 0; i < length(sets); ++i){
        auto const & matches = sets[i];
        auto matchIt = begin(matches, Standard());
        auto matchEnd = end(matches, Standard());
        while(matchIt != matchEnd){
            write(std::cout, *matchIt);
            ++matchIt;
        }
    }

    std::cout << "Print OSS Matches" << "\n\n\n";

    for(int i = 0; i < length(me.matchesSetByCoord); ++i){
        auto const & matches = me.matchesSetByCoord[i];
        auto matchIt = begin(matches, Standard());
        auto matchEnd = end(matches, Standard());
        while(matchIt != matchEnd){
            write(std::cout, *matchIt);
            ++matchIt;
        }
    }


    std::cout << "\n\n";
    int offSet = 0;
    auto matchesOSS = me.matchesSetByCoord[0];
    auto matchItOSS = begin(matchesOSS, Standard());
    auto matchEndOSS = end(matchesOSS, Standard());
    bool wrong = false;
    for(int i = 0; i < length(me.matchesSetByCoord); ++i){
        auto const matches = sets[i];

        auto matchIt = begin(matches, Standard());
        auto matchEnd = end(matches, Standard());


        while(matchIt != matchEnd){
            if(wrong)
                std::cout << "Something went wrong\n";

            bool same = isDuplicate(*matchIt, *matchItOSS, ContigBegin());
            if(!same){

                std::cout << "Need to verify this: " << "\n";
                write(std::cout, *matchIt);
                std::cout << "Compared to this one: " << "\n";
                write(std::cout, *matchItOSS);
                int nhits = testReadOcc( me.biIndex, me.contigs.seqs, *matchIt, maxError, 10); //TODO add Threshold as input option for me
                if(nhits > 10)
                    wrong = true;

                //start search occ
                //print results
                ++matchIt;
            }
            else
            {
                ++matchIt;
                ++matchItOSS;

                if(matchEndOSS == matchItOSS)
                {

                    ++offSet;
                    if(offSet < length(me.matchesSetByCoord))
                        --offSet;
                    matchesOSS = me.matchesSetByCoord[offSet];
                    matchItOSS = begin(matchesOSS, Standard());
                    matchEndOSS = end(matchesOSS, Standard());
                }
            }
        }
    }*/



    rankMatches(me, me.reads.seqs);
    if (me.options.verifyMatches)
        verifyMatches(me);
    alignMatches(me);
    copyMatches(mainMapper, me, disOptions);
    copyCigars(mainMapper, me, disOptions);
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


    typedef Delegate<TTraits>          TDelegate;
    TDelegate delegate(appender);

    delegate(ossContext, start, end, errors, readID);
    std::cout << "experimental Delegate call: "<< length(me.matchesByCoord) << "\n";

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
inline void mapReads(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig>  & mainMapper, DisOptions & disOptions)
{
    _mapReadsImpl(me, mainMapper, me.index, me.biIndex, me.reads.seqs, disOptions);

//     _mapReadsImplOSS(me, mainMapper, me.biIndex, me.reads.seqs, disOptions);
}




template <typename TSpec, typename TConfig, typename TMainConfig>
inline void runMapper(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig> & mainMapper, DisOptions & disOptions)
{
    // add load bitvectors here also add info to disoption about them
    loadFilteredReads(me, mainMapper, disOptions);
    if (empty(me.reads.seqs)) return;
    loadContigs(me);
//     loadContigsIndex(me);
    loadContigsBiIndex(me);
    mapReads(me, mainMapper, disOptions);
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
    typedef ReadMapperConfig<TThreading, TSequencing, TSeedsDistance, TContigsSize, TContigsLen, TContigsSum>  TConfig;
    Mapper<void, TConfig> mapper(options);
    runMapper(mapper, mainMapper, disOptions);
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
        spawnMapper<TContigsSize, TContigsLen, uint64_t>(options, mainMapper, disOptions, threading, sequencing, distance);
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


// ----------------------------------------------------------------------------
// Function loadAllContigs()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void loadAllContigs(Mapper<TSpec, TConfig> & mainMapper, DisOptions & disOptions)
{
    typedef typename MapperTraits<TSpec, TConfig>::TContigs          TContigs;

    start(mainMapper.timer);
    try
    {
        for (uint32_t i=0; i < disOptions.numberOfBins; ++i)
        {
            TContigs tmpContigs;
            CharString fileName;
            appendFileName(fileName, disOptions.IndicesDirectory, i);

            if (!open(tmpContigs, toCString(fileName), OPEN_RDONLY))
                throw RuntimeError("Error while opening reference file.");
            append(mainMapper.contigs.seqs, tmpContigs.seqs);
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
        std::cerr << "Loading reference:\t\t\t" << mainMapper.timer << std::endl;
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
        std::cerr << "Sorting time:\t\t\t" << me.timer << std::endl;

    // Update mapped reads.
    transform(me.ctx.mapped, me.primaryMatches, isValid<typename TTraits::TMatchSpec>, typename TTraits::TThreading());

    if (me.options.verbose > 0)
    {
        unsigned long mappedReads = count(me.ctx.mapped, true, typename TTraits::TThreading());
        me.stats.mappedReads += mappedReads;

        if (me.options.verbose > 1)
            std::cerr << "Mapped reads:\t\t\t" << mappedReads << std::endl;
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
        std::cerr << "Loading reference:\t\t\t" << mainMapper.timer << std::endl;

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
    if(disOptions.ossOff)
        aggregateMatches(mainMapper, mainMapper.reads.seqs);
    else
        aggregateMatchesOSS(mainMapper, mainMapper.reads.seqs);
    rankMatches2(mainMapper, mainMapper.reads.seqs);
    transferCigars(mainMapper, disOptions);

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
        if (mainMapper.options.verbose > 1) printRuler(std::cerr);
        loadReads(mainMapper);
        if (empty(mainMapper.reads.seqs)) break;
        prepairMainMapper(mainMapper, filter, disOptions);

        if(disOptions.filterType == NONE){
            for(uint32_t i = disOptions.startBin; i < disOptions.numberOfBins; ++i){ //weird interaction with InTextVerification
                std::cout << "\n\nIn bin Number: " << i << "\n";
                disOptions.currentBinNo = i;
                Options options = mainMapper.options;
                appendFileName(options.contigsIndexFile, disOptions.IndicesDirectory, i);
                appendFileName(options.mappabilitySubDirectory, disOptions.MappabilityDirectory, i);
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
    typedef ReadMapperConfig<TThreading, TSequencing, TSeedsDistance, TContigsSize, TContigsLen, TContigsSum>  TMainConfig;
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

        std::cerr << "\nFilter loading time:\t\t" << disOptions.loadFilter << " sec" << "\t\t" << disOptions.loadFilter / total << " %" << std::endl;
        std::cerr << "Reads filtering time:\t\t" << disOptions.filterReads << " sec" << "\t\t" << disOptions.filterReads / total << " %" << std::endl;
        std::cerr << "Reads copying time:\t\t" << disOptions.copyReads << " sec" << "\t\t" << disOptions.copyReads / total << " %" << std::endl;
        std::cerr << "Alignments copying time:\t" << disOptions.copyAlignments << " sec" << "\t\t" << disOptions.copyAlignments / total << " %" << std::endl;
        std::cerr << "Cigars moving time:\t\t" << disOptions.moveCigars << " sec" << "\t\t" << disOptions.moveCigars / total << " %" << std::endl;

        printStats(disMapper, timer);
		std::cerr << "Avg reads per bin:\t\t" << (double)disOptions.filteredReads / disOptions.numberOfBins << std::endl;
	}
}

#endif  // #ifndef APP_YARA_MAPPER_H_
