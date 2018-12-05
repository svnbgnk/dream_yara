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
    {
        ;
    }
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

    typedef MapperTraits<TSpec, TConfig>        TTraits;
    typedef typename TTraits::TReadsContext     TReadsContext;
    typedef typename TTraits::TMatchesAppender  TMatches;
    typedef typename TTraits::TMatch            TMatch;
    typedef typename TTraits::TContigSeqs       TContigSeqs;
    typedef typename TTraits::TReadSeqs         TReadSeqs;
    typedef typename TTraits::TBiIter           MySparseIter;
    typedef State<MySparseIter>                 TTState;


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
        clear(ctxOSS);
        resize(ctxOSS, readSeqs);

        if(!mScheme){
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

    template <size_t nbrBlocks>
    bool itvCondition(OptimalSearch<nbrBlocks> const & s,
                      uint8_t const blockIndex,
                      uint32_t ivalOne)
    {
        return(itv && ivalOne < (static_cast<int>(s.pi.size()) - blockIndex - 1 + normal.directsearchblockoffset) * normal.directsearch_th);
    }


    template <typename TText, typename TIndex, typename TIndexSpec,
              size_t nbrBlocks>
    bool itvConditionComp(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                      uint32_t const needleLeftPos,
                      uint32_t const needleRightPos,
                      uint8_t const errors,
                      OptimalSearch<nbrBlocks> const & s,
                      uint8_t const blockIndex)
    {
        return(itv && iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1 < (s.pi.size() - blockIndex - 1 + comp.directsearchblockoffset) * comp.directsearch_th);
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

    template <typename TContext, typename TContigsPos, typename TMatchErrors, typename TReadId>
    void operator() (TContext & ossContext, TContigsPos const & start, TContigsPos const & end, TMatchErrors errors, TReadId const needleId)
    {
        TMatch hit;
        setContigPosition(hit, start, end);
        hit.errors = errors;
        setReadId(hit, ossContext.readSeqs, needleId);
        setMapped(ossContext.ctx, needleId);
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
//             std::cout << "isMapped: "<< isMapped(ossContext.ctx, readId) << "\n";
//             std::cout << "MinErrors: "<< getMinErrors(ossContext.ctx, readId) << "\n";

            setReadId(hit, ossContext.readSeqs, needleId); // since reverse reads are at the end //TODO check this line (it determines if read  is rev) // if would take the original id than it can no longer determine if it is reversed or not
            setMapped(ossContext.ctx, readId);
            setMinErrors(ossContext.ctx, readId, errors);

//             std::cout << "isMapped: "<< isMapped(ossContext.ctx, readId) << "\n";
//             std::cout << "MinErrors: "<< (int)getMinErrors(ossContext.ctx, readId) << "\n";

    //         hit.readId = readId;
//TODO to du maybe use this form
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


template <typename TSpec, typename TConfig,
          typename TDelegatee, typename TDelegateDe,
          typename TIndex,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void myTfind(const int minErrors,
                    const int maxErrors,
                    OSSContext<TSpec, TConfig> & ossContext,
                    TDelegatee & delegate,
                    TDelegateDe & delegateDirect,
                    TIndex & index,
                    StringSet<TNeedle, TStringSetSpec> const & needles,
                    TDistanceTag const &)
{
    //TODO find out why i cant use lambda functions

//     switch (maxErrors)
//     {
//         case 1: find<0, 1>(ossContext, delegate, delegateDirect, index, needles, TDistanceTag());
//                 break;
//         case 2: find<0, 2>(ossContext, delegate, delegateDirect, index, needles, TDistanceTag());
//                 break;
//         case 3: find<0, 3>(ossContext, delegate, delegateDirect, index, needles, TDistanceTag());
//                 break;
//         case 4: find<0, 4>(ossContext, delegate, delegateDirect, index, needles, TDistanceTag());
//                 break;
//         default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
//                 exit(1);
//     }
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
    uint16_t maxError = me.options.errorRate * length(readSeqs[0]); //maybe include read length as input parameter
    uint16_t strata = disOptions.strataRate * length(readSeqs[0]);
//     ossContext.maxError = maxError;
//     ossContext.strata = strata;

    bool mscheme = false;
    ossContext.setReadContextOSS(maxError, strata, mscheme);

    std::cout << "Using 0 and " << maxError << " Scheme" << "\n";

//     myTfind(0, maxError, ossContext, delegate, delegateDirect, me.biIndex, readSeqs, EditDistance());
//     myfind(0, maxError, ossContext, delegate, delegateDirect, me.biIndex, readSeqs, HammingDistance());

    find(0, maxError, ossContext, delegate, delegateDirect, me.biIndex, readSeqs, EditDistance());

//     if(addMyOwnOption) //IsSameType<typename TConfig::TSeedsDistance, EditDistance>::VALUE
//         find(0, maxError, ossContext, delegate, delegateDirect, me.biIndex, readSeqs, EditDistance());
//     else
//         find(0, maxError, ossContext, delegate, delegateDirect, me.biIndex, readSeqs, HammingDistance());

/*
    std::cout << "Finished Searching: " << "\n";
    typedef typename TTraits::TMatch                             TMatch;
    for(int i = 0; i < length(me.matchesByCoord); ++i){
        std::cout << "Match: " << "\n";
        TMatch myMatch = me.matchesByCoord[i];
        std::cout << myMatch.contigBegin << "\n";
        std::cout << myMatch.contigEnd << "\n";
        std::cout << "ContigId: " << myMatch.contigId << "\n";
        std::cout << "Errors: " << myMatch.errors << "\n";
        std::cout << "ReadId: "<< myMatch.readId << "\n" << "\n";
    }
    std::cout << "End of Matches ____________________________________________________________________________________________________ " << "\n";*/

    aggregateMatchesOSS(me, readSeqs);
//     aggregateMatches(me, readSeqs);

    rankMatches(me, me.reads.seqs);
    if (me.options.verifyMatches)
        verifyMatches(me);
    alignMatches(me);
    copyMatches(mainMapper, me, disOptions);
    copyCigars(mainMapper, me, disOptions);
    appendStats(mainMapper, me);


//     Iter<TBiIndex, VSTree<TopDown<> > > iter(biindex);
//     Iter<TIndex, VSTree<TopDown<> > > iterUni(index);
//
//     Iter<TBiIndex, VSTree<TopDown<> > > iter2(biindex);
//     Iter<TIndex, VSTree<TopDown<> > > iterUni2(index);
//
// //                    1234567890123456
//     DnaString test =    "GGGTCGCGGTGCGCGGCGACGAAGG";
//     DnaString testrev = "GGAAGCAGCGGCGCGTGGCGCTGGG";

//     int k = 0;
//     while(k < length(testrev)){
//         std::cout << "k: " << k << "\n";
//         std::cout << iter.fwdIter.vDesc.range.i1 << ":" << iter.fwdIter.vDesc.range.i2 << "\n";
//         if (!goDown(iter, testrev[k], Fwd())){
//             std::cout << "Stop" << "\n";
//             std::cout << iter.fwdIter.vDesc.range.i1 << ":" << iter.fwdIter.vDesc.range.i2 << "\n";
//             break;
//         }
//         ++k;
//     }
//

//     std::cout << "Unidirectional Index" << "\n";
//     std::cout << iterUni.vDesc.range.i1 << ":" << iterUni.vDesc.range.i2 << "\n";
//     std::cout << goDown(iterUni, testrev) << "\n";
//     std::cout << iterUni.vDesc.range.i1 << ":" << iterUni.vDesc.range.i2 << "\n";
//
//     std::cout << "BidirectionalIndex" << "\n";
//     std::cout << iter.fwdIter.vDesc.range.i1 << ":" << iter.fwdIter.vDesc.range.i2 << "\n";
//     std::cout << goDown(iter, testrev, Rev()) << "\n";
//     std::cout << iter.fwdIter.vDesc.range.i1 << ":" << iter.fwdIter.vDesc.range.i2 << "\n";
//
//     std::cout << "Unidirectional Index" << "\n";
//     std::cout << iterUni2.vDesc.range.i1 << ":" << iterUni2.vDesc.range.i2 << "\n";
//     std::cout << goDown(iterUni2, test) << "\n";
//     std::cout << iterUni2.vDesc.range.i1 << ":" << iterUni2.vDesc.range.i2 << "\n";
//
//     std::cout << "BidirectionalIndex" << "\n";
//     std::cout << iter2.fwdIter.vDesc.range.i1 << ":" << iter2.fwdIter.vDesc.range.i2 << "\n";
//     std::cout << iter2.revIter.vDesc.range.i1 << ":" << iter2.revIter.vDesc.range.i2 << "\n";
//     std::cout << goDown(iter2, test[0], Fwd()) << "\n";
//     std::cout << goDown(iter2, test[0], Rev()) << "\n";
//     std::cout << iter2.fwdIter.vDesc.range.i1 << ":" << iter2.fwdIter.vDesc.range.i2 << "\n";
//     std::cout << iter2.revIter.vDesc.range.i1 << ":" << iter2.revIter.vDesc.range.i2 << "\n";


//     isMapped(me.ctx, readId)
//     getMinErrors(me.ctx, readId) + strata <= error i am currently searching for

/*
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
    aggregateMatches(me, readSeqs);
    rankMatches(me, me.reads.seqs);
    if (me.options.verifyMatches)
        verifyMatches(me);
    alignMatches(me);
    copyMatches(mainMapper, me, disOptions);
    copyCigars(mainMapper, me, disOptions);
    appendStats(mainMapper, me);*/

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
    aggregateMatches(mainMapper, mainMapper.reads.seqs);
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
            for(uint32_t i = 0; i < disOptions.numberOfBins; ++i){
                std::cout << "In bin Number: " << i << "\n";
                disOptions.currentBinNo = i;
                Options options = mainMapper.options;
                appendFileName(options.contigsIndexFile, disOptions.IndicesDirectory, i);
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
