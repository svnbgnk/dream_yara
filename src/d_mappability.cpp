// ==========================================================================
//                                 d_indexer.cpp
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

// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------
#include <string>
#include <vector>
#include <sys/stat.h>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <mutex>
#include <condition_variable>
#include <future>
#include <thread>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------
#include "store_seqs.h"
#include "misc_timer.h"
#include "misc_tags.h"
#include "misc_types.h"
#include "bits_matches.h"
#include "misc_options.h"
#include "d_misc_options.h"
#include "index_fm.h"

#include "mappability.h"
#include "createBit.h"

using namespace seqan;


#define YARA_MAPPABILITY



// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser, OptionsM const & options)
{
    setAppName(parser, "dream_yara_mappability");
    setShortDescription(parser, "DREAM-Yara Mappability and Bitvector calculation");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIREFERENCE INDEX DIRECTORY\\fP> <\\fISE-READS FILE\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "REFERENCE INDEX DIRECTORY"));
    setHelpText(parser, 0, "A directory containing multiple indices of reference genomes.");

    addOption(parser, ArgParseOption("b", "number-of-bins", "The number of bins (indices) for distributed mapper",
                                     ArgParseOption::INTEGER));
    setRequired(parser, "number-of-bins");
    setMinValue(parser, "number-of-bins", "1");
    setMaxValue(parser, "number-of-bins", "1024");

    addSection(parser, "Output OptionsM");

    addOption(parser, ArgParseOption("O", "output", "Path to output directory (error number, length and overlap will be appended to the output file)", ArgParseArgument::OUTPUT_FILE, "OUT"));
//     setRequired(parser, "output");

//    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "REFERENCE FILE DIR"));
//    //    setValidValues(parser, 0, SeqFileIn::getFileExtensions());
//    setHelpText(parser, 0, "A directory containing reference genome files.");

    addOption(parser, ArgParseOption("E", "errors", "Number of errors", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("K", "length", "Length of k-mers", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "length");

    addOption(parser, ArgParseOption("i", "indels", "Turns on indels (EditDistance). "
        "If not selected, only mismatches will be considered."));


    addOption(parser, ArgParseOption("triv", "trivial", "Use trivial Mappability Calculation"));

    addOption(parser, ArgParseOption("hi", "high", "Stores the mappability vector in 16 bit unsigned integers instead of 8 bit (max. value 65535 instead of 255)"));

    addOption(parser, ArgParseOption("o", "overlap", "Number of overlapping reads (o + 1 Strings will be searched at once beginning with their overlap region)", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("t", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", omp_get_max_threads());

    addSection(parser, "Additional Parameter for bitvectors");

    addOption(parser, ArgParseOption("T", "threshold", "Number of times a k-mer can occure and still be accepted as mappable", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("s", "strata", "Max errors allowed during mapping", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", "Verbose"));

}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(OptionsM & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;


    // Parse indexed genome input file.
    getArgumentValue(options.IndicesDirectory, parser, 0);

    // Append trailing slash if it doesn't exist.
    appendTrailingSlash(options.IndicesDirectory);

    if(isSet(parser, "output")){
        getOptionValue(options.output, parser, "output");
    }else{
        CharString path = options.IndicesDirectory;
        path += "mappability";
        options.output = path;
    }

    getOptionValue(options.numberOfBins, parser, "number-of-bins");

    getOptionValue(options.errors, parser, "errors");

    getOptionValue(options.k_length, parser, "length");

    getOptionValue(options.overlap, parser, "overlap");

    getOptionValue(options.threads, parser, "threads");

    getOptionValue(options.threshold, parser, "threshold");

    if(isSet(parser, "strata"))
        getOptionValue(options.strata, parser, "strata");
    else
        getOptionValue(options.strata, parser, "errors");


    options.indels = isSet(parser, "indels");
    options.trivial = isSet(parser, "trivial");
    options.high = isSet(parser, "high");
    options.verbose = isSet(parser, "verbose");

    if(!isSet(parser, "trivial")){
        if(!isSet(parser, "overlap")){
            std::cerr << "Overlap needs to be Set when using the non-trivial algorithm.\n";
            exit(1);
        }
    }


    return ArgumentParser::PARSE_OK;
}



template <typename TContigsSize,
          typename TContigsLen,
          typename TContigsSum>
inline void runMappability(OptionsM const & options,
                      OptionsM & disOptions)
{

    //load Text
    typedef SeqStore<void, YaraContigsConfig<Alloc<> > >              TContigs; //MMap<>
    TContigs text;
    try
    {
        if (!open(text, toCString(options.contigsIndexFile), OPEN_RDONLY))
            throw RuntimeError("Error while opening reference file.");
    }
    catch (BadAlloc const & ) // e
    {
        throw RuntimeError("Insufficient memory to load the reference.");
    }



    //load BidirectionalIndex
    typedef YaraFMConfig<TContigsSize, TContigsLen, TContigsSum, Alloc<> > TIndexConfig; //MMap<>
    typedef FMIndex<void, TIndexConfig>                             TIndexSpec;
    typedef BidirectionalIndex<TIndexSpec>                          TBiIndexSpec;
    typedef Index<typename TIndexConfig::Text, TIndexSpec>          TIndex;
    typedef Index<typename TIndexConfig::Text, TBiIndexSpec>        TBiIndex;

    TBiIndex biIndex;
    try
    {
        if (!open(biIndex, toCString(options.contigsIndexFile), OPEN_RDONLY))
            throw RuntimeError("Error while opening reference index file.");
    }
    catch (BadAlloc const &  e ) // e
    {
        throw RuntimeError("Insufficient memory to load the reference index.");
    }

    std::string mappability_path = toCString(options.output);
    mappability_path += "/mappability_" + to_string(options.errors) + "_" + to_string(options.k_length);
    mappability_path += ".gmapp" + string(options.high ? "16" : "8");

    //calculate mappability
    if(!file_exists(mappability_path))
    {
        calcMappa(biIndex, text.seqs, options);
    }else{
        std::cout << "Skip mappability Calculation" << "\n";
    }

    if(!file_exists(mappability_path))
    {
        std::cerr << "Cannot find mappability file" << "\n";
        exit(0);
    }

    std::string bitvectors_dir = toCString(options.output);
    bitvectors_dir += "/";

    //load mappability file
    if(options.verbose)
        std::cout << "Errors:" << options.errors << "\tStrata: " << options.strata << "\n"; //TODO comment
    vector<uint8_t> mappability = read(mappability_path);

    if(options.verbose)
        std::cout << "Loaded Mappability vector. Size: " << mappability.size() << "\n";
    bitvectors result = create_all_bit_vectors<TContigsSum>(mappability, options.k_length, options.threshold, options.errors, options.strata, options.threads, options.verbose);
    if(options.verbose)
        std::cout << "Finished bit vectors." << "\n";

    order_bit_vector<TContigsSize, TContigsLen, TContigsSum>(biIndex, text.seqs, result, options.threads, options.verbose);

    //save Bitvectors
    if(options.verbose)
        std::cout << "Finished sorting" << "\n";
    for(uint32_t i = 0; i < result.bv.size(); ++i){
        sdsl::store_to_file(result.bv[i], toCString(bitvectors_dir)  + result.names[i]);
    }

}

template <typename TContigsSize,
          typename TContigsLen>
void configureMappability(OptionsM const & options,
                     OptionsM & disOptions)
{
    if (options.contigsSum <= MaxValue<uint32_t>::VALUE)
    {
        runMappability<TContigsSize, TContigsLen, uint32_t>(options, disOptions);
    }
    else
    {
        runMappability<TContigsSize, TContigsLen, uint64_t>(options, disOptions);
    }
}



template <typename TContigsSize>
void configureMappability(OptionsM const & options,
                     OptionsM & disOptions)
{
    if (options.contigsMaxLength <= MaxValue<uint32_t>::VALUE)
    {
        configureMappability<TContigsSize, uint32_t>(options, disOptions);
    }
    else
    {
#ifdef DR_YARA_LARGE_CONTIGS
        configureMappability<TContigsSize, uint64_t>(options, disOptions);
#else
        throw RuntimeError("Maximum contig length exceeded. Recompile with -DDR_YARA_LARGE_CONTIGS=ON.");
#endif
    }
}

void configureMappability(OptionsM const & options,
                     OptionsM & disOptions)
{
    if (options.contigsSize <= MaxValue<uint8_t>::VALUE)
    {
        configureMappability<uint8_t>(options, disOptions);
    }
    else if (options.contigsSize <= MaxValue<uint16_t>::VALUE)
    {
        configureMappability<uint16_t>(options, disOptions);
    }
    else
    {
#ifdef DR_YARA_LARGE_CONTIGS
        configureMappability<uint32_t>(options, disOptions);
#else
        throw RuntimeError("Maximum number of contigs exceeded. Recompile with -DYARA_LARGE_CONTIGS=ON.");
#endif
    }
}



// ----------------------------------------------------------------------------
// Function runDisMapper()
// ----------------------------------------------------------------------------
inline void runDisMappability(OptionsM & options)
{
    for(uint32_t i = 0; i < options.numberOfBins; ++i){
        std::cout << "In bin Number: " << i << "\n";

        //create_dir
        CharString path = options.output;
        path += "/";
        path += to_string(i);
        std::string path_dir = toCString(path);

        const int dir_test = mkdir(path_dir.data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_test)
            std::cout << "Error creating directory!\n";

        options.currentBinNo = i;
        OptionsM binOption = options;
        binOption.output = path;
        appendFileName(binOption.contigsIndexFile, options.IndicesDirectory, i);
        if (!openContigsLimits(binOption))
            throw RuntimeError("Error while opening reference file.");
        configureMappability(binOption, options);
    }
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    OptionsM options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;


        if (options.overlap > options.k_length - 1)
    {
        cerr << "ERROR: overlap cannot be larger than K - 1.\n";
        exit(1);
    }

    if (!(options.k_length - options.overlap >= options.errors + 2))
    {
        cerr << "ERROR: overlap should be at least K - E - 2. (K - O >= E + 2 must hold since common overlap has length K - O and will be split into E + 2 parts).\n";
        exit(1);
    }

    if(options.strata > options.errors){
        cerr << "ERROR: strata should be smaller than maximum allowed Error.\n";
        exit(1);
    }

    options.overlap = options.k_length - options.overlap;


    //create_dir
    std::string path_dir = toCString(options.output);

    const int dir_test = mkdir(path_dir.data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (-1 == dir_test)
        std::cerr << "Dir already exist\n";

//    std::string comExt = commonExtension(options.contigsDir, options.numberOfBins);
//     typedef std::map<uint32_t,CharString>::iterator mapIter;


    try
    {
//         Timer<double>       timer;
        Timer<double>       globalTimer;
        start (globalTimer);

        runDisMappability(options);

        stop(globalTimer);
        std::cerr <<"\nFinshed in \t\t\t" << globalTimer << std::endl;
    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
