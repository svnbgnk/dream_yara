#ifndef APP_YARA_MAPPABILITY_H_
#define APP_YARA_MAPPABILITY_H


#include <vector>
#include <cstdint>
#include <limits>


using namespace std;
using namespace seqan;

#include "common.h"
#include "algo1.hpp"
#include "algo2.hpp"
#include "algo3.hpp"
#include "algo4.hpp"


// ----------------------------------------------------------------------------
// Class OptionsM
// ----------------------------------------------------------------------------

struct OptionsM
{

    CharString      IndicesDirectory;
    CharString      output;
    unsigned        numberOfBins;

    unsigned        errors;
    unsigned        strata;
    unsigned        k_length;
    unsigned        threshold;
    unsigned        overlap = 0;

    bool            indels;
    bool            high;

    unsigned        currentBinNo;
    bool            trivial;
    unsigned        threads;
    bool            verbose;

    //Parameters for single index
    CharString          contigsIndexFile;
    uint64_t            contigsSize;
    uint64_t            contigsMaxLength;
    uint64_t            contigsSum;
    std::vector<uint32_t>   contigOffsets;

    OptionsM() :
    numberOfBins(64),
    currentBinNo(0),
    threads(1),
    verbose(false),
    trivial(false)
    {}
};


string get_output_path(OptionsM const & opt)
{
    string output_path = toCString(opt.output);
    if(opt.indels)
        output_path += "/mappability_" + to_string(opt.errors) + "_" + to_string(opt.k_length - opt.errors)/* + "_" + to_string(opt.overlap)*/;
    else
        output_path += "/mappability_" + to_string(opt.errors) + "_" + to_string(opt.k_length)/* + "_" + to_string(opt.overlap)*/;
    output_path += ".gmapp" + string(opt.high ? "16" : "8");
    return output_path;
}

template <typename T>
inline void save(vector<T> const & c, string const & output_path, OptionsM const & opt)
{
    //edit many 0 to end
    if(opt.verbose && c.size() < 50000)
    {
        for(int i = 0; i < c.size(); ++i)
            std::cout << (int)c[i] << " ";
        std::cout << "\n";
    }

    ofstream outfile(output_path, ios::out | ios::binary);
    outfile.write((const char*) &c[0], c.size() * sizeof(T));
    outfile.close();
}

template <typename value_type, typename TIndex, typename TText, typename TSeqLengths, typename TDistanceTag>
inline void run(TIndex & index, TText const & text, TSeqLengths const & sL, OptionsM const & opt, TDistanceTag const &)
{
    // add Zeroes at the end to make up for smaller mappability for using a larger k_length
    uint8_t cmod = (opt.indels) ? (opt.errors) : 0;

    vector<value_type> c(length(text) - opt.k_length + cmod + 1, 0);
    if(opt.trivial){
        switch (opt.errors)
        {
            case 0:  runAlgoTrivial<0>(index, text, sL, c, opt);
                    break;
            case 1:  runAlgoTrivial<1>(index, text, sL, c, opt);
                    break;
            case 2:  runAlgoTrivial<2>(index, text, sL, c, opt);
                    break;
            case 3:  runAlgoTrivial<3>(index, text, sL, c, opt);
                    break;
            case 4:  runAlgoTrivial<4>(index, text, sL, c, opt);
                    break;
            default: cerr << "E = " << opt.errors << " not yet supported.\n";
                    exit(1);
        }
    }
    else
    {
        switch (opt.errors)
        {
            case 0:  runAlgo4<0>(index, text, sL, c, opt, TDistanceTag());
                    break;
            case 1:  runAlgo4<1>(index, text, sL, c, opt, TDistanceTag());
                    break;
            case 2:  runAlgo4<2>(index, text, sL, c, opt, TDistanceTag());
                    break;
            case 3:  runAlgo4<3>(index, text, sL, c, opt, TDistanceTag());
                    break;
            case 4:  runAlgo4<4>(index, text, sL, c, opt, TDistanceTag());
                    break;
            default: cerr << "E = " << opt.errors << " not yet supported.\n";
                    exit(1);
        }
    }

//     if (SearchParams::outputProgress)
//         std::cout << '\r';
//     std::cout << "Progress: 100.00%\n" << std::flush;
//     cout.flush();
    string output_path = get_output_path(opt);
    save(c, output_path, opt);
}


template <typename value_type, typename TIndex, typename TText>
inline void run(TIndex & index, TText const & text, OptionsM const & opt)
{
    std::vector<uint64_t> sequenceLengths = getSeqLengths<uint64_t, uint64_t>(text);
//     std::cout << "Number of Sequences: " << seqan::length(text) << "\n";


    if (opt.indels) {
        run<value_type>(index, concat(text), sequenceLengths, opt, EditDistance());
    }
    else
        run<value_type>(index, concat(text), sequenceLengths, opt, HammingDistance());
}


template <typename TIndex, typename TText>
inline void calcMappa(TIndex & index, TText const & text, OptionsM const & opt)
{
    if (opt.high) {
        run<uint16_t>(index, text, opt);
    }
    else
        run<uint8_t>(index, text, opt);
}

#endif
