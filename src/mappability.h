#ifndef APP_YARA_MAPPABILITY_H_
#define APP_YARA_MAPPABILITY_H


#include <vector>
#include <cstdint>
#include <limits>
#include <sys/stat.h>

// #include <seqan/arg_parse.h>
// #include <seqan/seq_io.h>
// #include <seqan/index.h>

using namespace std;
using namespace seqan;

#include "common.h"
#include "algo1.hpp"

// #include "algo2.hpp"
// #include "algo3.hpp"
// #include "algo4.hpp"

// #include "common_auxiliary.h"
// #include "find2_index_approx_extension.h"


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
    unsigned        overlap;

    bool            indels;
    bool            high;
    bool            mmap;


    unsigned        currentBinNo;
    unsigned        threads;
    bool            verbose;

    //Parameters for single index
    CharString      contigsIndexFile;
    uint64_t            contigsSize;
    uint64_t            contigsMaxLength;
    uint64_t            contigsSum;
    std::vector<uint32_t>   contigOffsets;

    OptionsM() :
    numberOfBins(64),
    currentBinNo(0),
    threads(1),
    verbose(false)
    {}
};


string get_output_path(OptionsM const & opt)
{
    string output_path = toCString(opt.output);
    output_path += "/mappability_" + to_string(opt.errors) + "_" + to_string(opt.k_length)/* + "_" + to_string(opt.overlap)*/;
    output_path += ".gmapp" + string(opt.high ? "16" : "8");
    return output_path;
}

template <typename T>
inline void save(vector<T> const & c, string const & output_path, OptionsM const & opt)
{

/*
    for(int i = 0; i < c.size(); ++i)
        std::cout << (int)c[i] << " ";
    std::cout << "\n";*/

    ofstream outfile(output_path, ios::out | ios::binary);
    outfile.write((const char*) &c[0], c.size() * sizeof(T));
    outfile.close();
}

template <typename TDistance, typename value_type, typename TIndex, typename TText>
inline void run(TIndex & index, TText const & text, OptionsM const & opt)
{
//     Dna5String const & infix_n = infix(text, 0, 100);

    vector<value_type> c(length(text) - opt.k_length + 1, 0);
//     if(opt.indels){
        switch (opt.errors)
        {
            case 0:  runAlgoTrivial<0>(index, text, c, opt);
                    break;
            case 1:  runAlgoTrivial<1>(index, text, c, opt);
                    break;
            case 2:  runAlgoTrivial<2>(index, text, c, opt);
                    break;
            case 3:  runAlgoTrivial<3>(index, text, c, opt);
                    break;
            case 4:  runAlgoTrivial<4>(index, text, c, opt);
                    break;
            default: cerr << "E = " << opt.errors << " not yet supported.\n";
                    exit(1);
        }

        /*
    }
    else
    {
        switch (opt.errors)
        {
            case 0:  runAlgo4<0>(index, text, c, opt);
                    break;
            case 1:  runAlgo4<1>(index, text, c, opt);
                    break;
            case 2:  runAlgo4<2>(index, text, c, opt);
                    break;
            case 3:  runAlgo4<3>(index, text, c, opt);
                    break;
            case 4:  runAlgo4<4>(index, text, c, opt);
                    break;
            default: cerr << "E = " << opt.errors << " not yet supported.\n";
                    exit(1);
        }
    }*/

//     if (SearchParams::outputProgress)
//         std::cout << '\r';
//     std::cout << "Progress: 100.00%\n" << std::flush;
//     cout.flush();
    string output_path = get_output_path(opt);
    save(c, output_path, opt);
}

template <typename TDistance, typename TIndex, typename TText>
inline void run(TIndex & index, TText const & text, OptionsM const & opt)
{
    if (opt.high) {
        run<TDistance, uint16_t>(index, concat(text), opt);
    }
    else
        run<TDistance, uint8_t>(index, concat(text), opt);
}

template<typename TIndex, typename TText>
inline void calcMappa(TIndex & index, TText const & text, OptionsM const & opt)
{
    if (opt.indels) {
        run<EditDistance>(index, text, opt);
    }
    else
        run<HammingDistance>(index, text, opt);
}



#endif