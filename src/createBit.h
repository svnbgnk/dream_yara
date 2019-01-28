#ifndef APP_YARA_CREATEBIT_H
#define APP_YARA_CREATEBIT_H


#include <vector>
#include <cstdint>
#include <limits>
#include "common.h"

// #include <seqan/arg_parse.h>
// #include <seqan/seq_io.h>
// #include <seqan/index.h>

using namespace std;
using namespace seqan;


struct bitvectors
{
    vector<bool> fwdd;
    vector<string> names;
    vector<sdsl::bit_vector> bv;
};

std::vector<int> getInt(std::string const& mappability_str)
{
  std::istringstream iss(mappability_str);
  return std::vector<int>{
    std::istream_iterator<int>(iss),
    std::istream_iterator<int>()
  };
}

vector<uint8_t> read(const string mappability_path){
    string mappability_str;

    vector<uint8_t> mappability_int;
    ifstream file(toCString(mappability_path), std::ios::binary);
    if (!file.eof() && !file.fail())
    {
        file.seekg(0, std::ios_base::end);
        std::streampos fileSize = file.tellg();
        mappability_int.resize(fileSize);
        file.seekg(0, std::ios_base::beg);
        file.read((char*)&mappability_int[0], fileSize);
        file.close();
        }
    return(mappability_int);
}

template<typename TVector, typename TElem>
bool checkForElem(TVector const & v, TElem const & e)
{
   for(uint16_t i = 0; i < v.size(); ++i){
       if(v[i].compare(e) == 0)
           return true;
   }
   return false;
}

template <typename TContigsSum>
bitvectors create_all_bit_vectors(const vector <uint8_t> & mappability,
                                  uint32_t const len, uint32_t const threshold, uint8_t const errors, uint8_t const strata, bool const indels,  uint32_t const mythreads, bool const verbose){

    //TODO switch left and right in the moment they discribe in which direction the k-mere is
    bitvectors b;
    sdsl::bit_vector righti (mappability.size() + len - 1, 0);
    sdsl::bit_vector lefti (mappability.size() + len - 1, 0);
    uint32_t anchored_shift = (indels) ? (len + errors) : (len);
    #pragma omp parallel for schedule(static) //num_threads(mythreads)
    for(TContigsSum i = 0; i < mappability.size(); ++i){
        lefti[i + anchored_shift - 1] = (mappability[i] <= threshold);
        righti[i] = (mappability[i] <= threshold);
    }
    if(verbose)
        cout << "Finished Default Bit Vectors.  Length: " << righti.size() << endl;


    b.bv.push_back(righti);
    b.names.push_back("left_anchored_bvector_" + to_string(len) + "_shift_0");
    b.fwdd.push_back(true);
    b.bv.push_back(lefti);
    b.names.push_back("right_anchored_bvector_" + to_string(len) + "_shift_0");
    b.fwdd.push_back(false);

    if(verbose)
        std::cout << "\nAdditonal bitvectors besides left_anchored and right_anchored bitvector: \n";

    vector<uint32_t> shift_r;
    vector<uint32_t> shift_l;

    for(uint8_t s = 0; s <= errors - strata; ++s){
        uint8_t se = errors - s;
        if(verbose)
            std::cout << "Bitvectors for Scheme <" << 0 << ", " << (int)se << ">: \n";
        auto blocklengths = loadBlockLimits(se, len);
        shift_r = blocklengths.first;
        shift_l = blocklengths.second;

        for(uint16_t i = 0; i < shift_r.size(); ++i){
            uint32_t shift = shift_r[i];
            bool skip = checkForElem(b.names, ("left_anchored_bvector_" + to_string(len) + "_shift_" + to_string(shift)));
            if(skip){
                if(verbose)
                    std::cout << "Bitvector already included" << "\n";
                continue;
            }

            b.names.push_back("left_anchored_bvector_" + to_string(len) + "_shift_" + to_string(shift));
            if(verbose)
                std::cout << b.names.back()<< "\n";

            sdsl::bit_vector newright(mappability.size() + len - 1, 0);
            #pragma omp parallel for schedule(static) //num_threads(mythreads)
            for(TContigsSum j = 0; j < righti.size(); ++j){
                if(j >= shift)
                    newright[j] = righti[j - shift];
            }
            b.bv.push_back(newright);
            b.fwdd.push_back(true);
        }

        for(uint16_t i = 0; i < shift_l.size(); ++i){
            uint32_t shift = shift_l[i];
            bool skip = checkForElem(b.names, ("right_anchored_bvector_" + to_string(len) + "_shift_" + to_string(shift)));
            if(skip){
                if(verbose)
                    std::cout << "Bitvector already included" << "\n";
                continue;
            }

            b.names.push_back("right_anchored_bvector_" + to_string(len) + "_shift_" + to_string(shift));
            if(verbose)
                std::cout << b.names.back() << "\n";

            sdsl::bit_vector newleft(mappability.size() + len - 1, 0);
            #pragma omp parallel for schedule(static) //num_threads(mythreads)
            for(TContigsSum j = 0; j < righti.size(); ++j){
                if(j + shift < lefti.size() - 1)
                    newleft[j] = lefti[j + shift];
            }
            b.bv.push_back(newleft);
            b.fwdd.push_back(false);
        }
    }

    if(verbose)
        std::cout << "Number of Bitvectors: " << b.names.size() << "\n\n";
    return(b);
}

template <typename TContigsSize, typename TContigsLen, typename TContigsSum, typename TIndex, typename TText>
void order_bit_vector(TIndex & index, TText const & text, bitvectors & b, uint32_t const mythreads, bool const verbose)
{
    TContigsSum indexSize = seqan::length(index.fwd.sa);
    if(verbose)
        cout << "Loaded Index. Size:" << indexSize << endl;
    vector<sdsl::bit_vector> bit_vectors_ordered (b.bv);

//     uint32_t number_of_indeces = countSequences(index);
    TContigsSize number_of_indeces = length(text);

    std::vector<TContigsLen> sequenceLengths = getSeqLengths<TContigsSize, TContigsLen>(text);
    if(verbose){

        std::cout << "cumSequenceLengths: \n";
        /*
        auto cumSequenceLengths = stringSetLimits(text);
        for(int i = 0; i < length(cumSequenceLengths); ++i){
            std::cout << cumSequenceLengths << "\t";

        }*/


        for(int i = 0; i < sequenceLengths.size(); ++i)
            std::cout << sequenceLengths[i] << "\t";
        std::cout << "\n";
    }

//     std::vector<uint32_t> sequenceLengths = getSeqLengths(index);
//     for(uint32_t i = 2; i < sequenceLengths.size(); ++i)
//         sequenceLengths[i] += sequenceLengths[i - 1];


    if(verbose){
        std::cout << "\nNumber of Sequences: " << (int)number_of_indeces << "\n";
        cout << "Start sorting bitvectors with " << mythreads << "threads" << endl;
    }

    //dynamic since tasks can take different amount of time (SA Sampling) critical to guarantee thread safety
    uint16_t bsize = b.bv.size();

//     #pragma omp parallel for schedule(dynamic) num_threads(mythreads)
    #pragma omp parallel for schedule(static, (indexSize/(mythreads*100))) num_threads(mythreads)
    for (TContigsSum j = 0; j < indexSize - number_of_indeces; ++j)
    {
        // skip sentinels
        Pair<TContigsSize, TContigsLen> sa_f = index.fwd.sa[j + number_of_indeces];
        Pair<TContigsSize, TContigsLen> sa_r = index.rev.sa[j + number_of_indeces];

//         Pair<uint16_t, uint32_t> sa_f = index.fwd.sa[j + number_of_indeces];
//         Pair<uint16_t, uint32_t> sa_r = index.rev.sa[j + number_of_indeces];
        TContigsSum fpos = getSeqOffset(sa_f) + sequenceLengths[getSeqNo(sa_f)];
        TContigsSum rpos = sequenceLengths[getSeqNo(sa_r) + 1] - getSeqOffset(sa_r) - 1;
        vector<bool> values(bsize);

        for(uint16_t i = 0; i < bsize; ++i){
            if(b.fwdd[i]){
                values[i] = b.bv[i][fpos];
            }
            else
            {
                values[i] = b.bv[i][rpos];
            }
        }
        #pragma omp critical
        {
        for(uint16_t i = 0; i < bsize; ++i)
            bit_vectors_ordered[i][j] = values[i];
        }
    }
    b.bv = bit_vectors_ordered;
}



#endif
