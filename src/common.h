#ifndef COMMON_H_
#define COMMON_H_

#include <sdsl/bit_vectors.hpp>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>


// reduce space consumption
// requires genome not to have more than ~ 4 gigabases
// multi-sequence fasta file must contain less than ~64k sequences of at most ~ 4 gigabases in total
namespace seqan {
/*
template <typename TChar, typename TAlloc, typename TOwner>
struct SAValue<StringSet<String<TChar, TAlloc>, TOwner> >
{
    typedef Pair<uint16_t, uint32_t> Type;
};

template <typename TChar, typename TAlloc>
struct SAValue<String<TChar, TAlloc> >
{
    typedef uint32_t Type;
};

template <typename TSpec = void, typename TLengthSum = size_t, unsigned LEVELS = 2, unsigned WORDS_PER_BLOCK = 1>
struct FastFMIndexConfigS1
{
     typedef TLengthSum                                                                         LengthSum;
     typedef Levels<void, LevelsPrefixRDConfig<LengthSum, Alloc<>, LEVELS, WORDS_PER_BLOCK> >   Bwt;
     typedef Levels<void, LevelsRDConfig<LengthSum, Alloc<>, LEVELS, WORDS_PER_BLOCK> >         Sentinels;
     static const unsigned SAMPLING = 1;
};

using TMyFastConfig = seqan::FastFMIndexConfig<void, uint32_t, 2, 1>;
using TIndexConfig = seqan::BidirectionalIndex<seqan::FMIndex<void, TMyFastConfig> >;
//TODO ask about this
// using TMyFastConfig = seqan::FastFMIndexConfig<void, uint32_t, 2, 1>;
// using TIndexConfig = seqan::BidirectionalIndex<seqan::FMIndex<RadixSortSACreateTag, TMyFastConfig> >;

// template <typename TText>
// using TIndex = seqan::Index<TText, TIndexConfig>;


// using TMyFastConfig = seqan::FastFMIndexConfigS1<void, uint32_t, 2, 1>;
// using TIndexConfig = seqan::BidirectionalIndex<seqan::FMIndex<void, TMyFastConfig> >;
//
//
// typedef FastFMIndexConfigS1<void, uint32_t, 2, 1> TMyFastConfig;
// typedef BidirectionalIndex<FMIndex<void, TMyFastConfig> > TIndexConfig;

typedef String<Dna, Alloc<>> TString;
typedef StringSet<TString, Owner<ConcatDirect<> > > TText;
typedef Index<TText, TIndexConfig> MyIndex;*/


typedef sdsl::bit_vector TBitvector;
typedef sdsl::rank_support_v<> TSupport;
// typedef sdsl::rank_support_v5<> TSupport;
/*
struct hit{
    bool rev;
    Pair <uint32_t, uint32_t> occ;
    Pair <uint32_t, uint32_t> occEnd;
    uint8_t errors;
    uint32_t readId;
    DnaString read;
};
*/

template <size_t nbrBlocks, size_t N>
inline void calcConstParameters(std::array<OptimalSearch<nbrBlocks>, N> & ss)
{
    uint8_t id = 0;
    for (OptimalSearch<nbrBlocks> & s : ss){
        s.id = id;
        ++id;
//         int bsize = s.pi.size();
        uint8_t min = s.pi[0];
        uint8_t max = s.pi[0];
        // maybe < N?
        for(int i = 0; i < nbrBlocks; ++i){
            if(min > s.pi[i])
                min = s.pi[i];
            if(max < s.pi[i])
                max = s.pi[i];
            s.min[i] = min;
            s.max[i] = max;
        }
        uint8_t lastValue = s.pi[nbrBlocks - 1];
        int k = (nbrBlocks > 1) ? nbrBlocks - 2 : 0;
        while(k >= 0){
            if(s.pi[k] == lastValue - 1 || s.pi[k] == lastValue + 1)
            {
                lastValue = s.pi[k];
                --k;
            }else{
                s.startUniDir = k + 1;
                break;
            }
        }
    }
}



template <size_t nbrBlocks, size_t N>
/*constexpr */inline void _optimalSearchSchemeComputeChronBlocklength(std::array<OptimalSearch<nbrBlocks>, N> & ss)
{

    int16_t blockSize = nbrBlocks;
    for (OptimalSearch<nbrBlocks> & s : ss){
        s.chronBL[s.pi[0] - 1]  = s.blocklength[0];
        for(int16_t j = 1; j < blockSize; ++j){
            s.chronBL[s.pi[j] - 1] = s.blocklength[j] -  s.blocklength[j - 1];
        }

        for(int16_t j = 1; j < blockSize; ++j){
            std::cerr << ""; //if blockSize == 2 compiler skips the loop??
            s.chronBL[j] += s.chronBL[j - 1];
        }

//         int32_t tmp = s.chronBL[1]; // if blockSize == 2 compiler skip previous loop
//         std::cout << tmp << "\n";


        if(blockSize > 1){
            s.revChronBL[s.pi[blockSize - 1] - 1]  = s.blocklength[blockSize - 1] - s.blocklength[blockSize - 2];
            for(int16_t i = blockSize - 2; i >= 0; --i){
                s.revChronBL[s.pi[i] - 1] = s.blocklength[i] - ((i > 0) ? s.blocklength[i - 1] : 0);
            }

            for(int16_t i = blockSize - 2; i >= 0; --i)
                s.revChronBL[i] += s.revChronBL[i + 1];
        }
        else
        {
            s.revChronBL[0] = s.blocklength[0];
        }
    }

    for (OptimalSearch<nbrBlocks> & s : ss){
        for (uint8_t j = 0; j < s.pi.size(); ++j)
        {
            s.blockStarts[j] = (s.pi[j] - 1 == 0) ? 0 : s.chronBL[s.pi[j] - 2];
            s.blockEnds[j] = s.chronBL[s.pi[j] - 1];
        }

        for(uint8_t j = 0; j < s.pi.size(); ++j){
            s.revblockStarts[j] = (s.pi[j] == s.pi.size()) ? 0 : s.revChronBL[s.pi[j]];
            s.revblockEnds[j] = s.revChronBL[s.pi[j] - 1];
        }
    }
}

template<typename TVector>
inline void filterCheckSortBlockLimits(TVector & r, TVector & l)
{
    sort(r.begin(), r.end());
    sort(l.rbegin(), l.rend());
/*
    if(len != r.back() || len != l.front())
    {
        std::cerr << "Search Schemes contain unexpected blocklengths\n";
        exit(1);
    }*/
    auto itr = r.end();
    auto itl = l.end();


    if(itr != unique(r.begin(), r.end()) || itl != unique(l.begin(), l.end())){
        std::cerr << "Search Schemes contain unexpected blocklengths\n";
        std::cerr << "r: \n";
        for(int i = 0; i < r.size(); ++i){
            std::cerr << r[i] << "\t";
        }
        std::cerr << "\nl: \n";
        for(int i = 0; i < l.size(); ++i){
            std::cerr << l[i] << "\t";
        }
        std::cerr << "\n";
        exit(0);
    }

    r.pop_back();
    l.erase(l.begin());
}

inline auto loadBlockLimits(uint8_t se, uint32_t const len)
{
    std::vector<uint32_t> r(1, 0);
    std::vector<uint32_t> l(1, 0);
    switch (se)
    {
        case 0:
        {
            auto scheme = OptimalSearchSchemes<0, 0>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
            _optimalSearchSchemeComputeChronBlocklength(scheme);
            auto s = scheme[0];
            for(uint16_t i = 0; i < s.pi.size(); ++i)
            {
                r.push_back(s.chronBL[i]);
                l.push_back(s.revChronBL[i]);
            }
            break;
        }
        case 1:
        {
            auto scheme = OptimalSearchSchemes<0, 1>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
            _optimalSearchSchemeComputeChronBlocklength(scheme);
            auto s = scheme[0];
            for(uint16_t i = 0; i < s.pi.size(); ++i)
            {
                r.push_back(s.chronBL[i]);
                l.push_back(s.revChronBL[i]);
            }
            break;
        }
        case 2:
        {
            auto scheme = OptimalSearchSchemes<0, 2>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
            _optimalSearchSchemeComputeChronBlocklength(scheme);
            auto s = scheme[0];

            for(uint16_t i = 0; i < s.pi.size(); ++i)
            {
                r.push_back(s.chronBL[i]);
                l.push_back(s.revChronBL[i]);
            }
            break;
        }
        case 3:
        {
            auto scheme = OptimalSearchSchemes<0, 3>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
            _optimalSearchSchemeComputeChronBlocklength(scheme);
            auto s = scheme[0];
            for(uint16_t i = 0; i < s.pi.size(); ++i)
            {
                r.push_back(s.chronBL[i]);
                l.push_back(s.revChronBL[i]);
            }
            break;
        }
        case 4:
        {
            auto scheme = OptimalSearchSchemes<0, 4>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
            _optimalSearchSchemeComputeChronBlocklength(scheme);
            auto s = scheme[0];
            for(uint16_t i = 0; i < s.pi.size(); ++i)
            {
                r.push_back(s.chronBL[i]);
                l.push_back(s.revChronBL[i]);
            }
            break;
        }
        default: std::cerr << "E = " << (int)se << " not yet supported.\n";
                exit(1);
    }

    filterCheckSortBlockLimits(r, l);

    return std::make_pair(r, l);
}

template<size_t nbrBlocks, size_t N>
inline auto loadBlockLimits(std::array<OptimalSearch<nbrBlocks>, N> const & ss)
{
    std::vector<uint32_t> r(1, 0);
    std::vector<uint32_t> l(1, 0);
    auto s = ss[0];
    for(uint32_t i = 0; i < s.pi.size(); ++i)
    {
        r.push_back(s.chronBL[i]);
        l.push_back(s.revChronBL[i]);
    }
    filterCheckSortBlockLimits(r, l);

    return std::make_pair(r, l);
}

template<typename TContext,
         size_t nbrBlocks, size_t N,
         typename TBitvectorPair>
inline void linkBitvectors(TContext & ossContext,
                    std::array<OptimalSearch<nbrBlocks>, N> ss,
                    std::vector<TBitvectorPair > & bitvectors,
                    std::vector<TBitvectorPair * > & bit_filtered)
{
    if(bitvectors.empty())
        return;

    _optimalSearchSchemeComputeFixedBlocklength(ss, ossContext.readLength);
    _optimalSearchSchemeComputeChronBlocklength(ss);


    auto blocklengths = loadBlockLimits(ss);
    std::vector<uint32_t> shift_r = blocklengths.first;
    std::vector<uint32_t> shift_l = blocklengths.second;

    std::vector<std::pair<uint32_t, bool> > & meta = ossContext.bitvectorsMeta;
    std::cout << "SearchScheme has: " << ss[0].pi.size() << " parts" << "\n";

    //test blocklengths of search scheme
    for(uint16_t i = 0; i < shift_r.size() - 1; ++i){
        if(shift_r[i] > shift_r[i + 1] || shift_l[i] > shift_l[i]){
            std::cerr << "blocklengths not in the right order" << "\n";
            exit(0);
        }
    }
    for(int16_t i = 0; i < shift_r.size(); ++i){
        for(uint16_t j = 0; j < meta.size(); ++j){
            if(shift_r[i] == meta[j].first && meta[j].second)
            {
                std::cout << "left anchored with shift: " << meta[j].first << "\n";
                bit_filtered.push_back(& bitvectors[j]);
                break;
            }
        }
    }

    for(int32_t i = shift_l.size() - 1; i >= 0; --i){
        for(uint16_t j = 0; j < meta.size(); ++j){
            if(shift_l[i] == meta[j].first && !meta[j].second)
            {
                std::cout << "right anchored with shift: " << meta[j].first << "\n";
                bit_filtered.push_back(& bitvectors[j]);
                break;
            }
        }
    }

    if((shift_r.size() + shift_l.size()) != bit_filtered.size()){
        std::cerr << "Unable to load appropiate bitvectors for scheme with " << ss[0].pi.size() << " pieces\n";
        exit(0);
    }
}

template <typename TContigsSize, typename TContigsLen, typename TText>
auto getSeqLengths(TText & text){

    std::vector<TContigsLen> sl;
    sl.push_back(0);
    for(TContigsSize i = 0; i < length(text); ++i){
        sl.push_back(seqan::length(text[i]) + sl.back());
    }
    return sl;
}


enum class ReturnCode {
	NOMAPPABILITY, DIRECTSEARCH, COMPMAPPABLE, ONEDIRECTION, MAPPABLE, FINISHED, UNIDIRECTIONAL, SUSPECTUNIDIRECTIONAL, FILTER, ERROR
};
/*
 * wrong type !!!! neede TSALength
template <typename TVector, typename TVSupport>
inline void getConsOnes(std::vector<std::pair<TVector, TVSupport>> & bitvectors,
                 Pair<uint8_t, Pair<uint32_t, uint32_t>> & inside_bit_interval,
                 int const intervalsize,
                 std::vector<std::pair<uint32_t, uint32_t>> & consOnesOutput);*/


/*
template <typename TText, typename TSpec, typename TConfig>
inline bool open(Index<TText, BidirectionalIndex<FMIndex<TSpec, TConfig> > > & index, const char * fileName, int openMode)
{
    String<char> name;

    {
        name = fileName;    append(name, ".txt");
        if (!open(getFibre(index.fwd, FibreText()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".sa");
        if (!open(getFibre(index.fwd, FibreSA()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".lf");
        if (!open(getFibre(index.fwd, FibreLF()), toCString(name), openMode)) return false;

        setFibre(getFibre(index.fwd, FibreSA()), getFibre(index.fwd, FibreLF()), FibreLF());
    }

    {
        // name = fileName;    append(name, ".rev.txt");
        // if (!open(getFibre(index.rev, FibreText()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".rev.sa");
        if (!open(getFibre(index.rev, FibreSA()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".rev.lf");
        if (!open(getFibre(index.rev, FibreLF()), toCString(name), openMode)) return false;

        setFibre(getFibre(index.rev, FibreSA()), getFibre(index.rev, FibreLF()), FibreLF());
    }

    return true;
}

template <typename TText, typename TSpec, typename TConfig>
inline bool save(Index<TText, BidirectionalIndex<FMIndex<TSpec, TConfig> > > const & index, const char * fileName, int openMode)
{
    String<char> name;

    {
        name = fileName;    append(name, ".txt");
        if (!save(getFibre(index.fwd, FibreText()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".sa");
        if (!save(getFibre(index.fwd, FibreSA()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".lf");
        if (!save(getFibre(index.fwd, FibreLF()), toCString(name), openMode)) return false;
    }

    {
        // name = fileName;    append(name, ".rev.txt");
        // if (!save(getFibre(index.rev, FibreText()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".rev.sa");
        if (!save(getFibre(index.rev, FibreSA()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".rev.lf");
        if (!save(getFibre(index.rev, FibreLF()), toCString(name), openMode)) return false;
    }

    return true;
}
*/

}
/*
struct SearchParams
{
    unsigned length;
    unsigned overlap;
    unsigned threads;
    // bool indels;
    static constexpr bool outputProgress = false;
};

std::string mytime()
{
    auto r = time(nullptr);
    auto c = ctime(&r);
    std::string buf(c);
    buf.insert(0, "[");
    buf.append("] ");
    buf.erase(remove(buf.begin(), buf.end(), '\n'), buf.end());
    return buf;
}


*/

template <typename TChar, typename TSeqLengths, typename TConfig, typename TMappVector>
void resetLimits(seqan::String<TChar, TConfig> const &, TSeqLengths const & , TMappVector const &, unsigned const)
{ }

template <typename TString, typename TSeqLengths, typename TConfig, typename TMappVector>
void resetLimits(seqan::StringSet<TString, TConfig> const & text, TSeqLengths const & sequenceLengths, TMappVector & c, unsigned const length)
{
//     auto const & limits = seqan::stringSetLimits(text);
    for (unsigned i = 1; i < sequenceLengths.size() - 1; ++i) //seqan::length(limits)
    {
        for (unsigned j = 1; j < length; ++j)
        {
            c[sequenceLengths[i] - j] = 0;
        }
    }
}

template <bool outputProgress>
inline void printProgress(uint64_t &, uint64_t const, uint64_t const);

template <>
inline void printProgress<false>(uint64_t &, uint64_t const, uint64_t const)
{ }

template <>
inline void printProgress<true>(uint64_t & progress_count, uint64_t const progress_step, uint64_t const progress_max)
{
    #pragma omp atomic
    ++progress_count; // TODO: necessary?
    if (omp_get_thread_num() == 0 && (progress_count & progress_step) == 0)
    {
        float progress = static_cast<float>(progress_count)/progress_max;
        std::cout << "\rProgress: " << (truncf(progress*10000)/100) << "%" << std::flush;
    }
}

template <bool outputProgress>
inline void initProgress(uint64_t &, uint64_t &, uint64_t &, uint64_t const, uint64_t const);

template <>
inline void initProgress<false>(uint64_t & progress_count, uint64_t & progress_step, uint64_t & progress_max, uint64_t const step_size, uint64_t const max_i)
{ }

template <>
inline void initProgress<true>(uint64_t & progress_count, uint64_t & progress_step, uint64_t & progress_max, uint64_t const step_size, uint64_t const max_i)
{
    progress_count = 0;
    progress_max = (max_i + step_size - 1) / step_size; // = ceil(max_i / step_size)
    progress_step = 2;
    uint64_t progress_step_tmp = std::max(progress_max/10000, static_cast<uint64_t>(1));
    while (progress_step < progress_step_tmp)
        progress_step <<= 1;
    --progress_step;
}


// template <typename TString, typename TPos, typename TLength>
// auto extractNeedle(TString const & text, TPos const pos, TLength const len)
// {
//     using namespace seqan;
//     return infix(text, pos, pos + len);
// }
//
// template <typename TString, typename TOwner, typename TPos, typename TLength>
// auto extractNeedle(seqan::StringSet<TString, TOwner> const & text, TPos const pos, TLength const len)
// {
//     using namespace seqan;
//     auto const & limits = stringSetLimits(text);
//     Pair<unsigned, unsigned> res;
//     posLocalize(res, pos, limits);
//     if (res.i2 + len > length(text[res.i1]))
//     {
//         TString r(suffix(text[res.i1], res.i2));
//         r += prefix(text[res.i1 + 2], length(text[res.i1]) - res.i2);
//         return r;
//     }
//     else
//     {
//         return infix(text[res.i1], res.i2, res.i2 + len);
//     }
// }

// template <typename T>
// struct is_int_vector
// {
//     static constexpr bool value = false;
// };
//
// template <uint8_t width_t>
// struct is_int_vector<sdsl::int_vector<width_t> >
// {
//     static constexpr bool value = true;
// };
//
// template <typename T>
// constexpr bool is_int_vector<T>::value;
//
// template <uint8_t width_t>
// constexpr bool is_int_vector<sdsl::int_vector<width_t> >::value;
#endif

