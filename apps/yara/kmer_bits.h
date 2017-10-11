// ==========================================================================
//                                 misc_functions.hpp
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

#ifdef YARA_DIST_NUM_OF_BINS
    uint32_t const NUM_OF_BINS = YARA_DIST_NUM_OF_BINS;
#else
    uint32_t const NUM_OF_BINS = 10;
#endif

#ifdef YARA_DIST_K_MER_LENGTH
    uint32_t const K_MER_LENGTH = YARA_DIST_K_MER_LENGTH;
#else
    uint32_t const K_MER_LENGTH = 10;
#endif

#ifdef YARA_DIST_INTIAL_ELEMENT_SIZE
uint64_t const INTIAL_ELEMENT_SIZE = YARA_DIST_INTIAL_ELEMENT_SIZE;
#else
uint64_t const INTIAL_ELEMENT_SIZE = 548573;
#endif

using namespace seqan;

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{

    CharString          fm_index_path;
    CharString          output_path;
    uint32_t            num_threads;

    uint64_t            contigsSize;
    uint64_t            contigsMaxLength;
    uint64_t            contigsSum;


    AppOptions() :
    num_threads(4)
    {}
};

inline void bits_to_dna(DnaString & result, uint32_t i, bool b0, bool b1)
{
    if(b0 && b1)
        result[i] = 'T';
    else if(b0)
        result[i] = 'G';
    else if(b1)
        result[i] = 'C';
}

inline void bit_vector_to_kmer(DnaString & result, std::vector<bool> bit_v)
{
    uint32_t i = 0;
    for (size_t current_bit = 0;  current_bit < bit_v.size(); current_bit += 2)
    {
        //        std::cout << result << "\t" << bit_v[current_bit] << ", " << bit_v[current_bit+1]  << "\n";
        bits_to_dna(result, i++, bit_v[current_bit], bit_v[current_bit+1]);
    }
}


inline Dna bits_to_dna(bool b0, bool b1)
{
    std::bitset<2> bits(0);

    if(b0)
        bits.set(1);
    if(b1)
        bits.set(0);

    return (Dna)(bits.to_ulong());
}



inline uint32_t log2_int32(uint32_t const int_num)
{
    uint32_t result = 0;
    std::bitset<32> num_in_bits(int_num);
    while (!num_in_bits.test(result))
        ++result;
    return result;
}


inline DnaString bit_vector_to_kmer(std::vector<bool> bit_v)
{
    size_t bit_v_l  = bit_v.size();

    DnaString result = "";

    for (size_t current_bit = 0;  current_bit < bit_v_l; current_bit += 2)
    {
//        std::cout << result << "\t" << bit_v[current_bit] << ", " << bit_v[current_bit+1]  << "\n";
        append(result, bits_to_dna(bit_v[current_bit], bit_v[current_bit+1]) );
    }
    return result;
}

inline bool next_bit_vector(std::vector<bool> && bit_v, uint32_t base)
{
    for (uint32_t i = bit_v.size(); i > base; --i)
    {
        bit_v[i-1] = !bit_v[i-1];
        if (bit_v[i-1] == true)
        {
            return true;
        }
    }
    return false;
}

inline bool next_bit_vector(std::vector<bool> && bit_v, uint32_t inc, uint32_t base)
{
    for (uint32_t i = bit_v.size() - inc + 1; i < bit_v.size(); ++i)
        bit_v[i] = true;
    return next_bit_vector(std::forward<std::vector<bool> >(bit_v), base);
}


std::vector<bool> get_start_bit_vector(uint32_t const & len, uint32_t curr, uint32_t base)
{
    std::bitset<32> curr_bits(curr);

//    std::cout << curr_bits << "\t" << base << "\n";
//
    std::vector<bool> result(len, 0);

     for (uint32_t i = base; i > 0; --i)
        result[i-1] = curr_bits.test(base-i);

    return result;
}



std::vector<bool> kmer_to_bitset(DnaString kmer)
{
    size_t k = length(kmer);
    std::vector<bool> bit_v;
    bit_v.reserve(k);

    for(auto c : kmer)
    {
        std::bitset<2> bits((int)c);
//        std::cout << kmer << "\t" << c << "\t " << x << "\t " << bits  << "\n";
        bit_v.push_back(bits.test(1));
        bit_v.push_back(bits.test(0));
    }
    return bit_v;
}

