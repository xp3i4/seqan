// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author:  Enrico Seiler <enrico.seiler@fu-berlin.de>
// ==========================================================================

#include <chrono>
#include <numeric>
#include <random>

#include <seqan/kmer.h>

// static global const values needed for bloom_filter
static const uint32_t filterMetadataSize = 256;
static const uint8_t INT_WIDTH = 0x40;
#include "../../../../apps/yara/bloom_filter.h"

using namespace seqan;

int main()
{
    std::mt19937 rng;
    StringSet<DnaString> input;
    std::vector<int64_t> ibfTime;
    std::vector<int64_t> sbfTime;
    std::vector<int64_t> daTime;
    reserve(input, 1000000);
    for (uint32_t seqNo = 0; seqNo < 1000000; ++seqNo)
    {
        DnaString tmp;
        for (uint16_t i = 0; i < 12; ++i)
        {
            appendValue(tmp, Dna(rng() % 4));
        }
        appendValue(input, tmp);
    }

    std::cout << "====================================================================\n"
              << "====================== Benchmarking addKmer() ======================\n"
              << "====================================================================\n";
    for (uint16_t r = 0; r < 10; r++)
    {
        KmerFilter<Dna, InterleavedBloomFilter> ibf (32, 3, 12, 16777472);
        SeqAnBloomFilter<DnaString> sbf (32, 3, 12, 16777472);
        KmerFilter<Dna, DirectAddressing> da (32, 12);

        auto start = std::chrono::high_resolution_clock::now();
        for(uint32_t i = 0; i < 1000000; ++i)
        {
            ibf.addKmer(input[i], i % 32);
        }
        auto elapsed = std::chrono::high_resolution_clock::now() - start;
        ibfTime.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());

        start = std::chrono::high_resolution_clock::now();
        for(uint32_t i = 0; i < 1000000; ++i)
        {
            sbf.addKmers(input[i], i % 32);
        }
        elapsed = std::chrono::high_resolution_clock::now() - start;
        sbfTime.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());

        start = std::chrono::high_resolution_clock::now();
        for(uint32_t i = 0; i < 1000000; ++i)
        {
            da.addKmer(input[i], i % 32);
        }
        elapsed = std::chrono::high_resolution_clock::now() - start;
        daTime.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
    }

    auto ibfAvg = accumulate(ibfTime.begin(), ibfTime.end(), 0)/ibfTime.size();
    auto sbfAvg = accumulate(sbfTime.begin(), sbfTime.end(), 0)/sbfTime.size();
    auto daAvg = accumulate(daTime.begin(), daTime.end(), 0)/daTime.size();

    std::cout << "Average InterleavedBloomFilter: " << ibfAvg << " ms.\n";
    std::cout << "Average SeqAnBloomFilter: " << sbfAvg << " ms.\n";
    std::cout << "Average DirectAddressing: " << daAvg << " ms.\n";

    ibfTime.clear();
    ibfTime.shrink_to_fit();
    sbfTime.clear();
    sbfTime.shrink_to_fit();
    daTime.clear();
    daTime.shrink_to_fit();


    std::cout << "====================================================================\n"
              << "===================== Benchmarking whichBins() =====================\n"
              << "====================================================================\n";

    KmerFilter<Dna, InterleavedBloomFilter> ibf (32, 3, 12, 16777472);
    SeqAnBloomFilter<DnaString> sbf (32, 3, 12, 16777472);
    KmerFilter<Dna, DirectAddressing> da (32, 12);

    for(uint32_t i = 0; i < 1000000; ++i)
    {
      ibf.addKmer(input[i], i % 32);
      sbf.addKmers(input[i], i % 32);
      da.addKmer(input[i], i % 32);
    }

    for (uint16_t r = 0; r < 10; r++)
    {
        auto start = std::chrono::high_resolution_clock::now();
        for(uint32_t i = 0; i < 1000000; ++i)
        {
            (void) whichBins(ibf, input[i], 1);
        }
        auto elapsed = std::chrono::high_resolution_clock::now() - start;
        ibfTime.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());

        start = std::chrono::high_resolution_clock::now();
        for(uint32_t i = 0; i < 1000000; ++i)
        {
            (void) sbf.whichBins(input[i], 1);
        }
        elapsed = std::chrono::high_resolution_clock::now() - start;
        sbfTime.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());

        start = std::chrono::high_resolution_clock::now();
        for(uint32_t i = 0; i < 1000000; ++i)
        {
            (void) whichBins(da, input[i], 1);
        }
        elapsed = std::chrono::high_resolution_clock::now() - start;
        daTime.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count());
    }

    ibfAvg = accumulate(ibfTime.begin(), ibfTime.end(), 0)/ibfTime.size();
    sbfAvg = accumulate(sbfTime.begin(), sbfTime.end(), 0)/sbfTime.size();
    daAvg = accumulate(daTime.begin(), daTime.end(), 0)/daTime.size();

    std::cout << "Average InterleavedBloomFilter: " << ibfAvg << " ms.\n";
    std::cout << "Average SeqAnBloomFilter: " << sbfAvg << " ms.\n";
    std::cout << "Average DirectAddressing: " << daAvg << " ms.\n";
}
