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
#include <random>
#include <benchmark/benchmark.h>
#include <seqan/kmer.h>

using namespace seqan;

uint64_t ipow(uint64_t base, uint64_t exp)
{
    uint64_t result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

template <typename TAlphabet>
static void addKmer_IBF(benchmark::State& state)
{
    KmerFilter<TAlphabet, InterleavedBloomFilter> ibf (state.range(0), state.range(3), state.range(1), state.range(2));
    std::mt19937 RandomNumber;
    String<TAlphabet> kmer("");
    for (uint8_t i = 0; i < state.range(1); ++i)
        appendValue(kmer, TAlphabet(RandomNumber() % ValueSize<TAlphabet>::VALUE));
    for (auto _ : state)
        addKmer(ibf, kmer, 0);
}

template <typename TAlphabet>
static void whichBins_IBF(benchmark::State& state)
{
    KmerFilter<TAlphabet, InterleavedBloomFilter> ibf (state.range(0), state.range(3), state.range(1), state.range(2));
    std::mt19937 RandomNumber;
    String<TAlphabet> kmer("");
    for (uint8_t i = 0; i < state.range(1); ++i)
        appendValue(kmer, TAlphabet(RandomNumber() % ValueSize<TAlphabet>::VALUE));
    for (auto _ : state)
        whichBins(ibf, kmer, 0);
}

template <typename TAlphabet>
static void addKmer_DA(benchmark::State& state)
{
    KmerFilter<TAlphabet, DirectAddressing> da (state.range(0), state.range(1));
    std::mt19937 RandomNumber;
    String<TAlphabet> kmer("");
    for (uint8_t i = 0; i < state.range(1); ++i)
        appendValue(kmer, TAlphabet(RandomNumber() % ValueSize<TAlphabet>::VALUE));
    for (auto _ : state)
        addKmer(da, kmer, 0);
}

template <typename TAlphabet>
static void whichBins_DA(benchmark::State& state)
{
    KmerFilter<TAlphabet, DirectAddressing> da (state.range(0), state.range(1));
    std::mt19937 RandomNumber;
    String<TAlphabet> kmer("");
    for (uint8_t i = 0; i < state.range(1); ++i)
        appendValue(kmer, TAlphabet(RandomNumber() % ValueSize<TAlphabet>::VALUE));
    for (auto _ : state)
        whichBins(da, kmer, 0);
}

static void IBFArguments(benchmark::internal::Benchmark* b)
{
    for (int32_t binNo = 1; binNo <= 8192; binNo *= 2)
    {
        if (binNo > 1 && binNo < 64)
            continue;
        for (int32_t k = 20; k <= 20; ++k)
        {
            for (int32_t bits = 1<<16; bits <= 1<<21; bits <<= 1 )
            {
                for (int32_t hashNo = 3; hashNo < 4; ++hashNo)
                {
                    b->Args({binNo, k, bits+256, hashNo});
                }
            }
        }
    }
}

static void DAArguments(benchmark::internal::Benchmark* b)
{
    for (int32_t binNo = 1; binNo <= 128; binNo *= 2)
    {
        for (int32_t k = 12; k <= 13; ++k)
        {
            b->Args({binNo, k});
        }
    }
}

// 1<<16 = 65536

// 0=bins 1=k 2=bits 3=noOfHashes

BENCHMARK_TEMPLATE(addKmer_IBF, Dna)->Apply(IBFArguments);
BENCHMARK_TEMPLATE(addKmer_DA, Dna)->Apply(DAArguments);
BENCHMARK_TEMPLATE(whichBins_IBF, Dna)->Apply(IBFArguments);
BENCHMARK_TEMPLATE(whichBins_DA, Dna)->Apply(DAArguments);

BENCHMARK_MAIN();
