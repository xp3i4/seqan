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
// Author: Enrico Seiler <enrico.seiler@fu-berlin.de>
// ==========================================================================
// kmer lib
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_KMER_KMER_BASE_H_
#define SEQAN_INCLUDE_SEQAN_KMER_KMER_BASE_H_

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Tag k-mer Filter Tags
// --------------------------------------------------------------------------

struct InterleavedBloomFilter_;
typedef Tag<InterleavedBloomFilter_> InterleavedBloomFilter;

struct DirectAddressing_;
typedef Tag<DirectAddressing_> DirectAddressing;

// --------------------------------------------------------------------------
// Class KmerFilter
// --------------------------------------------------------------------------

template<typename TValue = Dna, typename TSpec = DirectAddressing>
class KmerFilter;

// --------------------------------------------------------------------------
// Class KmerFilter using direct addressing
// --------------------------------------------------------------------------

template<typename TValue>
class KmerFilter<TValue, DirectAddressing>
{
public:
    typename Value<KmerFilter>::Type    noBins;
    typename Value<KmerFilter>::Type    kmerSize;
    typename Value<KmerFilter>::Type    noBits;
    sdsl::bit_vector                    filterVector;
    typedef Shape<TValue, SimpleShape>  TShape;

    KmerFilter():
        noBins(0),
        kmerSize(0),
        noBits(0),
        filterVector(sdsl::bit_vector(0, 0)) {}

    KmerFilter(n_bins, kmer_size, vec_size):
        noBins(n_bins),
        kmerSize(kmer_size),
        noBits(vec_size),
        filterVector(sdsl::bit_vector(vec_size, 0))
    {
            //TODO _init()
            _init();
    }

    KmerFilter(KmerFilter<TValue, DirectAddressing> const & other)
    {
        *this = other;
    }

    KmerFilter & operator=(KmerFilter const & other)
    {
        noBins = other.noBins;
        kmerSize = other.kmerSize;
        noBits = other.noBits;
        filterVector = other.filterVector;
        _init();
        return *this;
        // TODO what about the stuff set by _init?
    }

}
// TODO load from file? -> via open?
// TODO static const uint32_t filterMetadataSize = 256;
// TODO static const uint8_t INT_WIDTH = 0x40;


// --------------------------------------------------------------------------
// Class KmerFilter using Interleaved Bloom Filter
// --------------------------------------------------------------------------

template<typename TValue>
class KmerFilter<TValue, InterleavedBloomFilter>
{
public:
    typename Value<KmerFilter>::Type noBins;
    typename Value<KmerFilter>::Type noHash;
    typename Value<KmerFilter>::Type kmerSize;
    typename Value<KmerFilter>::Type vecSize;
    typedef Shape<TValue, SimpleShape> TShape;

    KmerFilter(n_bins, n_hash, kmer_size, vec_size):
        noBins(n_bins),
        noHash(n_hash),
        kmerSize(kmer_size),
        vecSize(vec_size)
    {

    }
}

// ==========================================================================
// Metafunctions
// ==========================================================================

template<typename TValue, typename TSpec>
struct Value<KmerFilter<TValue, TSpec> >
{
    typedef uint64_t Type;
}

// --------------------------------------------------------------------------
// Metafunction MetafunctionName
// --------------------------------------------------------------------------

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function functionName()
// --------------------------------------------------------------------------

}  // namespace seqan

#endif  // INCLUDE_SEQAN_BASIC_ITERATOR_BASE_H_
