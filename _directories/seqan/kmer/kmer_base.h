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

//#ifndef SEQAN_INCLUDE_SEQAN_KMER_KMER_BASE_H_
//#define SEQAN_INCLUDE_SEQAN_KMER_KMER_BASE_H_
#include <sdsl/bit_vectors.hpp>
#include <valarray>
#include <algorithm>
#include <future>

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

// ==========================================================================
// Metafunctions
// ==========================================================================

template<typename TValue, typename TSpec>
struct Value<KmerFilter<TValue, TSpec> >
{
    typedef uint64_t Type;
};

// --------------------------------------------------------------------------
// Metafunction MetafunctionName
// --------------------------------------------------------------------------

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function functionName()
// --------------------------------------------------------------------------

template<typename TValue, typename TSpec, typename TString, typename TInt>
inline void addKmer(KmerFilter<TValue, TSpec> & me, TString const & kmer, TInt && binNo)
{
    me.addKmer(kmer, binNo);
}

template<typename TValue, typename TSpec, typename TInt1, typename TInt2>
inline void clearBins(KmerFilter<TValue, TSpec> &  me, std::vector<TInt1> & bins, TInt2&& threads)
{
    me.clearBins(bins, static_cast<uint64_t>(threads));
}

template<typename TValue, typename TSpec, typename TInt>
inline void addFastaFile(KmerFilter<TValue, TSpec> &  me, const char * fastaFile, TInt && binNo)
{
    CharString id;
    String<TValue> seq;

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, fastaFile))
    {
        CharString msg = "Unable to open contigs file: ";
        append(msg, CharString(fastaFile));
        throw toCString(msg);
    }
    while(!atEnd(seqFileIn))
    {
        readRecord(id, seq, seqFileIn);
        if(length(seq) < me.kmerSize)
            continue;
        addKmer(me, seq, binNo);
    }
    close(seqFileIn);
}

template<typename TValue, typename TSpec, typename TInt>
inline void whichBins(KmerFilter<TValue, TSpec> &  me, std::vector<bool> & selected, String<TValue> const & text, TInt threshold)
{
    me.whichBins(selected, text, threshold);
}

template<typename TValue, typename TSpec, typename TInt>
inline std::vector<bool> whichBins(KmerFilter<TValue, TSpec> &  me, String<TValue> const & text, TInt && threshold)
{
    std::vector<bool> selected(me.noOfBins, false);
    whichBins(me, selected, text, threshold);
    return selected;
}

template<typename TValue, typename TSpec>
inline typename Value<KmerFilter<TValue, TSpec> >::Type getNumberOfBins(KmerFilter<TValue, TSpec> &  me)
{
    return me.noBins;
}

template<typename TValue, typename TSpec>
inline typename Value<KmerFilter<TValue, TSpec> >::Type getKmerSize(KmerFilter<TValue, TSpec> &  me)
{
    return me.kmerSize;
}

template<typename TValue, typename TSpec>
inline void getMetadata(KmerFilter<TValue, TSpec> &  me)
{
    typedef typename Value<KmerFilter<TValue, TSpec> >::Type THValue;

    //-------------------------------------------------------------------
    //| n_bins | n_hash_func | kmer_size |              bf              |
    //-------------------------------------------------------------------
    me.noOfBits = me.filterVector.bit_size();

    THValue metadataStart = me.noOfBits - me.filterMetadataSize;
    me.noOfBins = me.filterVector.get_int(metadataStart);
    me.kmerSize = me.filterVector.get_int(metadataStart+128);
}

template<typename TValue, typename TSpec>
inline void setMetadata(KmerFilter<TValue, TSpec> &  me)
{
    typedef typename Value<KmerFilter<TValue, TSpec> >::Type THValue;

    //-------------------------------------------------------------------
    //| n_bins | n_hash_func | kmer_size |              bf              |
    //-------------------------------------------------------------------

    THValue metadataStart = me.noOfBits - me.filterMetadataSize;

    me.filterVector.set_int(metadataStart, me.noOfBins);
    me.filterVector.set_int(metadataStart+128, me.kmerSize);
}

template<typename TValue, typename TSpec, typename THValue>
inline double size(KmerFilter<TValue, TSpec> &  me)
{
    return sdsl::size_in_mega_bytes(me.filterVector);
}

template<typename TValue, typename TSpec>
inline bool store(KmerFilter<TValue, TSpec> &  me, const char * fileName)
{
    setMetadata(me);
    return sdsl::store_to_file(me.filterVector, fileName);
}

template<typename TValue, typename TSpec>
inline bool retrieve(KmerFilter<TValue, TSpec> &  me, const char * fileName)
{
    if (!sdsl::load_from_file(me.filterVector, fileName))
    {
        std::cerr << "File \"" << fileName << "\" could not be read." << std::endl;
        exit(1);
    }
    getMetadata(me);
    me.init();
    return true;
}

template<typename TValue, typename TSpec, typename TInt1, typename TInt2>
inline bool isBitSet(KmerFilter<TValue, TSpec> &  me, TInt1&& num, TInt2&& bit)
{
    return 1 == ( (num >> bit) & 1);
}

}  // namespace seqan

//#endif  // INCLUDE_SEQAN_BASIC_ITERATOR_BASE_H_
