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
// Author:  Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
//          Enrico Seiler <enrico.seiler@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_KMER_KMER_IBF_H_
#define INCLUDE_SEQAN_KMER_KMER_IBF_H_

// --------------------------------------------------------------------------
// Class KmerFilter using an interleaved bloom filter
// --------------------------------------------------------------------------
namespace seqan{

template<typename TValue>
class KmerFilter<TValue, InterleavedBloomFilter>
{
public:
    typedef typename Value<KmerFilter>::Type    THValue;
    THValue    noOfBins;
    THValue    noOfHashFunc;
    THValue    kmerSize;
    THValue    noOfBits;
    THValue    noOfBlocks;

    std::vector<THValue>   preCalcValues;
    static const THValue   shiftValue = 27;
    static const THValue   seedValue = 0x90b45d39fb6da1fa;

    sdsl::bit_vector                    filterVector;

    static const uint32_t               filterMetadataSize{256};

    typedef Shape<TValue, SimpleShape>  TShape;

    // Default constructor
    KmerFilter():
        noOfBins(0),
        noOfHashFunc(0),
        kmerSize(0),
        noOfBits(0),
        filterVector(sdsl::bit_vector(0, 0)) {}

    // Default constructor
    KmerFilter(THValue n_bins, THValue n_hash_func, THValue kmer_size, THValue vec_size):
        noOfBins(n_bins),
        noOfHashFunc(n_hash_func),
        kmerSize(kmer_size),
        noOfBits(vec_size),
        filterVector(sdsl::bit_vector(vec_size, 0))
    {
            init();
    }

    // Copy constructor
    KmerFilter(KmerFilter<TValue, InterleavedBloomFilter> const & other)
    {
        *this = other;
    }

    // Copy assignment
    KmerFilter<TValue, InterleavedBloomFilter> & operator=(KmerFilter<TValue, InterleavedBloomFilter> const & other)
    {
        noOfBins = other.noOfBins;
        noOfHashFunc = other.noOfHashFunc;
        kmerSize = other.kmerSize;
        noOfBits = other.noOfBits;
        filterVector = other.filterVector;
        init();
        return *this;
    }

    // Move constrcutor
    KmerFilter(KmerFilter<TValue, InterleavedBloomFilter> && other)
    {
        *this = std::move(other);
    }

    // Move assignment
    KmerFilter<TValue, InterleavedBloomFilter> & operator=(KmerFilter<TValue, InterleavedBloomFilter> && other)
    {
        noOfBins = std::move(other.noOfBins);
        noOfHashFunc = std::move(other.noOfHashFunc);
        kmerSize = std::move(other.kmerSize);
        noOfBits = std::move(other.noOfBits);
        filterVector = std::move(other.filterVector);
        init();
        return *this;
    }

    // Destructor
    ~KmerFilter<TValue, InterleavedBloomFilter>() = default;

    template<typename TInt>
    void clearBins(std::vector<THValue> const & bins, TInt&& threads)
    {
        std::vector<std::future<void>> tasks;

        // We have so many blocks that we want to distribute to so many threads
        uint64_t batchSize = noOfBlocks / threads;
        if(batchSize * threads < noOfBlocks) ++batchSize;

        for (uint32_t taskNo = 0; taskNo < threads; ++taskNo)
        {
            // hashBlock is the number of the block the thread will work on. Each block contains binNo bits that
            // represent the individual bins. Each thread has to work on batchSize blocks. We can get the position in
            // our filterVector by multiplying the hashBlock with noOfBins. Then we just need to add the respective
            // binNo. We have to make sure that the vecPos we generate is not out of bounds, only the case in the last
            // thread if the blocks could not be evenly distributed, and that we do not clear a bin that is assigned to
            // another thread.
            tasks.emplace_back(std::async([=] {
                for (uint64_t hashBlock=taskNo*batchSize;
                    hashBlock < noOfBlocks && hashBlock < (taskNo +1) * batchSize;
                    ++hashBlock)
                {
                    uint64_t vecPos = hashBlock * noOfBins;
                    for(uint32_t binNo : bins)
                    {
                        filterVector[vecPos + binNo] = false;
                    }
                }
            }));
        }
        for (auto &&task : tasks)
        {
            task.get();
        }
    }

    template<typename TString>
    void whichBins(std::vector<uint64_t> & counts, TString const & text) const
    {
        uint8_t possible = length(text) - kmerSize + 1;

        std::vector<uint64_t> kmerHashes(possible, 0);

        TShape kmerShape;
        resize(kmerShape, kmerSize);
        hashInit(kmerShape, begin(text));
        auto it = begin(text);
        for (uint32_t i = 0; i < possible; ++i)
        {
            kmerHashes[i] = hashNext(kmerShape, it);
            ++it;
        }

        for (uint64_t kmerHash : kmerHashes)
        {
            std::vector<uint64_t> vecIndices = preCalcValues;
            for(uint8_t i = 0; i < noOfHashFunc ; i++)
            {
                vecIndices[i] *= kmerHash;
                // Move to first bit representing the hash kmerHash for bin 0,
                // the next bit would be for bin 1, and so on
                getHashValue(vecIndices[i]);
            }

            // get_int(idx, len) returns the integer value of the binary string of length len starting at position idx.
            // I.e. len+idx-1|_______|idx, Vector is right to left.
            uint64_t tmp = filterVector.get_int(vecIndices[0], noOfBins);

            for(uint8_t i = 1; i < noOfHashFunc;  i++)
            {
                tmp &= filterVector.get_int(vecIndices[i], noOfBins);
            }

            uint64_t binNo = 0;

            // Magic thingie that prevents segfault
            // TODO Learn more about magic thingie
            if (tmp ^ (1ULL<<(noOfBins-1)))
            {
                // As long as any bit is set
                while (tmp > 0)
                {
                    // sdsl::bits::lo calculates the position of the rightmost 1-bit in the 64bit integer x if it exists.
                    // For example, for 8 = 1000 it would return 3
                    uint64_t step = sdsl::bits::lo(tmp);
                    // Adjust our bins
                    binNo += step;
                    // Remove up to next 1
                    ++step;
                    tmp >>= step;
                    // Count
                    ++counts[binNo];
                    // ++binNo because step is 0-based, e.g., if we had a hit with the next bit we would otherwise count it
                    // for binNo=+ 0
                    ++binNo;
                }
            }
            else
            {
                ++counts[binNo + noOfBins - 1];
            }
        }
    }

    template<typename TString, typename TInt>
    inline void whichBins(std::vector<bool> & selected, TString const & text, TInt && threshold) const
    {
        std::vector<uint64_t> counts(noOfBins, 0);
        whichBins(counts, text);
        for(uint32_t binNo=0; binNo < noOfBins; ++binNo)
        {
            if(counts[binNo] >= threshold)
                selected[binNo] = true;
        }
    }

    inline void getHashValue(uint64_t & vecIndex) const
    {
        // We do something
        vecIndex ^= vecIndex >> shiftValue;
        // Bring it back into our vector range (noOfBlocks = possible hash values)
        vecIndex %= noOfBlocks;
        // Since each block needs binNo bits, we multiply to get the correct location
        vecIndex *= noOfBins;
    }

    template<typename TString, typename TInt>
    inline void addKmer(TString const & text, TInt && binNo)
    {

        TShape kmerShape;
        resize(kmerShape, kmerSize);
        hashInit(kmerShape, begin(text));

        for (uint64_t i = 0; i < length(text) - length(kmerShape) + 1; ++i)
        {
            uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);

            for(uint8_t i = 0; i < noOfHashFunc ; i++)
            {
                uint64_t vecIndex = preCalcValues[i] * kmerHash;
                getHashValue(vecIndex);
                vecIndex += binNo;
                filterVector[vecIndex] = 1;
            }
        }
    }

    inline void init()
    {
        noOfBlocks = (noOfBits - filterMetadataSize) / noOfBins;

        preCalcValues.resize(noOfHashFunc);
        for(uint64_t i = 0; i < noOfHashFunc ; i++)
            preCalcValues[i] = i ^  (kmerSize * seedValue);
    }
};
}

#endif  // INCLUDE_SEQAN_KMER_KMER_IBF_H_
