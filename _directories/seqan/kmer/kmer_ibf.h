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

// --------------------------------------------------------------------------
// Class KmerFilter using Interleaved Bloom Filter
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

    THValue    binIntWidth;
    THValue    blockBitSize;
    THValue    noOfHashPos;

    sdsl::bit_vector                    filterVector;

    static const uint32_t               filterMetadataSize{256};
    static const uint8_t                INT_WIDTH{0x40};

    std::vector<THValue>                preCalcValues;
    THValue const                       shiftValue{27};
    THValue const                       seedValue{0x90b45d39fb6da1fa};

    typedef Shape<TValue, SimpleShape>  TShape;

    KmerFilter():
        noOfBins(0),
        noOfHashFunc(0),
        kmerSize(0),
        noOfBits(0),
        filterVector(sdsl::bit_vector(0, 0)) {}

    KmerFilter(THValue n_bins, THValue n_hash_func, THValue kmer_size, THValue vec_size):
        noOfBins(n_bins),
        noOfHashFunc(n_hash_func),
        kmerSize(kmer_size),
        noOfBits(vec_size),
        filterVector(sdsl::bit_vector(vec_size, 0))
    {
            init(*this);
            seed();
            std::cout << noOfBins << '\n';
            std::cout << noOfHashFunc << '\n';
            std::cout << kmerSize << '\n';
            std::cout << noOfBits << '\n';
            std::cout << binIntWidth << '\n';
            std::cout << blockBitSize << '\n';
            std::cout << noOfHashPos << '\n';
    }

    KmerFilter(KmerFilter<TValue, DirectAddressing> const & other)
    {
        *this = other;
    }

    KmerFilter<TValue, DirectAddressing> & operator=(KmerFilter<TValue, DirectAddressing> const & other)
    {
        noOfBins = other.noOfBins;
        noOfHashFunc = other.noOfHashFunc;
        kmerSize = other.kmerSize;
        noOfBits = other.noOfBits;
        filterVector = other.filterVector;
        init(*this);
        seed();
        return *this;
    }

    inline void seed()
    {
        preCalcValues.resize(noOfHashFunc);
        for(uint8_t i = 0; i < noOfHashFunc ; i++)
            preCalcValues[i] = i ^  (kmerSize * seedValue);
    }

    void clearBins(std::vector<TValue> & bins, THValue & threads)
    {
        std::vector<std::future<void>> tasks;

        uint64_t batchSize = noOfHashPos/threads;
        if(batchSize * threads < noOfHashPos) ++batchSize;

        for (uint32_t taskNo = 0; taskNo < threads; ++taskNo)
        {
            tasks.emplace_back(std::async([=] {
                for (uint64_t hashBlock=taskNo*batchSize; hashBlock < noOfHashPos && hashBlock < (taskNo +1) * batchSize; ++hashBlock)
                {
                    uint64_t vecPos = hashBlock * blockBitSize;
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
    inline void whichBins(std::vector<bool> & selected, TString const & text, THValue const & threshold) const
    {
        uint8_t possible = length(text) - kmerSize + 1;

        std::vector<uint8_t> counts(noOfBins, 0);
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
                getHashValue(vecIndices[i]);
            }
            uint32_t binNo = 0;
            for (uint8_t batchNo = 0; batchNo < binIntWidth; ++batchNo)
            {
                binNo = batchNo * INT_WIDTH;
                uint64_t tmp = filterVector.get_int(vecIndices[0], INT_WIDTH);
                for(uint8_t i = 1; i < noOfHashFunc;  i++)
                {
                    tmp &= _filterVector.get_int(vecIndices[i], INT_WIDTH);
                }

                if (tmp ^ (1ULL<<(INT_WIDTH-1)))
                {
                    while (tmp > 0)
                    {
                        uint64_t step = sdsl::bits::lo(tmp);
                        binNo += step;
                        ++step;
                        tmp >>= step;
                        ++counts[binNo];
                        ++binNo;
                    }
                }
                else
                {
                    ++counts[binNo + INT_WIDTH - 1];
                }
                for(uint8_t i = 0; i < noOfHashFunc ; i++)
                {
                    vecIndices[i] += INT_WIDTH;
                }
            }
        }

        for(uint32_t binNo=0; binNo < noOfBins; ++binNo)
        {
            if(counts[binNo] >= threshold)
                selected[binNo] = true;
        }
    }

    inline void insertKmer(THValue & kmerHash, THValue const & batchOffset)
    {
        for(uint8_t i = 0; i < noOfHashFunc ; i++)
        {
            uint64_t vecIndex = preCalcValues[i] * kmerHash;
            getHashValue(vecIndex);
            vecIndex += batchOffset;
            filterVector[vecIndex] = 1;
        }
    }

    template<typename TString>
    inline void addKmer(TString const & text, THValue const & binNo)
    {
        TShape kmerShape;
        resize(kmerShape, kmerSize);
        hashInit(kmerShape, begin(text));

        for (uint32_t i = 0; i < length(text) - length(kmerShape) + 1; ++i)
        {
            uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);
            insertKmer(kmerHash, binNo);
        }
    }
};
}
