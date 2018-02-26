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
// Class KmerFilter using direct addressing
// --------------------------------------------------------------------------
namespace seqan{

template<typename TValue>
class KmerFilter<TValue, DirectAddressing>
{
public:
    typedef typename Value<KmerFilter>::Type    THValue;
    THValue    noOfBins;
    THValue    kmerSize;
    THValue    noOfBits;
    THValue    noOfBlocks;

    sdsl::bit_vector                    filterVector;

    static const uint32_t               filterMetadataSize{256};
    THValue                noOfHashFunc{1};

    typedef Shape<TValue, SimpleShape>  TShape;

    KmerFilter():
        noOfBins(0),
        kmerSize(0)
        {}

    KmerFilter(THValue n_bins, THValue kmer_size):
        noOfBins(n_bins),
        kmerSize(kmer_size)
    {
        init();
    }

    KmerFilter(KmerFilter<TValue, DirectAddressing> const & other)
    {
        *this = other;
    }

    KmerFilter<TValue, DirectAddressing> & operator=(KmerFilter<TValue, DirectAddressing> const & other)
    {
        noOfBins = other.noOfBins;
        kmerSize = other.kmerSize;
        noOfBits = other.noOfBits;
        noOfBlocks = other.noOfBlocks;
        filterVector = other.filterVector;
        return *this;
    }

    THValue ipow(THValue base, THValue exp)
    {
        THValue result = 1;
        while (exp)
        {
            if (exp & 1)
                result *= base;
            exp >>= 1;
            base *= base;
        }
        return result;
    }

    template<typename TInt>
    void clearBins(std::vector<THValue> & bins, TInt&& threads)
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
            // Move to first bit representing the hash kmerHash for bin 0, the next bit would be for bin 1, and so on
            kmerHash *= noOfBins;

            // get_int(idx, len) returns the integer value of the binary string of length len starting at position idx.
            // I.e. len+idx-1|_______|idx, Vector is right to left.
            uint64_t tmp = filterVector.get_int(kmerHash, noOfBins);

            uint64_t binNo = 0;

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

    template<typename TString, typename TInt>
    inline void addKmer(TString const & text, TInt && binNo)
    {
        TShape kmerShape;
        resize(kmerShape, kmerSize);
        hashInit(kmerShape, begin(text));

        for (uint64_t i = 0; i < length(text) - length(kmerShape) + 1; ++i)
        {
            uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);
            filterVector[noOfBins * kmerHash + binNo] = 1;
        }
    }

    inline void init()
    {
        noOfBlocks = ipow(ValueSize<TValue>::VALUE, kmerSize);
        noOfBits = noOfBins * noOfBlocks + filterMetadataSize;
        std::cout << "Direct Addressing will need " << noOfBits << " bits." << '\n';
        filterVector = sdsl::bit_vector(noOfBits, 0);
    }
};
}
