// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
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

#ifndef APP_YARA_MISC_OPTIONS_DIS_H_
#define APP_YARA_MISC_OPTIONS_DIS_H_

using namespace seqan;

namespace seqan
{
    template<uint8_t KMER_SIZE, uint8_t N_HASH, uint32_t SIZE, typename TString=Dna5String>
    class SeqAnBloomFilter
    {
    public:

        typedef Shape<Dna5, UngappedShape<KMER_SIZE> > TShape;

        bool open(const char *fileName)
        {
            return seqan::open(_filterFile, fileName, OPEN_RDONLY);
        }

        SeqAnBloomFilter(){
            _initPreCalcValues();
        }

        SeqAnBloomFilter(const char *fileName)
        {
            openNewFile(fileName);
        }

        SeqAnBloomFilter(const char *fileName,
                         TString & text)
        {
            openNewFile(fileName);
            if(length(text) >= KMER_SIZE)
                addKmers(text);
        }

        void addKmers(TString const & text)
        {
            _addKmers(text);
        }

        bool containsNKmers(TString const & text, uint8_t const & threshold) const
        {
            TShape kmerShape;
            hashInit(kmerShape, begin(text));

            uint8_t left = threshold;
            uint8_t possible = length(text) - length(kmerShape) + 1;
            auto it = begin(text);
            uint64_t kmerHash = 0;
            for (uint8_t i = 0; left > 0 && (possible - i) > left ; ++i)
            {
                kmerHash = hashNext(kmerShape, it);
//                getNHashValues(hashValues, kmerHash);
                if(containsKmer(kmerHash))
                    --left;
                ++it;
            }
            return (left == 0);
        }


        bool containsNKmers(TString const & fwd, TString const & rev, uint32_t const & threshold) const
        {
            if(containsNKmers(fwd, threshold))
                return true;
            else if(containsNKmers(rev, threshold))
                return true;
            return false;
        }

        uint32_t kmerCount(TString const & text)
        {
            TShape kmerShape;
            hashInit(kmerShape, begin(text));
            uint32_t found = 0;
            for (unsigned i = 0; i < length(text) - length(kmerShape) + 1 ; ++i)
            {
                uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);
                std::vector<uint64_t> hashValues(N_HASH);
                getNHashValues(hashValues, kmerHash);
                if(containsKmer(kmerHash))
                    ++found;
            }
            return found ;
        }


    private:
        void getNHashValues(std::vector<uint64_t> & hashValues, uint64_t & kmerHash)
        {
            uint64_t tmp = 0;
            for(uint8_t i = 0; i < N_HASH ; i++)
            {
                tmp = kmerHash * (i ^ KMER_SIZE * _seedValue);
                tmp ^= tmp >> _shiftValue;
                hashValues[i] = tmp;
            }
        }

        void insertKmer(std::vector<uint64_t> & hashValues)
        {
            for(uint8_t i = 0; i < N_HASH ; i++)
                assignValue(_filterFile, (hashValues[i] % SIZE), true);
        }

        void insertKmer(uint64_t & kmerHash)
        {
            uint64_t tmp = 0;
            for(uint8_t i = 0; i < N_HASH ; i++)
            {
                tmp = kmerHash * (i ^ KMER_SIZE * _seedValue);
                tmp ^= tmp >> _shiftValue;
                assignValue(_filterFile, (tmp % SIZE), true);
            }
        }

        inline bool containsKmer(uint64_t & kmerHash) const
        {
            uint64_t tmp = 0;
            for(uint8_t i = 0; i < N_HASH ; i++)
            {
                tmp = kmerHash * (_preCalcValues[i]);
                tmp ^= tmp >> _shiftValue;
                if(!_filterFile[(tmp % SIZE)])
                    return false;
            }
            return true;
        }

        bool containsKmer(std::vector<uint64_t> & hashValues)
        {
            for(uint8_t i = 0; i < N_HASH ; i++)
                if(!_filterFile[(hashValues[i] % SIZE)])
                    return false;
            return true;
        }

        void _addKmers(TString const & text)
        {
            TShape kmerShape;
            hashInit(kmerShape, begin(text));

            for (uint32_t i = 0; i < length(text) - length(kmerShape) + 1; ++i)
            {
                uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);
                std::vector<uint64_t> hashValues(N_HASH);
                getNHashValues(hashValues, kmerHash);
                insertKmer(hashValues);
            }
        }

        inline void _initPreCalcValues()
        {
            for(uint8_t i = 0; i < N_HASH ; i++)
            {
                _preCalcValues.push_back(i ^ KMER_SIZE * _seedValue);
            }
        }
        void openNewFile(const char *fileName)
        {
            seqan::open(_filterFile, fileName, OPEN_RDWR | OPEN_CREATE);
            reserve(_filterFile, SIZE, Exact());
            resize(_filterFile, SIZE, false);
        }
        std::vector<uint64_t> _preCalcValues = {};
        uint64_t const _shiftValue = 27;
        uint64_t const _seedValue = 0x90b45d39fb6da1fa;
        String<bool, Packed<MMap<> > > _filterFile;
    };
}

// ============================================================================
// Functions
// ============================================================================
// ----------------------------------------------------------------------------
// Function appendFileName()
// ----------------------------------------------------------------------------
inline void appendFileName(CharString & target, CharString const & source, uint32_t const i)
{
    target = source;
    append(target, std::to_string(i));
}

inline void appendFileName(CharString & target, uint32_t const i)
{
    CharString source = target;
    appendFileName(target, source, i);
}

#endif  // #ifndef APP_YARA_MISC_OPTIONS_DIS_H_
