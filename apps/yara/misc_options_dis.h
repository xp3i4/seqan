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
    static const uint8_t bitsPerChar = 0x08;
    static const unsigned char bitMask[0x08] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20,
        0x40, 0x80 };

    template<uint8_t BINS_SIZE, uint8_t KMER_SIZE, uint8_t N_HASH, uint32_t SIZE, typename TString=Dna5String>
    class SeqAnBloomFilter
    {
    public:

        typedef Shape<Dna5, UngappedShape<KMER_SIZE> > TShape;

        bool save(const char *fileName) const
        {
            std::ofstream myFile(fileName, std::ios::out | std::ios::binary);

            std::cerr << "Storing filter. Filter is " << m_sizeInBytes << " bytes." << std::endl;
            assert(myFile);

            //write out each block
            myFile.write(reinterpret_cast<char*>(_filterFile), m_sizeInBytes);

            myFile.close();
            assert(myFile);
            return true;
        }

        SeqAnBloomFilter()
        {
            _initPreCalcValues();
            initSize();
            memset(_filterFile, 0, m_sizeInBytes);
        }

        SeqAnBloomFilter(const char *fileName)
        {
            _initPreCalcValues();
            initSize();
            FILE *file = fopen(fileName, "rb");
            if (file == NULL)
            {
                std::cerr << "file \"" << fileName << "\" could not be read." << std::endl;
                exit(1);
            }

            long int lCurPos = ftell(file);
            fseek(file, 0, 2);
            size_t fileSize = ftell(file);
            fseek(file, lCurPos, 0);
            if (fileSize != m_sizeInBytes)
            {
                std::cerr << "Error: " << fileName
                << " does not match size given by its header. Size: "
                << fileSize << " vs " << m_sizeInBytes << " bytes." << std::endl;
                exit(1);
            }

            size_t countRead = std::fread(_filterFile, fileSize, 1, file);
            if (countRead != 1 && fclose(file) != 0)
            {
                std::cerr << "file \"" << fileName << "\" could not be read." << std::endl;
                exit(1);
            }
        }

        void addKmers(TString const & text)
        {
            _addKmers(text);
        }
        void addKmers(TString const & text, uint8_t const & binNo)
        {
            _addKmers(text, binNo);
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
                if(containsKmer(kmerHash))
                    --left;
                ++it;
            }
            return (left == 0);
        }

        std::vector<bool> whichBins(TString const & text, uint8_t const & threshold) const
        {
            TShape kmerShape;
            hashInit(kmerShape, begin(text));

            std::vector<bool> selected(BINS_SIZE);

            uint8_t possible = length(text) - length(kmerShape) + 1;
            std::vector<std::vector<uint8_t>> profileMatrix(possible, std::vector<uint8_t>(m_binSizeInChars, 255));

            auto it = begin(text);
            std::vector<uint64_t> kmerHashes(possible, 0);
            for (uint8_t i = 0; i < possible; ++i)
            {
                kmerHashes[i] = hashNext(kmerShape, it);
                ++it;
            }

            _whichBinsImpl(profileMatrix, kmerHashes);
            for (uint8_t binNo = 0; binNo < BINS_SIZE; ++binNo)
            {
                uint8_t count = 0;
                for (uint8_t i = 0; i < possible; ++i)
                {
                    uint8_t bit = bitMask[binNo % bitsPerChar];
                    if ((profileMatrix[i][binNo/bitsPerChar] & bit) == bit)
                    {
                        ++count;
                        if(count >= threshold)
                        {
                            selected[binNo] = true;
                            break;
                        }
                    }
                }
            }
            return selected;
        }

        inline void _whichBinsImpl(std::vector<std::vector<uint8_t> >  & profileMatrix,
                                   std::vector<uint64_t> const & kmerHashes) const
        {
            for(uint8_t k = 0; k < kmerHashes.size() ; k++)
            {
                uint64_t tmp = 0;
                for(uint8_t i = 0; i < N_HASH ; i++)
                {
                    tmp = kmerHashes[k] * (_preCalcValues[i]);
                    tmp ^= tmp >> _shiftValue;
                    uint64_t normalizedValue = (tmp % m_sizeInHashes) * BINS_SIZE;
                    for(uint8_t j = 0; j < m_binSizeInChars ; j++)
                    {
                        profileMatrix[k][j] &= _filterFile[normalizedValue / bitsPerChar + j];
                    }
                }
            }
        }

        inline std::vector<uint8_t> _whichBinsImpl(uint64_t & kmerHash) const
        {
            uint64_t tmp = 0;
            std::vector<uint8_t> mt;
            mt.resize(m_binSizeInChars, 255);
            for(uint8_t i = 0; i < N_HASH ; i++)
            {
                tmp = kmerHash * (_preCalcValues[i]);
                tmp ^= tmp >> _shiftValue;
                uint64_t normalizedValue = (tmp % m_sizeInHashes) * BINS_SIZE;
                for(uint8_t j = 0; j < m_binSizeInChars ; j++)
                {
                    mt[j] &= _filterFile[normalizedValue / bitsPerChar + j];
                }
            }
            return mt;
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

        void initSize()
        {
            if (SIZE % 8 != 0)
            {
                std::cerr << "ERROR: Filter Size \"" << SIZE << "\" is not a multiple of 8." << std::endl;
                exit(1);
            }
            m_sizeInBytes = SIZE / bitsPerChar;
            m_sizeInHashes = SIZE / BINS_SIZE;
            m_binSizeInChars = BINS_SIZE/bitsPerChar;
            _filterFile = new unsigned char[m_sizeInBytes];
        }


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

        void insertKmer(uint64_t & kmerHash)
        {
            uint64_t tmp = 0;
            for(uint8_t i = 0; i < N_HASH ; i++)
            {
                tmp = kmerHash * (i ^ KMER_SIZE * _seedValue);
                tmp ^= tmp >> _shiftValue;
                uint64_t normalizedValue = tmp % SIZE;
                __sync_or_and_fetch(&_filterFile[normalizedValue / bitsPerChar],
                                    bitMask[normalizedValue % bitsPerChar]);
            }
        }

        void insertKmer(uint64_t & kmerHash, uint8_t const & binNo)
        {
            uint64_t tmp = 0;
            for(uint8_t i = 0; i < N_HASH ; i++)
            {
                tmp = kmerHash * (_preCalcValues[i]);
                tmp ^= tmp >> _shiftValue;
                uint64_t normalizedValue = (tmp % m_sizeInHashes) * BINS_SIZE + binNo;
                __sync_or_and_fetch(&_filterFile[normalizedValue / bitsPerChar],
                                    bitMask[normalizedValue % bitsPerChar]);
            }
        }

        inline bool containsKmer(uint64_t & kmerHash) const
        {
            uint64_t tmp = 0;
            for(uint8_t i = 0; i < N_HASH ; i++)
            {
                tmp = kmerHash * (_preCalcValues[i]);
                tmp ^= tmp >> _shiftValue;
                size_t normalizedValue = tmp % SIZE;
                unsigned char bit = bitMask[normalizedValue % bitsPerChar];
                if ((_filterFile[normalizedValue / bitsPerChar] & bit) != bit)
                    return false;
            }
            return true;
        }

        inline void whichBinContainsKmer(char *setBins, uint64_t & kmerHash) const
        {
            uint64_t tmp = 0;
            for(uint8_t i = 0; i < N_HASH ; i++)
            {

                tmp = kmerHash * (_preCalcValues[i]);
                tmp ^= tmp >> _shiftValue;
                size_t normalizedValue = (tmp % m_sizeInHashes) * BINS_SIZE;
                setBins = _filterFile[normalizedValue / bitsPerChar] & setBins;
            }
        }

        bool containsKmer(std::vector<uint64_t> & hashValues)
        {
            for(uint8_t i = 0; i < N_HASH ; i++)
            {
                size_t normalizedValue = hashValues[i] % SIZE;
                unsigned char bit = bitMask[normalizedValue % bitsPerChar];
                if ((_filterFile[normalizedValue / bitsPerChar] & bit) != bit)
                    return false;
            }
            return true;
        }

        void _addKmers(TString const & text)
        {
            TShape kmerShape;
            hashInit(kmerShape, begin(text));

            for (uint32_t i = 0; i < length(text) - length(kmerShape) + 1; ++i)
            {
                uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);
                insertKmer(kmerHash);
            }
        }
        void _addKmers(TString const & text, uint8_t const & binNo)
        {
            TShape kmerShape;
            hashInit(kmerShape, begin(text));

            for (uint32_t i = 0; i < length(text) - length(kmerShape) + 1; ++i)
            {
                uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);
                insertKmer(kmerHash, binNo);
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

        size_t                  m_sizeInBytes;
        size_t                  m_sizeInHashes;
        size_t                  m_binSizeInChars;
        uint8_t*                _filterFile;
        std::vector<uint64_t>   _preCalcValues = {};
        uint64_t const          _shiftValue = 27;
        uint64_t const          _seedValue = 0x90b45d39fb6da1fa;
    };
}

//namespace seqan
//{
//    template<uint8_t KMER_SIZE, uint8_t N_HASH, uint32_t SIZE, typename TString=Dna5String>
//    class SeqAnBloomFilter
//    {
//    public:
//
//        typedef Shape<Dna5, UngappedShape<KMER_SIZE> > TShape;
//
//        bool open(const char *fileName)
//        {
//            return seqan::open(_filterFile, fileName, OPEN_RDONLY);
//        }
//
//        SeqAnBloomFilter(){
//            _initPreCalcValues();
//        }
//
//        SeqAnBloomFilter(const char *fileName)
//        {
//            openNewFile(fileName);
//        }
//
//        SeqAnBloomFilter(const char *fileName,
//                         TString & text)
//        {
//            openNewFile(fileName);
//            if(length(text) >= KMER_SIZE)
//                addKmers(text);
//        }
//
//        void addKmers(TString const & text)
//        {
//            _addKmers(text);
//        }
//
//        bool containsNKmers(TString const & text, uint8_t const & threshold) const
//        {
//            TShape kmerShape;
//            hashInit(kmerShape, begin(text));
//
//            uint8_t left = threshold;
//            uint8_t possible = length(text) - length(kmerShape) + 1;
//            auto it = begin(text);
//            uint64_t kmerHash = 0;
//            for (uint8_t i = 0; left > 0 && (possible - i) > left ; ++i)
//            {
//                kmerHash = hashNext(kmerShape, it);
////                getNHashValues(hashValues, kmerHash);
//                if(containsKmer(kmerHash))
//                    --left;
//                ++it;
//            }
//            return (left == 0);
//        }
//
//
//        bool containsNKmers(TString const & fwd, TString const & rev, uint32_t const & threshold) const
//        {
//            if(containsNKmers(fwd, threshold))
//                return true;
//            else if(containsNKmers(rev, threshold))
//                return true;
//            return false;
//        }
//
//        uint32_t kmerCount(TString const & text)
//        {
//            TShape kmerShape;
//            hashInit(kmerShape, begin(text));
//            uint32_t found = 0;
//            for (unsigned i = 0; i < length(text) - length(kmerShape) + 1 ; ++i)
//            {
//                uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);
//                std::vector<uint64_t> hashValues(N_HASH);
//                getNHashValues(hashValues, kmerHash);
//                if(containsKmer(kmerHash))
//                    ++found;
//            }
//            return found ;
//        }
//
//
//    private:
//        void getNHashValues(std::vector<uint64_t> & hashValues, uint64_t & kmerHash)
//        {
//            uint64_t tmp = 0;
//            for(uint8_t i = 0; i < N_HASH ; i++)
//            {
//                tmp = kmerHash * (i ^ KMER_SIZE * _seedValue);
//                tmp ^= tmp >> _shiftValue;
//                hashValues[i] = tmp;
//            }
//        }
//
//        void insertKmer(std::vector<uint64_t> & hashValues)
//        {
//            for(uint8_t i = 0; i < N_HASH ; i++)
//                assignValue(_filterFile, (hashValues[i] % SIZE), true);
//        }
//
//        void insertKmer(uint64_t & kmerHash)
//        {
//            uint64_t tmp = 0;
//            for(uint8_t i = 0; i < N_HASH ; i++)
//            {
//                tmp = kmerHash * (i ^ KMER_SIZE * _seedValue);
//                tmp ^= tmp >> _shiftValue;
//                assignValue(_filterFile, (tmp % SIZE), true);
//            }
//        }
//
//        inline bool containsKmer(uint64_t & kmerHash) const
//        {
//            uint64_t tmp = 0;
//            for(uint8_t i = 0; i < N_HASH ; i++)
//            {
//                tmp = kmerHash * (_preCalcValues[i]);
//                tmp ^= tmp >> _shiftValue;
//                if(!_filterFile[(tmp % SIZE)])
//                    return false;
//            }
//            return true;
//        }
//
//        bool containsKmer(std::vector<uint64_t> & hashValues)
//        {
//            for(uint8_t i = 0; i < N_HASH ; i++)
//                if(!_filterFile[(hashValues[i] % SIZE)])
//                    return false;
//            return true;
//        }
//
//        void _addKmers(TString const & text)
//        {
//            TShape kmerShape;
//            hashInit(kmerShape, begin(text));
//
//            for (uint32_t i = 0; i < length(text) - length(kmerShape) + 1; ++i)
//            {
//                uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);
//                std::vector<uint64_t> hashValues(N_HASH);
//                getNHashValues(hashValues, kmerHash);
//                insertKmer(hashValues);
//            }
//        }
//
//        inline void _initPreCalcValues()
//        {
//            for(uint8_t i = 0; i < N_HASH ; i++)
//            {
//                _preCalcValues.push_back(i ^ KMER_SIZE * _seedValue);
//            }
//        }
//        void openNewFile(const char *fileName)
//        {
//            seqan::open(_filterFile, fileName, OPEN_RDWR | OPEN_CREATE);
//            reserve(_filterFile, SIZE, Exact());
//            resize(_filterFile, SIZE, false);
//        }
//        std::vector<uint64_t> _preCalcValues = {};
//        uint64_t const _shiftValue = 27;
//        uint64_t const _seedValue = 0x90b45d39fb6da1fa;
//        String<bool, Packed<MMap<> > > _filterFile;
//    };
//}

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
