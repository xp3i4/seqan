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

std::mutex mtx;

using namespace seqan;
class Semaphore
{
    std::mutex m;
    std::condition_variable cv;
    int count;

public:
    Semaphore(int n) : count{n} {}
    void notify()
    {
        std::unique_lock<std::mutex> l(m);
        ++count;
        cv.notify_one();
    }
    void wait()
    {
        std::unique_lock<std::mutex> l(m);
        cv.wait(l, [this]{ return count!=0; });
        --count;
    }
};

class Critical_section
{
    Semaphore &s;
public:
    Critical_section(Semaphore &ss) : s{ss} { s.wait(); }
    ~Critical_section() { s.notify(); }
};


namespace seqan
{
    static const uint8_t bitsPerChar = 0x08;
    static const uint8_t bitMask[0x08] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

    template<uint8_t BINS_SIZE, uint8_t KMER_SIZE, uint8_t N_HASH, uint64_t SIZE, typename TString=Dna5String>
    class SeqAnBloomFilter
    {
    public:

        typedef Shape<Dna, UngappedShape<KMER_SIZE> > TShape;

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

            uint64_t lCurPos = ftell(file);
            fseek(file, 0, 2);
            uint64_t fileSize = ftell(file);
            fseek(file, lCurPos, 0);
            if (fileSize != m_sizeInBytes)
            {
                std::cerr << "Error: " << fileName
                << " does not match size given by its header. Size: "
                << fileSize << " vs " << m_sizeInBytes << " bytes." << std::endl;
                exit(1);
            }

            uint64_t countRead = std::fread(_filterFile, fileSize, 1, file);
            if (countRead != 1 && fclose(file) != 0)
            {
                std::cerr << "file \"" << fileName << "\" could not be read." << std::endl;
                exit(1);
            }
        }

        void addKmers(TString const & text, uint8_t const & binNo)
        {
            _addKmers(text, binNo);
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

            _filProfileMatrix(profileMatrix, kmerHashes);
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

        inline void _filProfileMatrix(std::vector<std::vector<uint8_t> >  & profileMatrix,
                                      std::vector<uint64_t> const & kmerHashes) const
        {

            Semaphore thread_limiter(8);
            std::vector<std::future<void>> tasks;

            for(uint8_t k = 0; k < kmerHashes.size() ; k++)
            {
                tasks.emplace_back(std::async([=, &thread_limiter, &profileMatrix] {
                    uint64_t tmp = 0;
                    std::vector<uint8_t> tmpProfile(m_binSizeInChars, 255);
                    for(uint8_t i = 0; i < N_HASH ; i++)
                    {
                        tmp = kmerHashes[k] * (_preCalcValues[i]);
                        tmp ^= tmp >> _shiftValue;
                        uint64_t normalizedValue = (tmp % m_sizeInHashes) * BINS_SIZE;
                        for(uint8_t j = 0; j < m_binSizeInChars ; j++)
                        {
                            tmpProfile[j] &= _filterFile[normalizedValue / bitsPerChar + j];
                        }
                    }
                    profileMatrix[k] = std::move(tmpProfile);
                }));
            }
            for (auto &&task : tasks)
            {
                task.get();
            }
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
            _filterFile = new uint8_t[m_sizeInBytes];
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

        void _addKmers(TString const & text, uint8_t const & binNo)
        {
            TShape kmerShape;
            hashInit(kmerShape, begin(text));

            uint8_t possible = length(text) - length(kmerShape) + 1;
            auto it = begin(text);
            for (uint8_t i = 0; i < possible; ++i)
            {
                uint64_t kmerHash = hashNext(kmerShape, it);
                insertKmer(kmerHash, binNo);
                ++it;
            }
        }


        inline void _initPreCalcValues()
        {
            for(uint8_t i = 0; i < N_HASH ; i++)
            {
                _preCalcValues.push_back(i ^ KMER_SIZE * _seedValue);
            }
        }


        uint64_t                  m_sizeInBytes;
        uint64_t                  m_sizeInHashes;
        uint64_t                  m_binSizeInChars;
        uint8_t*                _filterFile;
        std::vector<uint64_t>   _preCalcValues = {};
        uint64_t const          _shiftValue = 27;
        uint64_t const          _seedValue = 0x90b45d39fb6da1fa;
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


template<typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v) {
    char comma[3] = {'\0', ' ', '\0'};
    for (const auto& e : v) {
        s << comma << e;
        comma[0] = ',';
    }
    return s;
}

#endif  // #ifndef APP_YARA_MISC_OPTIONS_DIS_H_
