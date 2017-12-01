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

namespace seqan
{
    static const uint8_t bitsPerChar = 0x08;
    static const uint8_t bitMask[0x08] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

    template<typename TString = Dna5String>
    class SeqAnBloomFilter
    {
    public:

        typedef Shape<Dna, SimpleShape> TShape;

        bool save(const char *fileName) const
        {
            std::ofstream outStream(fileName, std::ios::out | std::ios::binary);

            std::cerr << "Storing filter. Filter is " << _noOfBytes << " bytes." << std::endl;
            assert(outStream);

            outStream.write(reinterpret_cast<const char*>(&_filterVector[0]), _filterVector.size()*sizeof(uint8_t));

            outStream.close();
            assert(outStream);
            return true;
        }

        SeqAnBloomFilter(uint8_t n_bins, uint8_t n_hash_func, uint8_t kmer_size, uint64_t vec_size):
                        _noOfBins(n_bins),
                        _noOfHashFunc(n_hash_func),
                        _kmerSize(kmer_size),
                        _noOfBytes(vec_size)
        {
            _init();
        }

        SeqAnBloomFilter(const char *fileName, uint8_t n_bins, uint8_t n_hash_func, uint8_t kmer_size, uint64_t vec_size):
                        _noOfBins(n_bins),
                        _noOfHashFunc(n_hash_func),
                        _kmerSize(kmer_size),
                        _noOfBytes(vec_size)
        {
            _init();
            std::ifstream inStream(fileName, std::ios::binary);
            inStream.read(reinterpret_cast<char*>(&_filterVector[0]), vec_size*sizeof(uint8_t));
        }

        void addKmers(TString const & text, uint8_t const & binNo)
        {
            _addKmers(text, binNo);
        }

        void whichBins(std::vector<bool> & selected, TString const & text, uint8_t const & threshold) const
        {
            TShape kmerShape;
            resize(kmerShape, _kmerSize);
            hashInit(kmerShape, begin(text));

            std::vector<uint8_t> counts(_noOfBins, 0);

            uint8_t possible = length(text) - length(kmerShape) + 1;

            std::vector<uint64_t> kmerHashes(possible, 0);
            auto it = begin(text);
            for (uint8_t i = 0; i < possible; ++i)
            {
                kmerHashes[i] = hashNext(kmerShape, it);
                ++it;
            }


//            for (uint64_t kmerHash : kmerHashes)
//            {
//                for (uint8_t binNo = 0; binNo < _noOfBins; ++binNo)
//                {
//                    if (threshold - counts[binNo] > possible || selected[binNo])
//                        continue;
//
//                    if (containsKmer(kmerHash, binNo))
//                    {
//                        ++counts[binNo];
//                        if(counts[binNo] >= threshold)
//                            selected[binNo] = true;
//                    }
//                }
//                --possible;
//            }
            
            for (uint64_t kmerHash : kmerHashes)
            {
                for (uint8_t batchNo = 0; batchNo < _binWidthInBytes; ++batchNo)
                {
                    uint8_t batchRes = containsKmerBatch(kmerHash, batchNo);
                    if(batchRes == 0) continue;
                    for(uint8_t offset=0; offset<bitsPerChar; ++offset)
                    {
                        uint8_t binNo = batchNo * bitsPerChar + offset;
                        if (threshold - counts[binNo] > possible || selected[binNo]) continue;

                        if (IsBitSet(batchRes, offset))
                        {
                            ++counts[binNo];
                            if(counts[binNo] >= threshold)
                                selected[binNo] = true;
                        }
                    }
                }
                --possible;
            }

        }

        std::vector<bool> whichBins(TString const & text, uint8_t const & threshold) const
        {
            std::vector<bool> selected(_noOfBins, false);
            whichBins(selected, text, threshold);
            return selected;
        }

        std::vector<bool> whichBins(TString const & text_fwd, TString const & text_rev, uint8_t const & threshold) const
        {
            std::vector<bool> selected(_noOfBins, false);
            whichBins(selected, text_fwd, threshold);
            whichBins(selected, text_rev, threshold);
            return selected;
        }

    private:

        void _init()
        {
            _initPreCalcValues();
            _binWidthInBytes = _noOfBins/bitsPerChar;
            _noOfHashePos = _noOfBytes/_binWidthInBytes;
            _filterVector.resize(_noOfBytes);
        }

        void insertKmer(uint64_t & kmerHash, uint8_t const & binNo)
        {
            uint64_t tmp = 0;
            for(uint8_t i = 0; i < _noOfHashFunc ; i++)
            {
                tmp = kmerHash * (_preCalcValues[i]);
                tmp ^= tmp >> _shiftValue;
                uint64_t normalizedValue = (tmp % _noOfHashePos) * _noOfBins + binNo;
                __sync_or_and_fetch(&_filterVector[normalizedValue / bitsPerChar],
                                    bitMask[normalizedValue % bitsPerChar]);
            }
        }
        bool IsBitSet(uint8_t num, uint8_t bit) const
        {
            return 1 == ( (num >> bit) & 1);
        }

        uint8_t containsKmerBatch(uint64_t & kmerHash, uint8_t const & batch) const
        {
            uint8_t res = 255;
            uint64_t tmp = 0;
            uint64_t normalizedValue;
            for(uint8_t i = 0; i < _noOfHashFunc ; i++)
            {
                tmp = kmerHash * (_preCalcValues[i]);
                tmp ^= tmp >> _shiftValue;
                normalizedValue = (tmp % _noOfHashePos) * _noOfBins;
                res = res & _filterVector[normalizedValue / bitsPerChar + batch];
            }
            return res;
        }

        bool containsKmer(uint64_t & kmerHash, uint8_t const & binNo) const
        {
            uint64_t tmp = 0;
            for(uint8_t i = 0; i < _noOfHashFunc ; i++)
            {
                tmp = kmerHash * (_preCalcValues[i]);
                tmp ^= tmp >> _shiftValue;
                uint64_t normalizedValue = (tmp % _noOfHashePos) * _noOfBins + binNo;
                if (!IsBitSet(_filterVector[normalizedValue / bitsPerChar], (binNo % bitsPerChar)))
                    return false;
            }
            return true;
        }

        void _addKmers(TString const & text, uint8_t const & binNo)
        {
            TShape kmerShape;
            resize(kmerShape, _kmerSize);
            hashInit(kmerShape, begin(text));

            for (uint32_t i = 0; i < length(text) - length(kmerShape) + 1; ++i)
            {
                uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);
                insertKmer(kmerHash, binNo);
            }
        }


        inline void _initPreCalcValues()
        {
            for(uint8_t i = 0; i < _noOfHashFunc ; i++)
            {
                _preCalcValues.push_back(i ^ _kmerSize * _seedValue);
            }
        }


        uint8_t                  _noOfBins;
        uint8_t                  _noOfHashFunc;
        uint8_t                  _kmerSize;
        uint8_t                  _binWidthInBytes;

        //sizes in diferent units
        size_t                  _noOfBytes;
        size_t                  _noOfHashePos;
        std::vector<uint8_t>    _filterVector;
        std::vector<uint64_t>   _preCalcValues = {};
        uint64_t const          _shiftValue = 27;
        uint64_t const          _seedValue = 0x90b45d39fb6da1fa;
    };
}
