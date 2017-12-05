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
#include <sdsl/bit_vectors.hpp>

namespace seqan
{
    static const uint8_t uInt64Width = 0x40;
    static const uint64_t bitMask[0x40] = {
        0x0000000000000001, 0x0000000000000002, 0x0000000000000004, 0x0000000000000008,
        0x0000000000000010, 0x0000000000000020, 0x0000000000000040, 0x0000000000000080,
        0x0000000000000100, 0x0000000000000200, 0x0000000000000400, 0x0000000000000800,
        0x0000000000001000, 0x0000000000002000, 0x0000000000004000, 0x0000000000008000,
        0x0000000000010000, 0x0000000000020000, 0x0000000000040000, 0x0000000000080000,
        0x0000000000100000, 0x0000000000200000, 0x0000000000400000, 0x0000000000800000,
        0x0000000001000000, 0x0000000002000000, 0x0000000004000000, 0x0000000008000000,
        0x0000000010000000, 0x0000000020000000, 0x0000000040000000, 0x0000000080000000,
        0x0000000100000000, 0x0000000200000000, 0x0000000400000000, 0x0000000800000000,
        0x0000001000000000, 0x0000002000000000, 0x0000004000000000, 0x0000008000000000,
        0x0000010000000000, 0x0000020000000000, 0x0000040000000000, 0x0000080000000000,
        0x0000100000000000, 0x0000200000000000, 0x0000400000000000, 0x0000800000000000,
        0x0001000000000000, 0x0002000000000000, 0x0004000000000000, 0x0008000000000000,
        0x0010000000000000, 0x0020000000000000, 0x0040000000000000, 0x0080000000000000,
        0x0100000000000000, 0x0200000000000000, 0x0400000000000000, 0x0800000000000000,
        0x1000000000000000, 0x2000000000000000, 0x4000000000000000, 0x8000000000000000};

    template<typename TString = Dna5String>
    class SeqAnBloomFilter
    {
    public:

        typedef Shape<Dna, SimpleShape> TShape;

        SeqAnBloomFilter(uint32_t n_bins, uint8_t n_hash_func, uint8_t kmer_size, uint64_t vec_size):
                        _noOfBins(n_bins),
                        _noOfHashFunc(n_hash_func),
                        _kmerSize(kmer_size),
                        _noOfBits(vec_size),
                        _filterVector(sdsl::bit_vector(vec_size, 0));

        {
            _init();
        }

//        SeqAnBloomFilter(const char *fileName, uint32_t n_bins, uint8_t n_hash_func, uint8_t kmer_size, uint64_t vec_size):
//                        _noOfBins(n_bins),
//                        _noOfHashFunc(n_hash_func),
//                        _kmerSize(kmer_size),
//                        _noOfBits(vec_size)
//        {
//            _init();
//            std::ifstream inStream(fileName, std::ios::binary);
//            inStream.read(reinterpret_cast<char*>(&_filterVector[0]), _filterVector.size()*sizeof(uint64_t));
//        }

        SeqAnBloomFilter(const char *fileName, uint32_t n_bins, uint8_t n_hash_func, uint8_t kmer_size, uint64_t vec_size):
                        _noOfBins(n_bins),
                        _noOfHashFunc(n_hash_func),
                        _kmerSize(kmer_size),
                        _noOfBits(vec_size)
        {
            _init();
            if (!sdsl::load_from_file(_filterVector, fileName))
            {
                std::cerr << "File \"" << fileName << "\" could not be read." << std::endl;
                exit(1);
            }
            if (vec_size != _filterVector.bit_size())
            {
                std::cerr << "Size mismatch: \n\t Loaded file \t" << _filterVector.bit_size()
                << " bits\n\t expected\t" << vec_size << std::endl;
                exit(1);
            }
        }

        void addKmers(TString const & text, uint32_t const & binNo)
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
//                    if (_containsKmer(kmerHash, binNo))
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
                //if the number of bu
                for (uint8_t batchNo = 0; batchNo < _binIntWidth; ++batchNo)
                {
                    uint32_t binNo = batchNo * uInt64Width;
                    std::bitset<64> bitSet = _containsKmerBatch(kmerHash, batchNo);;
                    if(bitSet.none()) continue;
                    for(uint8_t offset=0; binNo < _noOfBins && offset < uInt64Width; ++offset,++binNo)
                    {
                        if (!selected[binNo] && bitSet.test(offset))
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

        //        bool save(const char *fileName) const
        //        {
        //            std::ofstream outStream(fileName, std::ios::out | std::ios::binary);
        //
        //            std::cerr << "Storing filter. Filter is " << _noOfBits << " bytes." << std::endl;
        //            assert(outStream);
        //
        //            outStream.write(reinterpret_cast<const char*>(&_filterVector[0]), _filterVector.size()*sizeof(uint64_t));
        //
        //            outStream.close();
        //            assert(outStream);
        //            return true;
        //        }

        //save case sdsl
        bool save(const char *fileName)
        {
            std::cerr << "Storing filter. Filter is " << size_in_mega_bytes(_filterVector) << " MB." << std::endl;
            return sdsl::store_to_file(_filterVector, fileName);
        }

    private:

        void _init()
        {
            _initPreCalcValues();
            _binIntWidth = std::ceil((float)_noOfBins / uInt64Width);
            _noOfHashPos = _noOfBits / (uInt64Width * _binIntWidth);
//            _filterVector.resize(_noOfBits / uInt64Width);
        }


        bool _isBitSet(uint8_t num, uint8_t bit) const
        {
            return 1 == ( (num >> bit) & 1);
        }

//        std::bitset<64> _containsKmerBatch(uint64_t & kmerHash, uint8_t const & batch) const
//        {
//            uint64_t tmp = kmerHash * (_preCalcValues[0]);
//            tmp ^= tmp >> _shiftValue;
//            uint64_t vectIndex = (tmp % _noOfHashPos) * _binIntWidth + batch;
//
//            std::bitset<64> res(_filterVector[vectIndex]);
//
//            for(uint8_t i = 1; i < _noOfHashFunc ; i++)
//            {
//                tmp = kmerHash * (_preCalcValues[i]);
//                tmp ^= tmp >> _shiftValue;
//                vectIndex = (tmp % _noOfHashPos) * _binIntWidth + batch;
//                res &= _filterVector[vectIndex];
//            }
//            return res;
//        }

        //sdsl case
        std::bitset<64> _containsKmerBatch(uint64_t & kmerHash, uint8_t const & batch) const
        {
            uint64_t tmp = kmerHash * (_preCalcValues[0]);
            tmp ^= tmp >> _shiftValue;
            uint64_t vectIndex = (tmp % _noOfHashPos) * _binIntWidth + batch;

            std::bitset<64> res(_filterVector.get_int(vectIndex));

            for(uint8_t i = 1; i < _noOfHashFunc ; i++)
            {
                tmp = kmerHash * (_preCalcValues[i]);
                tmp ^= tmp >> _shiftValue;
                vectIndex = (tmp % _noOfHashPos) * _binIntWidth + batch;
                res &= _filterVector.get_int(vectIndex);
            }
            return res;
        }

        bool _containsKmer(uint64_t & kmerHash, uint32_t const & binNo) const
        {
            uint64_t tmp = 0;
            uint8_t binOfset = binNo%uInt64Width;
            for(uint8_t i = 0; i < _noOfHashFunc ; i++)
            {
                tmp = kmerHash * (_preCalcValues[i]);
                tmp ^= tmp >> _shiftValue;
                uint64_t vectIndex = (tmp % _noOfHashPos) * _binIntWidth + binNo/uInt64Width;
                std::bitset<64> res(_filterVector[vectIndex]);
                if (!res.test(binOfset))
                    return false;
            }
            return true;
        }

//        void _insertKmer(uint64_t & kmerHash, uint32_t const & binNo)
//        {
//            uint64_t tmp = 0;
//            uint8_t binOfset = binNo%uInt64Width;
//            for(uint8_t i = 0; i < _noOfHashFunc ; i++)
//            {
//                tmp = kmerHash * (_preCalcValues[i]);
//                tmp ^= tmp >> _shiftValue;
//                uint64_t vectIndex = (tmp % _noOfHashPos) * _binIntWidth + binNo/uInt64Width;
//                __sync_or_and_fetch(&_filterVector[vectIndex], bitMask[binOfset]);
//            }
//        }

        //sdsl vector
        void _insertKmer(uint64_t & kmerHash, uint32_t const & binNo)
        {
            uint64_t tmp = kmerHash * (_preCalcValues[0]);
            tmp ^= tmp >> _shiftValue;

            uint64_t vectIndex = (tmp % _noOfHashPos) * _binIntWidth + binNo;
            _filterVector[vectIndex] = true;

            for(uint8_t i = 1; i < _noOfHashFunc ; i++)
            {
                tmp = kmerHash * (_preCalcValues[i]);
                tmp ^= tmp >> _shiftValue;
                vectIndex = (tmp % _noOfHashPos) * _binIntWidth + binNo;
                _filterVector[vectIndex] = true;
            }
        }

        void _addKmers(TString const & text, uint32_t const & binNo)
        {
            TShape kmerShape;
            resize(kmerShape, _kmerSize);
            hashInit(kmerShape, begin(text));

            for (uint32_t i = 0; i < length(text) - length(kmerShape) + 1; ++i)
            {
                uint64_t kmerHash = hashNext(kmerShape, begin(text) + i);
                _insertKmer(kmerHash, binNo);
            }
        }


        inline void _initPreCalcValues()
        {
            for(uint8_t i = 0; i < _noOfHashFunc ; i++)
            {
                _preCalcValues.push_back(i ^ _kmerSize * _seedValue);
            }
        }


        uint32_t                 _noOfBins;
        uint8_t                  _noOfHashFunc;
        uint8_t                  _kmerSize;
        uint8_t                  _binIntWidth;

        //sizes in diferent units
        size_t                  _noOfBits;
        size_t                  _noOfHashPos;
        sdsl::bit_vector        _filterVector;
//        std::vector<uint64_t>   _filterVector;
        std::vector<uint64_t>   _preCalcValues = {};
        uint64_t const          _shiftValue = 27;
        uint64_t const          _seedValue = 0x90b45d39fb6da1fa;
    };
}
