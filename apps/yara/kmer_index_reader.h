// ==========================================================================
//                                 kmer_index_reader.h
// ==========================================================================
// Copyright (c) 2017-2022, Temesgen H. Dadi, FU Berlin
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
//     * Neither the name of Temesgen H. Dadi or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL TEMESGEN H. DADI OR THE FU BERLIN BE LIABLE
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

template <typename TShapeValue, typename TMapValue, typename TOptions>
inline void read_occurance_table(TOptions const & options, std::unordered_map <TShapeValue, TMapValue> &kmer_occurance_table)
{
    std::ifstream istrm;
    if (!open(istrm,toCString(options.kmer_index_file), OPEN_RDONLY))
    {
        std::cerr << "File: " << options.kmer_index_file << std::endl;
        throw "failed to open!";
    }

    cereal::BinaryInputArchive iarchive(istrm);
    iarchive(kmer_occurance_table );

    istrm.close();
}


template <typename TShapeValue, typename TMapValue, typename TOptions>
inline void read_occurance_table_raw(TOptions const & options, std::unordered_map <TShapeValue, TMapValue> &kmer_occurance_table)
{
    // open file
    std::ifstream istrm(toCString(options.kmer_index_file), std::ios::binary);
    if (!istrm.is_open())
    {
        std::cerr << "File: " << options.kmer_index_file << std::endl;
        throw "failed to open!";
    }

    // check if lengths are equal
    std::string kmer, occurance;
    while (std::getline (istrm, kmer, '\t'))
    {
        std::getline (istrm, occurance, '\n');
        Shape<Dna, UngappedShape<K_MER_LENGTH> > kmer_shape;
        DnaString kmer2 = kmer;

        TShapeValue k = hash(kmer_shape, begin(kmer2));
        std::bitset<NUM_OF_BINS> occurance_bits(occurance);
        kmer_occurance_table.emplace_hint(kmer_occurance_table.end(), std::make_pair(k, occurance_bits));// [kmer_to_int(kmer)] = occurance_bits;
    }
    // close file
    istrm.close();
}
