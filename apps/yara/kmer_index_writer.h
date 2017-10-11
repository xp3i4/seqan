// ==========================================================================
//                                 kmer_index_writer.h
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

template <typename TShapeValue, typename TMapValue>
inline void write_kmer_index(AppOptions & options,
                              std::unordered_map <TShapeValue, TMapValue> &kmer_occurance_table)
{

    std::ofstream os(toCString(options.output_path), std::ios::binary);
    cereal::BinaryOutputArchive archive( os );

    archive( kmer_occurance_table );
}

// without filter old raw text
template <typename TShapeValue, typename TMapValue>
inline void write_kmer_index_raw(AppOptions & options,
                              std::unordered_map <TShapeValue, TMapValue> &kmer_occurance_table)
{

    //output file
    std::ofstream output_file;
    output_file.open (toCString(options.output_path));


    for(auto it = kmer_occurance_table.begin(); it != kmer_occurance_table.end(); )
    {
        DnaString needle;
        unhash(needle, it->first, K_MER_LENGTH);
        output_file << needle << "\t" << it->second  << std::endl;
        it = kmer_occurance_table.erase(it);
    }
    output_file.close();

}

template <typename TShapeValue, typename TMapValue>
inline void purify_kmer_index(std::unordered_map <TShapeValue, TMapValue > &kmer_occurance_table)
{
    for(auto it = kmer_occurance_table.begin(); it != kmer_occurance_table.end(); )
        if(!(it->second).any())
            it = kmer_occurance_table.erase(it);
        else
            ++it;
}
