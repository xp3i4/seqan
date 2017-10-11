// ==========================================================================
//                                 classify_reads.h
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

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class SuperOptions
// --------------------------------------------------------------------------
using namespace seqan;

struct ClassifyReadsOption
{
    bool                singleEnd;
    Pair<CharString>    reads_file;
    CharString          kmer_index_file;

    uint32_t            max_errors = 3;
    uint32_t            error_diff;

    void reload()
    {
        error_diff =  K_MER_LENGTH - 1 + (max_errors * K_MER_LENGTH);
    }

};

// ==========================================================================
// Functions
// ==========================================================================

template <typename TShapeValue, typename TMapValue>
inline std::bitset<NUM_OF_BINS> which_bins(DnaString const & seq,
                                           ClassifyReadsOption const & options,
                                           std::unordered_map <TShapeValue, TMapValue> const & kmer_occurance_table)
{


    // q-gram based threshold
    uint32_t threshold = static_cast<uint32_t>(length(seq)) - options.error_diff;

    // open file
    std::bitset<NUM_OF_BINS> candidate_bins(0);
    std::vector<uint32_t>  sum(NUM_OF_BINS, 0);

    Shape<Dna, UngappedShape<K_MER_LENGTH> > kmer_shape;
    hashInit(kmer_shape, begin(seq));


    for (unsigned i = 0; i < length(seq) - K_MER_LENGTH + 1; ++i)
    {
        TShapeValue key = hashNext(kmer_shape, begin(seq) + i);

        auto iter = kmer_occurance_table.find(key);
        if(iter != kmer_occurance_table.end())
        {
            for (uint32_t bin_no=0; bin_no < NUM_OF_BINS; ++bin_no)
            {
                if (iter->second[bin_no])
                    ++sum[bin_no];
            }
        }
    }
    for (uint32_t bin_no=0; bin_no < NUM_OF_BINS; ++bin_no)
    {
        if (sum[bin_no] >= threshold)
            candidate_bins.set(bin_no);
    }

    return candidate_bins;
}

template <typename TKey, typename TReads>
inline void get_candidate_bins (TReads const & reads,
                                ClassifyReadsOption const & options,
                                std::vector<std::bitset<NUM_OF_BINS> > & candidate_bins,
                                std::unordered_map <TKey, std::bitset<NUM_OF_BINS> > const  & kmer_occurance_table)
{
    // parallel
    for (uint32_t seq_no = 0; seq_no < length(reads.seqs); ++seq_no)
    {
        candidate_bins[seq_no] = which_bins(reads.seqs.strings[seq_no], options, kmer_occurance_table);
    }
}

template <typename TShapeValue, typename TMapValue>
inline void print_read_distribution(ClassifyReadsOption const & options,
                                    std::unordered_map <TShapeValue, TMapValue> const & kmer_occurance_table,
                                    CharString const & readPath)
{

    // ----------------------------------------------------------------------------
    // temporary print which bins are candidate for individual reads
    // ----------------------------------------------------------------------------

    std::vector<uint64_t> read_count(NUM_OF_BINS, 0);
    std::vector<uint64_t> bin_count_freq(NUM_OF_BINS+1, 0);

    SeqFileIn read_file;
    open(read_file , toCString(readPath));

    String<Dna5> seq;
    CharString id;

    while(!atEnd(read_file))
    {
        uint32_t bin_count = 0;
        readRecord(id, seq, read_file);
        std::bitset<NUM_OF_BINS> bin_bits = which_bins(seq, options, kmer_occurance_table);
        for (uint32_t bin_no = 0; bin_no< NUM_OF_BINS; ++bin_no)
        {
            if (bin_bits.test(bin_no))
            {
                ++bin_count;
                ++read_count[bin_no];
            }
        }
        ++bin_count_freq[bin_count];
    }
    close(read_file);

    uint32_t first_dot_pos = 0;
    for (; first_dot_pos < length(readPath); ++first_dot_pos)
    {
        if( readPath[first_dot_pos] == '.')
            break;
    }

    CharString dist_file_base = prefix(readPath, first_dot_pos);
    append(dist_file_base, "_k");
    append(dist_file_base, std::to_string(K_MER_LENGTH));
    append(dist_file_base, "_e");
    append(dist_file_base, std::to_string(options.max_errors));

    CharString freq_file_path = dist_file_base;
    CharString count_file_path = dist_file_base;
    append(freq_file_path, "_num_bins_freq.tsv");
    append(count_file_path, "_num_reads.tsv");

    std::ofstream freq_output_file;
    std::ofstream count_output_file;
    freq_output_file.open (toCString(freq_file_path));
    count_output_file.open (toCString(count_file_path));

    for(size_t bin_no = 0; bin_no < NUM_OF_BINS; ++bin_no)
    {
        count_output_file <<"bin_" << bin_no  << "\t" << read_count[bin_no] << "\n";
        freq_output_file << bin_no  << "\t" << bin_count_freq[bin_no] << "\n";
    }
    //finish the last bin
    freq_output_file << NUM_OF_BINS  << "\t" << bin_count_freq[NUM_OF_BINS] << "\n";
    count_output_file.close();
    freq_output_file.close();
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
}

template <typename TShapeValue, typename TMapValue>
inline void annotate_reads(ClassifyReadsOption const & options,
                           std::unordered_map <TShapeValue, TMapValue> const & kmer_occurance_table,
                           CharString const & readPath)
{

    // ----------------------------------------------------------------------------
    // temporary print which bins are candidate for individual reads
    // ----------------------------------------------------------------------------

    SeqFileIn read_file;
    SeqFileOut annotated_read_file;
    open(read_file , toCString(readPath));

    uint32_t first_dot_pos = 0;
    for (; first_dot_pos < length(readPath); ++first_dot_pos)
    {
        if(readPath[first_dot_pos] == '.')
            break;
    }
    CharString annotated_read_file_path = prefix(readPath, first_dot_pos);
    append(annotated_read_file_path, "_k");
    append(annotated_read_file_path, std::to_string(K_MER_LENGTH));
    append(annotated_read_file_path, "_e");
    append(annotated_read_file_path, std::to_string(options.max_errors));
    append(annotated_read_file_path, suffix(readPath, first_dot_pos));

    open(annotated_read_file , toCString(annotated_read_file_path));

    String<Dna5> seq;
    CharString id, qual;

    while(!atEnd(read_file))
    {
        readRecord(id, seq, qual, read_file);
        append(id, "|");
        std::bitset<NUM_OF_BINS> bin_bits = which_bins(seq, options, kmer_occurance_table);
        for (uint32_t bin_no = 0; bin_no < NUM_OF_BINS; ++bin_no)
        {
            if (bin_bits.test(bin_no))
            {
                append(id, std::to_string(bin_no));
                append(id, "|");
            }
        }
        writeRecord(annotated_read_file, id, seq, qual);
    }
    close(read_file);
    close(annotated_read_file);
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
}



template <typename TKey>
inline void purify_kmer_index(std::unordered_map <TKey, std::bitset<NUM_OF_BINS> > &kmer_occurance_table)
{
    for(auto it = kmer_occurance_table.begin(); it != kmer_occurance_table.end(); )
        if(!(it->second).any())
            it = kmer_occurance_table.erase(it);
        else
            ++it;
}
