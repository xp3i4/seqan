// ==========================================================================
//                                 classify-reads
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

#ifndef APP_YARA_DIS_MAPPER_H_
#define APP_YARA_DIS_MAPPER_H_

using namespace seqan;

typedef Mapper<Owner<ConcatDirect<> >, ReadMapperConfig<> > TDefaultMapper;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class SuperOptions
// --------------------------------------------------------------------------

struct SuperOptions : public Options
{
    CharString          superContigsIndicesFile;
    CharString          superOutputFile;
    uint32_t            NUM_OF_BINS = 64;

//    Pair<CharString>    superReadsFile;

//    uint32_t            max_errors = 3;
//    CharString          kmer_index_file;
//    uint32_t error_diff;
//
//    void reload()
//    {
//        error_diff =  K_MER_LENGTH - 1 + (max_errors * K_MER_LENGTH);
//    }
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function set_current_index_file()
// ----------------------------------------------------------------------------

inline void set_current_index_file(Options & yaraOptions, SuperOptions const & options, uint32_t const file_no)
{
    // Get the current file name
    yaraOptions.contigsIndexFile = options.superContigsIndicesFile;
//    if (file_no < 10)
//        append(yaraOptions.contigsIndexFile, "0");
    append(yaraOptions.contigsIndexFile, std::to_string(file_no));
//    append(yaraOptions.contigsIndexFile, "-genomes");
}

template <typename TOptions>
inline void set_current_index_file(TOptions & options, uint32_t const file_no)
{
    // Get the current file name
    options.contigsIndexFile = options.superContigsIndicesFile;
//    if (file_no < 10)
//        append(options.contigsIndexFile, "0");
    append(options.contigsIndexFile, std::to_string(file_no));
//    append(options.contigsIndexFile, "-genomes");
}

inline void set_output_file(Options & yaraOptions, SuperOptions const & options, uint32_t const file_no)
{
    uint32_t first_dot_pos = 0;
    for (; first_dot_pos < length(options.superOutputFile); ++first_dot_pos)
    {
        if(options.superOutputFile[first_dot_pos] == '.')
            break;
    }
    // Get the current file name
    yaraOptions.outputFile = prefix(options.superOutputFile, first_dot_pos);
    append(yaraOptions.outputFile, "_");
    if (file_no < 10)
        append(yaraOptions.outputFile, "0");
    append(yaraOptions.outputFile, std::to_string(file_no));
    append(yaraOptions.outputFile, suffix(options.superOutputFile, first_dot_pos));
}

template <typename TOptions>
inline void set_output_file(TOptions & options, uint32_t const file_no)
{
    uint32_t first_dot_pos = 0;
    for (; first_dot_pos < length(options.superOutputFile); ++first_dot_pos)
    {
        if(options.superOutputFile[first_dot_pos] == '.')
            break;
    }
    // Get the current file name
    options.outputFile = prefix(options.superOutputFile, first_dot_pos);
    append(options.outputFile, "_");
    if (file_no < 10)
        append(options.outputFile, "0");
    append(options.outputFile, std::to_string(file_no));
    append(options.outputFile, suffix(options.superOutputFile, first_dot_pos));
}

// ----------------------------------------------------------------------------
// Function writeMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TOutputFile>
inline void writeMatches(Mapper<TSpec, TConfig> & me, TOutputFile & outputFile)
{
    typedef MapperTraits<TSpec, TConfig>        TTraits;
    typedef MatchesWriter<TSpec, TTraits>       TMatchesWriter;

    start(me.timer);
    TMatchesWriter writer(outputFile,
                          me.suboptimalMatchesSet,
                          me.primaryMatches, me.primaryMatchesProbs, me.cigarSet,
                          me.ctx, me.reads,
                          me.options);
    stop(me.timer);
    me.stats.writeMatches += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cerr << "Output time:\t\t\t" << me.timer << std::endl;
}
// ----------------------------------------------------------------------------
// Function copyMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void copyMatches(TDefaultMapper & target, Mapper<TSpec, TConfig> & source, uint32_t const & contig_offset)
{
    typedef TDefaultMapper::Mapper::Traits         TTraits;

    uint32_t matchCount = length(source.matchesByCoord);
    reserve(target.matchesByCoord, matchCount);
//    resize(target.matchesByCoord, matchCount);
    for (uint32_t i = 0; i<matchCount; ++i)
    {
        TTraits::TMatch currentMatch;
        currentMatch.readId        =source.matchesByCoord[i].readId;
        currentMatch.contigId      =source.matchesByCoord[i].contigId + contig_offset;
        currentMatch.isRev         =source.matchesByCoord[i].isRev;
        currentMatch.contigBegin   =source.matchesByCoord[i].contigBegin;
        currentMatch.contigEnd     =source.matchesByCoord[i].contigEnd;
        currentMatch.errors        =source.matchesByCoord[i].errors;
        appendValue(target.matchesByCoord, currentMatch);
        setSeedErrors(target.ctx, currentMatch.readId, currentMatch.errors);
        setMinErrors(target.ctx, currentMatch.readId, currentMatch.errors);
        setMapped(target.ctx, currentMatch.readId);
    }
}
// ----------------------------------------------------------------------------
// Function mapReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void disMapReads(Mapper<TSpec, TConfig> & me, TDefaultMapper & main_mapper, uint32_t const & contig_offset)
{
    swap(me.reads.seqs, main_mapper.reads.seqs);
    swap(me.reads.names, main_mapper.reads.names);
    me.stats.loadedReads += getReadsCount(me.reads.seqs);
    _disMapReadsImpl(me, main_mapper, me.reads.seqs, contig_offset);
}

template <typename TSpec, typename TConfig, typename TReadSeqs>
inline void _disMapReadsImpl(Mapper<TSpec, TConfig> & me, TDefaultMapper & main_mapper, TReadSeqs & readSeqs, uint32_t const & contig_offset)
{
    initReadsContext(me, readSeqs);
    initSeeds(me, readSeqs);

    collectSeeds<0>(me, readSeqs);
    findSeeds<0>(me, 0);
    classifyReads(me);
    collectSeeds<1>(me, readSeqs);
    collectSeeds<2>(me, readSeqs);
    findSeeds<0>(me, 1);
    findSeeds<0>(me, 2);
    rankSeeds(me);
    reserveMatches(me);
    extendHits<0>(me, 0);
    extendHits<0>(me, 1);
    extendHits<0>(me, 2);
    clearSeeds(me);
    clearHits(me);

    initSeeds(me, readSeqs);
    collectSeeds<1>(me, readSeqs);
    findSeeds<1>(me, 1);
    collectSeeds<2>(me, readSeqs);
    findSeeds<1>(me, 2);
    rankSeeds(me);
    // TODO(esiragusa): filter out hits with distance < 1.
    extendHits<1>(me, 1);
    extendHits<1>(me, 2);
    clearSeeds(me);
    clearHits(me);

    if (me.options.sensitivity > LOW)
    {
        initSeeds(me, readSeqs);
        collectSeeds<2>(me, readSeqs);
        findSeeds<2>(me, 2);
        rankSeeds(me);
        // TODO(esiragusa): filter out hits with distance < 2.
        extendHits<2>(me, 2);
        clearHits(me);
        clearSeeds(me);
    }

    aggregateMatches(me, readSeqs);
//    rankMatches(me, readSeqs);
//    if (me.options.verifyMatches)
//        verifyMatches(me);
//    alignMatches(me);
    swap(me.reads.seqs, main_mapper.reads.seqs);
    swap(me.reads.names, main_mapper.reads.names);
    copyMatches(main_mapper, me, contig_offset);

//    clearMatches(me);
//    clearAlignments(me);
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void runDisMapper(Mapper<TSpec, TConfig> & me, TDefaultMapper & main_mapper, uint32_t const & contig_offset)
{

    configureThreads(me);
    if (me.options.verbose > 1) printRuler(std::cerr);

    loadContigs(me);
    loadContigsIndex(me);

    // Write on the main output file.
    disMapReads(me, main_mapper, contig_offset);
    clearReads(me);
}

template <typename TSpec, typename TConfig>
void copyReads(SeqStore<TSpec, TConfig> & target, SeqStore<TSpec, TConfig> const & source)
{
    typedef SeqStore<TSpec, TConfig>        TSeqStore;
    typedef typename TSeqStore::TSeqs       TSeqs;
    typedef typename Size<TSeqs>::Type      TSeqId;
    typedef typename TSeqStore::TSeqNames   TSeqNames;
    typedef typename Size<TSeqNames>::Type  TNameId;

    TSeqId seqsCount = length(source.seqs);

    reserve(target.seqs, seqsCount, Exact());
    reserve(concat(target.seqs), lengthSum(source.seqs), Exact());
    for (TSeqId seqId = 0; seqId < seqsCount; ++seqId)
    {
        appendValue(target.seqs, source.seqs[seqId]);
    }

    TNameId namesCount = length(source.names);
    reserve(target.names, namesCount, Exact());
    reserve(concat(target.names), lengthSum(source.names), Exact());
    for (TNameId nameId = 0; nameId < namesCount; ++nameId)
    {
        appendValue(target.names, source.names[nameId]);
    }
}


// ----------------------------------------------------------------------------
// Function spawnMapper()
// ----------------------------------------------------------------------------

template <typename TContigsSize, typename TContigsLen, typename TContigsSum,
typename TThreading, typename TSequencing, typename TSeedsDistance>
inline void spawnDisMapper(Options const & options,
                           TDefaultMapper & main_mapper,
                           uint32_t const & contig_offset,
                           TThreading const & /* tag */,
                           TSequencing const & /* tag */,
                           TSeedsDistance const & /* tag */)
{
    typedef ReadMapperConfig<TThreading, TSequencing, TSeedsDistance, TContigsSize, TContigsLen, TContigsSum>  TConfig;

    Mapper<Owner<ConcatDirect<> >, TConfig> mapper(options);
    runDisMapper(mapper, main_mapper, contig_offset);
}

#endif  // #ifndef APP_YARA_DIS_MAPPER_H_
