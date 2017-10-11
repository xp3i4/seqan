// ==========================================================================
//                       Yara - Yet Another Read Aligner (Distributed)
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

#define YARA_MAPPER

// ============================================================================
// Forwards
// ============================================================================

struct Options;

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------

#include <random>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/parallel.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "basic_alphabet.h"
#include "file_pair.h"
#include "file_prefetched.h"
#include "store_seqs.h"
#include "misc_timer.h"
#include "misc_tags.h"
#include "misc_types.h"
#include "index_fm.h"
#include "bits_reads.h"
#include "bits_hits.h"
#include "bits_context.h"
#include "bits_matches.h"
#include "bits_seeds.h"
#include "bits_bucket.h"
#include "find_verifier.h"
#include "find_extender.h"
#include "misc_options.h"
#include "mapper_collector.h"
#include "mapper_classifier.h"
#include "mapper_ranker.h"
#include "mapper_filter.h"
#include "mapper_extender.h"
#include "mapper_verifier.h"
#include "mapper_aligner.h"
#include "mapper_writer.h"
#include "mapper.h"

#include "mapper_dis.h"

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------
void setupArgumentParser(ArgumentParser & parser, SuperOptions const & options)
{
    setAppName(parser, "yara_mapper_dis");
    setShortDescription(parser, "Yara Mapper (distributed)");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    // Setup mandatory arguments.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIREFERENCE INDEX PREFIX\\fP> <\\fIKMER INDEX KDX FILE\\fP> <\\fISE-READS FILE\\fP>");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIREFERENCE INDEX PREFIX\\fP> <\\fIKMER INDEX KDX FILE\\fP> <\\fIPE-READS FILE 1\\fP> <\\fIPE-READS FILE 2\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "REFERENCE INDEX PREFIX"));
    setHelpText(parser, 0, "Prefix to multiple indices of reference genomes.");

//    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "KMER INDEX KDX FILE"));
//    setHelpText(parser, 1, "contains a mapping of a kmer and which bins it is pressent at");
//    setValidValues(parser, 1, "kdx tsv");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "READS FILE", true));
    setValidValues(parser, 1, SeqFileIn::getFileExtensions());
    setHelpText(parser, 1, "Either one single-end or two paired-end / mate-pair read files.");

    addOption(parser, ArgParseOption("v", "verbose", "Displays global statistics."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Displays extensive statistics for each batch of reads."));

    // Setup output options.
    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "output-file", "Specify an output file. Default: write the file to standard output.",
                                     ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output-file", BamFileOut::getFileExtensions());

    addOption(parser, ArgParseOption("f", "output-format", "Specify an output format. Note: when specifying the option \
                                     --output-file, the output format is taken from the filename \
                                     extension.", ArgParseOption::STRING));
    setValidValues(parser, "output-format", getExtensionsWithoutLeadingDot(BamFileOut::getFileExtensions()));
    setDefaultValue(parser, "output-format", "sam");

#if SEQAN_HAS_ZLIB
    addOption(parser, ArgParseOption("u", "uncompressed-bam", "Turn off compression of BAM written to standard output."));
    hideOption(getOption(parser, "uncompressed-bam"));
#endif

    addOption(parser, ArgParseOption("rg", "read-group", "Specify a read group for all records in the SAM/BAM file.",
                                     ArgParseOption::STRING));
    setDefaultValue(parser, "read-group", options.readGroup);

    addOption(parser, ArgParseOption("sa", "secondary-alignments", "Specify whether to output secondary alignments in \
                                     the XA tag of the primary alignment, as separate \
                                     secondary records, or to omit them.",
                                     ArgParseOption::STRING));
    setValidValues(parser, "secondary-alignments", options.secondaryAlignmentsList);
    setDefaultValue(parser, "secondary-alignments", options.secondaryAlignmentsList[options.secondaryAlignments]);

    addOption(parser, ArgParseOption("ra", "rabema-alignments", "Output alignments compatible with the \
                                     Read Alignment BEnchMArk (RABEMA)."));

    // Setup mapping options.
    addSection(parser, "Mapping Options");

    addOption(parser, ArgParseOption("e", "error-rate", "Consider alignments within this percentual number of errors. \
                                     Increase this threshold to increase the number of mapped reads. \
                                     Decrease this threshold to decrease the runtime.",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "error-rate", "0");
    setMaxValue(parser, "error-rate", "10");
    setDefaultValue(parser, "error-rate", 100.0 * options.errorRate);

    addOption(parser, ArgParseOption("s", "strata-rate", "Consider suboptimal alignments within this percentual number \
                                     of errors from the optimal alignment. Increase this threshold to increase \
                                     the number of alternative alignments at the expense of runtime.",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "strata-rate", "0");
    setMaxValue(parser, "strata-rate", "10");
    setDefaultValue(parser, "strata-rate", 100.0 * options.strataRate);

    addOption(parser, ArgParseOption("y", "sensitivity", "Sensitivity with respect to edit distance. \
                                     Full sensitivity guarantees to find all considered alignments \
                                     but increases runtime, low sensitivity decreases runtime by \
                                     breaking such guarantee.",
                                     ArgParseOption::STRING));
    setValidValues(parser, "sensitivity", options.sensitivityList);
    setDefaultValue(parser, "sensitivity", options.sensitivityList[options.sensitivity]);

    // Setup paired-end mapping options.
    addSection(parser, "Paired-End Mapping Options");

    addOption(parser, ArgParseOption("ll", "library-length", "Expected library length. Default: autodetected.",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "library-length", "1");

    addOption(parser, ArgParseOption("ld", "library-deviation", "Deviation from the expected library length. \
                                     Default: autodetected.", ArgParseOption::INTEGER));
    setMinValue(parser, "library-deviation", "1");

    addOption(parser, ArgParseOption("i", "indel-rate", "Rescue unaligned ends within this percentual number of indels.",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "indel-rate", "0");
    setMaxValue(parser, "indel-rate", "50");
    setDefaultValue(parser, "indel-rate", 100.0 * options.indelRate);

    addOption(parser, ArgParseOption("ni", "no-indels", "Turn off the rescue of unaligned ends containing indels."));

    //    addOption(parser, ArgParseOption("lo", "library-orientation", "Expected orientation of the segments in the library.",
    //                                     ArgParseOption::STRING));
    //    setValidValues(parser, "library-orientation", options.libraryOrientationList);
    //    setDefaultValue(parser, "library-orientation", options.libraryOrientationList[options.libraryOrientation]);

    //    addOption(parser, ArgParseOption("la", "anchor", "Anchor one read and verify its mate."));

    // Setup performance options.
    addSection(parser, "Performance Options");

    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
#ifdef _OPENMP
    setMaxValue(parser, "threads", "2048");
#else
    setMaxValue(parser, "threads", "1");
#endif
    setDefaultValue(parser, "threads", options.threadsCount);

    addOption(parser, ArgParseOption("rb", "reads-batch", "Specify the number of reads to process in one batch.",
                                     ArgParseOption::INTEGER));

    addOption(parser, ArgParseOption("x", "max-errors", "the number of errors that are alowed in the alignment",
                                     ArgParseArgument::INTEGER, "INT"));

    setMinValue(parser, "reads-batch", "1000");
    setMaxValue(parser, "reads-batch", "1000000");
    setDefaultValue(parser, "reads-batch", options.readsCount);
    hideOption(getOption(parser, "reads-batch"));
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(SuperOptions & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Parse indexed genome input file.
    getArgumentValue(options.superContigsIndicesFile, parser, 0);

    // Parse kmer index input file.
//    getArgumentValue(options.kmer_index_file, parser, 1);

    // Parse read input files.
    switch (getArgumentValueCount(parser, 1))
    {
        case 1:
            getArgumentValue(options.readsFile.i1, parser, 1, 0);
            options.singleEnd = true;
            break;
        case 2:
            getArgumentValue(options.readsFile.i1, parser, 1, 0);
            getArgumentValue(options.readsFile.i2, parser, 1, 1);
            options.singleEnd = false;
            break;
        default:
            std::cerr << getAppName(parser) << ": Too many arguments!" << std::endl;
            return ArgumentParser::PARSE_ERROR;
    }

    // Parse output file.
    getOptionValue(options.superOutputFile, parser, "output-file");

    // Parse output format.
    CharString outputFormat;
    if (getOptionValue(outputFormat, parser, "output-format"))
    {
        addLeadingDot(outputFormat);
        guessFormatFromFilename(outputFormat, options.outputFormat);
    }
    else
        assign(options.outputFormat, Sam());

#if SEQAN_HAS_ZLIB
    getOptionValue(options.uncompressedBam, parser, "uncompressed-bam");
#endif

    // Parse output options.
    getOptionValue(options.readGroup, parser, "read-group");
    getOptionValue(options.secondaryAlignments, parser, "secondary-alignments", options.secondaryAlignmentsList);
    getOptionValue(options.rabema, parser, "rabema-alignments");

    // Parse mapping options.
    unsigned errorRate;
    if (getOptionValue(errorRate, parser, "error-rate"))
        options.errorRate = errorRate / 100.0;

    unsigned strataRate;
    if (getOptionValue(strataRate, parser, "strata-rate"))
        options.strataRate = strataRate / 100.0;

    getOptionValue(options.sensitivity, parser, "sensitivity", options.sensitivityList);

    // Parse paired-end mapping options.
    getOptionValue(options.libraryLength, parser, "library-length");
    getOptionValue(options.libraryDev, parser, "library-deviation");
    //    getOptionValue(options.libraryOrientation, parser, "library-orientation", options.libraryOrientationList);

    unsigned indelRate;
    if (getOptionValue(indelRate, parser, "indel-rate"))
        options.indelRate = indelRate / 100.0;

    options.verifyMatches = !isSet(parser, "no-indels");

    // Parse performance options.
    getOptionValue(options.threadsCount, parser, "threads");
    getOptionValue(options.readsCount, parser, "reads-batch");


//    getOptionValue(options.max_errors, parser, "max-errors");

    if (isSet(parser, "verbose")) options.verbose = 1;
    if (isSet(parser, "very-verbose")) options.verbose = 2;

    // Get version.
    options.version = getVersion(parser);

    // Get command line.
    for (int i = 0; i < argc; i++)
    {
        append(options.commandLine, argv[i]);
        appendValue(options.commandLine, ' ');
    }
    eraseBack(options.commandLine);
    
    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function configureMapper()
// ----------------------------------------------------------------------------

template <typename TContigsSize, typename TContigsLen, typename TThreading, typename TSequencing, typename TSeedsDistance>
void configureMapper(Options const & options,
                     TDefaultMapper & main_mapper,
                     uint32_t const & contig_offset,
                     TThreading const & threading,
                     TSequencing const & sequencing,
                     TSeedsDistance const & distance)
{
    if (options.contigsSum <= MaxValue<uint32_t>::VALUE)
    {
        spawnDisMapper<TContigsSize, TContigsLen, uint32_t>(options, main_mapper, contig_offset, threading, sequencing, distance);
    }
    else
    {
        spawnDisMapper<TContigsSize, TContigsLen, uint64_t>(options, main_mapper, contig_offset, threading, sequencing, distance);
    }
}

template <typename TContigsSize, typename TThreading, typename TSequencing, typename TSeedsDistance>
void configureMapper(Options const & options,
                     TDefaultMapper & main_mapper,
                     uint32_t const & contig_offset,
                     TThreading const & threading,
                     TSequencing const & sequencing,
                     TSeedsDistance const & distance)
{
    if (options.contigsMaxLength <= MaxValue<uint32_t>::VALUE)
    {
        configureMapper<TContigsSize, uint32_t>(options, main_mapper, contig_offset, threading, sequencing, distance);
    }
    else
    {
#ifdef YARA_LARGE_CONTIGS
        configureMapper<TContigsSize, uint64_t>(options, main_mapper, contig_offset, threading, sequencing, distance);
#else
        throw RuntimeError("Maximum contig length exceeded. Recompile with -DYARA_LARGE_CONTIGS=ON.");
#endif
    }
}

template <typename TThreading, typename TSequencing, typename TSeedsDistance>
void configureMapper(Options const & options,
                     TDefaultMapper & main_mapper,
                     uint32_t const & contig_offset,
                     TThreading const & threading,
                     TSequencing const & sequencing,
                     TSeedsDistance const & distance)
{
    if (options.contigsSize <= MaxValue<uint8_t>::VALUE)
    {
        configureMapper<uint8_t>(options, main_mapper, contig_offset, threading, sequencing, distance);
    }
    else if (options.contigsSize <= MaxValue<uint16_t>::VALUE)
    {
        configureMapper<uint16_t>(options, main_mapper, contig_offset, threading, sequencing, distance);
    }
    else
    {
#ifdef YARA_LARGE_CONTIGS
        configureMapper<uint32_t>(options, main_mapper, contig_offset, threading, sequencing, distance);
#else
        throw RuntimeError("Maximum number of contigs exceeded. Recompile with -DYARA_LARGE_CONTIGS=ON.");
#endif
    }
}

template <typename TThreading, typename TSequencing>
void configureMapper(Options const & options,
                     TDefaultMapper & main_mapper,
                     uint32_t const & contig_offset,
                     TThreading const & threading,
                     TSequencing const & sequencing)
{
    if (options.sensitivity == FULL)
        return configureMapper(options, main_mapper, contig_offset, threading, sequencing, EditDistance());
    else
        return configureMapper(options, main_mapper, contig_offset, threading, sequencing, HammingDistance());
}

template <typename TThreading>
void configureMapper(Options const & options,
                     TDefaultMapper & main_mapper,
                     uint32_t const & contig_offset,
                     TThreading const & threading)
{
    if (options.singleEnd)
        configureMapper(options, main_mapper, contig_offset, threading, SingleEnd());
    else
        configureMapper(options, main_mapper, contig_offset, threading, PairedEnd());
}

void configureMapper(Options const & options,
                     TDefaultMapper & main_mapper,
                     uint32_t const & contig_offset)
{
#ifdef _OPENMP
    if (options.threadsCount > 1)
        configureMapper(options, main_mapper, contig_offset, Parallel());
    else
#endif
        configureMapper(options, main_mapper, contig_offset, Serial());
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------
//
int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    SuperOptions options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    try
    {

        options.outputFile = options.superOutputFile;
        TDefaultMapper main_mapper(options);

        openReads(main_mapper);

        std::vector<uint32_t> contig_offsets(options.NUM_OF_BINS, 0);
        Timer<double> timer;

        start(timer);

        for (uint32_t i=0; i < options.NUM_OF_BINS; ++i)
        {
            set_current_index_file(options, i);
            SeqStore<void, YaraContigsConfig<> >  tempContigs;

            if (!open(tempContigs, toCString(main_mapper.options.contigsIndexFile), OPEN_RDONLY))
                throw RuntimeError("Error while opening reference file.");
            contig_offsets[i] = length(main_mapper.contigs.names);
            append(main_mapper.contigs.names, tempContigs.names);
            append(main_mapper.contigs.seqs, tempContigs.seqs);
        }
        // Open output file and write header.
        openOutputFile(main_mapper);


        // Process reads in blocks.
        // load reads here

        // classify reads here
        // create mappers and run them on subsets
        while (true)
        {
            if (main_mapper.options.verbose > 1) printRuler(std::cerr);
            loadReads(main_mapper);
            initReadsContext(main_mapper, main_mapper.reads.seqs);

            if (empty(main_mapper.reads.seqs)) break;
            for (uint32_t i=0; i < options.NUM_OF_BINS; ++i)
            {
                Options yaraOptions = options;
                set_current_index_file(yaraOptions, options, i);
                if (!openContigsLimits(yaraOptions))
                    throw RuntimeError("Error while opening reference file.");
                configureMapper(yaraOptions, main_mapper, contig_offsets[i]);
            }

            aggregateMatches(main_mapper, main_mapper.reads.seqs);
            rankMatches(main_mapper, main_mapper.reads.seqs);
            if (main_mapper.options.verifyMatches)
                verifyMatches(main_mapper);
            alignMatches(main_mapper);
            writeMatches(main_mapper);
            clearMatches(main_mapper);
            clearAlignments(main_mapper);

            clearReads(main_mapper);
        }
        closeReads(main_mapper);
        closeOutputFile(main_mapper);
        stop(timer);
        if (main_mapper.options.verbose > 0)
            printStats(main_mapper, timer);
     }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
