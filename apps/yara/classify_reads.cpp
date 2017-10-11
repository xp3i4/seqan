// ==========================================================================
//                                 classify_reads.cpp
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
#include <unordered_map>
#include <fstream>

// ----------------------------------------------------------------------------
// cereal headers
// ----------------------------------------------------------------------------

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/bitset.hpp>
#include <cereal/archives/binary.hpp>
// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------
#include "misc_timer.h"
#include "kmer_bits.h"
#include "kmer_index_reader.h"
#include "classify_reads.h"

using namespace seqan;


void setupArgumentParser(ArgumentParser & parser)
{
    setAppName(parser, "classify_reads");
    setShortDescription(parser, "classify reads based on a kmer-index");
    setCategory(parser, "Read Mapping");

    // Setup mandatory arguments.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIKMER INDEX KDX FILE\\fP> <\\fISE-READS FILE\\fP>");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIKMER INDEX KDX FILE\\fP> <\\fIPE-READS FILE 1\\fP> <\\fIPE-READS FILE 2\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "KMER INDEX KDX FILE"));
    setHelpText(parser, 0, "contains a mapping of a kmer and which bins it is pressent at");
    setValidValues(parser, 0, "kdx tsv");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "READS FILE", true));
    setValidValues(parser, 1, SeqFileIn::getFileExtensions());
    setHelpText(parser, 1, "Either one single-end or two paired-end / mate-pair read files.");


    addOption(parser, ArgParseOption("x", "max-errors", "the number of errors that are alowed in the alignment",
                                     ArgParseArgument::INTEGER, "INT"));
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(ClassifyReadsOption & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Parse kmer index input file.
    getArgumentValue(options.kmer_index_file, parser, 0);

    // Parse read input files.
    switch (getArgumentValueCount(parser, 1))
    {
        case 1:
            getArgumentValue(options.reads_file.i1, parser, 1, 0);
            options.singleEnd = true;
            break;
        case 2:
            getArgumentValue(options.reads_file.i1, parser, 1, 0);
            getArgumentValue(options.reads_file.i2, parser, 1, 1);
            options.singleEnd = false;
            break;
        default:
            std::cerr << getAppName(parser) << ": Too many arguments!" << std::endl;
            return ArgumentParser::PARSE_ERROR;
    }
    getOptionValue(options.max_errors, parser, "max-errors");
    options.reload();

    return ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function run_read_classifier()
// --------------------------------------------------------------------------
template <typename TShapeValue, typename TMapValue>
inline int run_read_classifier(ClassifyReadsOption & options)
{
    Timer<double>                       timer;
    start(timer);
    std::unordered_map <TShapeValue, TMapValue> kmer_occurance_table;
    kmer_occurance_table.reserve(INTIAL_ELEMENT_SIZE);

    read_occurance_table(options, kmer_occurance_table);
    stop(timer);

    std::cout<<"finished loading kmer_occurance_table in " << getValue(timer) <<" secs"  << std::endl << std::endl;

    start(timer);
//    print_read_distribution(options, kmer_occurance_table, options.reads_file.i1);
    annotate_reads(options, kmer_occurance_table, options.reads_file.i1);
    stop(timer);

    std::cout<<"finished processing reads in " << getValue(timer) <<" secs"  << std::endl << std::endl;

    if (!options.singleEnd)
    {
        start(timer);
//        print_read_distribution(options, kmer_occurance_table, options.reads_file.i2);
        annotate_reads(options, kmer_occurance_table, options.reads_file.i2);
        stop(timer);

        std::cout<<"finished processing reads (2nd Pair) in " << getValue(timer) <<" secs"  << std::endl << std::endl;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------
//
int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    setupArgumentParser(parser);
    ClassifyReadsOption options;
    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    typedef Shape<Dna, UngappedShape<K_MER_LENGTH> >    TShape;
    typedef Value<TShape>::Type                         TShapeValue;

    return run_read_classifier <TShapeValue, std::bitset<NUM_OF_BINS> >(options);
}
