// ==========================================================================
//                                 kmer_indexer.cpp
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

struct Options;

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------

#include <mutex>
#include <condition_variable>
#include <future>
#include <thread>
#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <math.h>

#include <limits>
// ----------------------------------------------------------------------------
// cereal headers
// ----------------------------------------------------------------------------

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/bitset.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/archives/binary.hpp>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/parallel.h>
#include <seqan/file.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------
#include "store_seqs.h"
#include "misc_timer.h"
#include "misc_tags.h"
#include "misc_types.h"
#include "bits_matches.h"
#include "misc_options.h"
#include "index_fm.h"
#include "kmer_bits.h"
#include "kmer_index_writer.h"
#include "kmer_indexer.h"

using namespace seqan;

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("kmer_indexer");
    // Set short description, version, and date.
    setShortDescription(parser, "generates a kmer-occurance table from FM-Indices");
    setVersion(parser, "0.1");
    setDate(parser, "June 2017");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIIN\\fP\" \"\\fIOUT\\fP\"");
    addDescription(parser, "generates a kmer-occurance table from FM-Indices");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_PREFIX, "REFERENCE INDICES PREFIX"));
    setHelpText(parser, 0, "An prifix to indexed reference genomes.");
    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUT_PREFIX, "OUTPUT FILE PREFIX"));


    addOption(parser, ArgParseOption("t", "threads", "the number of threads to run in parellel",
                                     ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "32");


    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    if (isSet(parser, "threads"))
        getOptionValue(options.num_threads, parser, "threads");
    getArgumentValue(options.fm_index_path, parser, 0);
    getArgumentValue(options.output_path, parser, 1);

    //modify output file
    append(options.output_path, "_k");
    append(options.output_path, std::to_string(K_MER_LENGTH));
    append(options.output_path, ".kdx");

    return ArgumentParser::PARSE_OK;
}
// --------------------------------------------------------------------------
// Function run_kmer_indexer()
// --------------------------------------------------------------------------
template <typename TShapeValue, typename TMapValue>
inline int run_kmer_indexer(AppOptions & options)
{
    Timer<double>                       timer;
    start(timer);
    std::unordered_map <TShapeValue, TMapValue> kmer_occurance_table;
    kmer_occurance_table = get_kmer_occurance_table<TShapeValue, TMapValue>(options);
    stop(timer);

    std::cout<<"=========================================" << std::endl;
    std::cout<<"finished computing kmer_occurance_table in " << getValue(timer) <<" secs"  << std::endl;
    std::cout<<"Total number of kmers = " << kmer_occurance_table.size() << std::endl;
    std::cout<<"=========================================" << std::endl;


    start(timer);
    write_kmer_index(options, kmer_occurance_table);
    stop(timer);

    std::cout<<"finished writing kmer_occurance_table in " << getValue(timer) <<" secs"  << std::endl ;

    return 0;
}
// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.
int main(int argc, char const ** argv)
{

    // Parse the command line.
    ArgumentParser parser;
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;
    typedef Shape<Dna, UngappedShape<K_MER_LENGTH> >    TShape;
    typedef Value<TShape>::Type                         TShapeValue;

    return run_kmer_indexer <TShapeValue, std::bitset<NUM_OF_BINS> >(options);
}
