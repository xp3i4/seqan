// ==========================================================================
//                                 update_filter.cpp
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

#define UPDATE_FILTER
// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------
#include <string>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <future>
#include <thread>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------
#include <seqan/index.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------
#include "store_seqs.h"
#include "misc_timer.h"
#include "misc_types.h"
#include "bits_matches.h"
#include "misc_options.h"
#include "misc_options_dis.h"
#include "bloom_filter.h"

using namespace seqan;

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    CharString                      filterFile;
    std::map<uint32_t, CharString>  binContigs;

    uint32_t                        kmerSize;
    uint32_t                        numberOfBins;
    uint64_t                        bloomFilterSize;
    uint32_t                        numberOfHashes;

    unsigned                        threadsCount;
    bool                            verbose;

    Options() :
    kmerSize(20),
    numberOfBins(64),
    bloomFilterSize(8589934592), // 1GB
    numberOfHashes(4),
    threadsCount(1),
    verbose(false)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function getBinNoFromFile()
// ----------------------------------------------------------------------------
bool getBinNoFromFile(uint32_t & binNo, CharString const & currentFile)
{
    CharString fileName = getFilename(currentFile);
    CharString fileNameNoext = trimExtension(fileName);
//    std::string _binNo = toCString(fileNameNoext);

    char* p;
    binNo = std::strtol(toCString(fileNameNoext), &p, 10);
    return *p == 0;
}

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "yara_update_filter");
    setShortDescription(parser, "Update Bloom Filter for Distributed Yara");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIBLOOM-FILTER FILE \\fP> <\\fI4.fna\\fP> <\\fI7.fna\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "BLOOM FILTER"));
    setValidValues(parser, 0, "bf");
    setHelpText(parser, 0, "The path of the bloom filter to be updated.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "FASTA FILES", true));
    setValidValues(parser, 1, SeqFileIn::getFileExtensions());
    setHelpText(parser, 1, "The fasta files of the bins to updated. File names should be exactly the same us bin number (0-indexing). e.g. 0.fna");


    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use (valid for bloom filter only).", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "2048");
    setDefaultValue(parser, "threads", options.threadsCount);

    addOption(parser, ArgParseOption("v", "verbose", "Displays verbose output."));
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Parse verbose output option.
    getOptionValue(options.verbose, parser, "verbose");

    // Parse bloom filter path.
    getArgumentValue(options.filterFile, parser, 0);

    if (isSet(parser, "threads")) getOptionValue(options.threadsCount, parser, "threads");

    // std::map<uint32_t, CharString>  binContigs;
    // Parse read input files.
    uint32_t updateCount = getArgumentValueCount(parser, 1);
    for (uint32_t i = 0; i < updateCount; ++i)
    {
        CharString  currentFile;
        uint32_t    currentBinNo;

        getArgumentValue(currentFile, parser, 1, i);
        
        if (getBinNoFromFile(currentBinNo, currentFile))
            options.binContigs[currentBinNo] = currentFile;
        else
        {
            std::cerr << "File: " << currentFile << "\ndoesn't have a valid name\n";
            exit(1);
        }

    }
    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function verifyFnaFiles()
// ----------------------------------------------------------------------------
inline bool verifyFnaFiles(std::map<uint32_t,CharString> const & fileList)
{
    for (auto fastaFile : fileList)
    {
        if(!verifyFnaFile(fastaFile.second))
            return false;
    }
    return true;
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // verify all the new fasta files
    if (!verifyFnaFiles(options.binContigs))
        return 1;

    try
    {
        SeqAnBloomFilter<> bf(toCString(options.filterFile));

        // clear the bins to updated;
        std::vector<uint32_t> bins2update = {};
        typedef std::map<uint32_t,CharString>::iterator mapIter;
        for(mapIter iter = options.binContigs.begin(); iter != options.binContigs.end(); ++iter)
        {
            bins2update.push_back(iter->first);
        }

        bf.clearBins(bins2update, options.threadsCount);

        Semaphore thread_limiter(options.threadsCount);
        std::vector<std::future<void>> tasks;

        Timer<double>       timer;
        Timer<double>       globalTimer;
        start (timer);
        start (globalTimer);

        // add the new kmers from the new files
        //iterate over the maps
        for(mapIter iter = options.binContigs.begin(); iter != options.binContigs.end(); ++iter)
        {
            tasks.emplace_back(std::async([=, &thread_limiter, &bf] {
                Critical_section _(thread_limiter);
                Timer<double>       binTimer;

                start (binTimer);
                bf.addFastaFile(iter->second, iter->first);
                stop(binTimer);

                if (options.verbose)
                {
                    mtx.lock();
                    std::cerr <<"[bin " << iter->first << "] updated using " << iter->second << "!\t\t" << binTimer << std::endl;
                    mtx.unlock();
                }
            }));
        }
        for (auto &&task : tasks)
        {
            task.get();
        }
        stop(timer);
        if (options.verbose)
            std::cerr <<"All given bins are updated!\t\t" << timer << std::endl;

        start(timer);
        bf.save(toCString(options.filterFile));
        stop(timer);
        if (options.verbose)
            std::cerr <<"Done saving filter (" << bf.size_mb() <<" MB)\t\t" << timer << std::endl;
        stop(globalTimer);
        std::cerr <<"\nFinshed in \t\t\t" << globalTimer << std::endl;
   }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
