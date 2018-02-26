#include <seqan/seq_io.h>
#include <vector>
#include <algorithm>
#include <cassert>
static const uint32_t filterMetadataSize = 256;
static const uint8_t INT_WIDTH = 0x40;

#include "seqan/kmer/kmer_base.h"
#include "seqan/kmer/kmer_direct.h"
#include "seqan/kmer/kmer_ibf.h"

using namespace seqan;

int main()
{
    //SeqAnBloomFilter<> filter (10, 1, 12, 16777216);
    uint64_t threads{3};
    uint64_t noBins{10};
    uint64_t kmerSize{3};
    uint64_t hashFunc{3};
    uint64_t bits{11829};

    // ==========================================================================
    // Test constructors
    // ==========================================================================
    std::cout << "Testing ctors" << '\n';

    //typedef InterleavedBloomFilter TSpec;
    typedef DirectAddressing       TSpec;

    // Empty default constructor
    KmerFilter<Dna, TSpec> ctor_empty;
    // Default constructor
    // KmerFilter<Dna, TSpec> ctor_default (noBins, hashFunc, kmerSize, bits);
    /*
    KmerFilter<Dna, TSpec> ctor_default (noBins, hashFunc, kmerSize, bits);
    KmerFilter<Dna, TSpec> ctor_default_helper1 (noBins, hashFunc, kmerSize, bits);
    KmerFilter<Dna, TSpec> ctor_default_helper2 (noBins, hashFunc, kmerSize, bits);
    */

    KmerFilter<Dna, TSpec> ctor_default (noBins, kmerSize);
    KmerFilter<Dna, TSpec> ctor_default_helper1 (noBins, kmerSize);
    KmerFilter<Dna, TSpec> ctor_default_helper2 (noBins, kmerSize);

    // Copy constructor
    KmerFilter<Dna, TSpec> ctor_copy (ctor_default);
    // Copy assignment
    KmerFilter<Dna, TSpec> assignment_copy;
    assignment_copy = ctor_default;
    // Move constructor
    KmerFilter<Dna, TSpec> ctor_move(std::move(ctor_default_helper1));
    // Move assignment
    KmerFilter<Dna, TSpec> assignment_move;
    assignment_move = std::move(ctor_default_helper2);

    // ==========================================================================
    // Test addKmer()
    // ==========================================================================
    std::cout << "Testing addKmer" << '\n';

    for (uint64_t i = 0; i < 10; ++i)
    {
        CharString fasta("../test/");
        append(fasta, CharString(std::to_string(i)));
        append(fasta, CharString(".fasta.gz"));
        if (i % 2)
            addFastaFile(ctor_default, toCString(fasta), i);
        if (i == 0)
        {
            addFastaFile(ctor_copy, toCString(fasta), i);
            addFastaFile(assignment_copy, toCString(fasta), i);
            addFastaFile(ctor_move, toCString(fasta), i);
            addFastaFile(assignment_move, toCString(fasta), i);
        }
    }

    // ==========================================================================
    // Test storing and retrieving
    // ==========================================================================
    std::cout << "Testing store/retrieve" << '\n';

    CharString store1("file.dat");
    store(ctor_default, toCString(store1));
    retrieve(ctor_empty, toCString(store1));

    // ==========================================================================
    // Test whichBins()
    // ==========================================================================
    std::cout << "Testing whichBins" << '\n';

    std::vector<uint64_t> ctor_default_set;

    std::vector<bool> which = whichBins(ctor_default, DnaString("AAA"), 1);
    (void) whichBins(ctor_default, DnaString("AAA"));
    for (uint64_t i = 0; i < which.size(); ++i)
    {
        if (i % 2)
        {
            assert(which[i]);
            ctor_default_set.push_back(i);
        }
        else
            assert(!which[i]);
    }

    // ==========================================================================
    // Test clearBins()
    // ==========================================================================
    std::cout << "Testing clearBins" << '\n';

    // Check if any elements are set in the filters.
    bool ctor_default_any = false;
    for (uint64_t i = 0; i < ctor_default.noOfBlocks * ctor_default.noOfBins; ++i)
    {
        if (ctor_default.filterVector[i])
        {
            ctor_default_any = true;
            break;
        }
    }
    assert(ctor_default_any == true);

    // Reset the filter vectors.
    clearBins(ctor_default, ctor_default_set, threads);

    // Check if filter Vectors are empty.
    ctor_default_any = false;
    for (uint64_t i = 0; i < ctor_default.noOfBlocks * ctor_default.noOfBins; ++i)
    {
        if (ctor_default.filterVector[i])
        {
            ctor_default_any = true;
            break;
        }
    }
    assert(ctor_default_any == false);

    return 0;
}
