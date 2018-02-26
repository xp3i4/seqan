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

    KmerFilter<Dna, InterleavedBloomFilter> ctor_empty ();
    KmerFilter<Dna, InterleavedBloomFilter> ctor_param (noBins, hashFunc, kmerSize, bits);
    KmerFilter<Dna, InterleavedBloomFilter> ctor_insta (ctor_param);
    KmerFilter<Dna, InterleavedBloomFilter> ctor_assig;
    KmerFilter<Dna, InterleavedBloomFilter> from_file;
/*
    KmerFilter<Dna, DirectAddressing> ctor_empty ();
    KmerFilter<Dna, DirectAddressing> ctor_param (noBins, kmerSize);
    KmerFilter<Dna, DirectAddressing> ctor_insta (ctor_param);
    KmerFilter<Dna, DirectAddressing> ctor_assig;
    KmerFilter<Dna, DirectAddressing> from_file;
*/

    ctor_assig = ctor_insta;

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
            addFastaFile(ctor_param, toCString(fasta), i);
        else
            addFastaFile(ctor_insta, toCString(fasta), i);
        if ( i==1 || i==7)
            addFastaFile(ctor_assig, toCString(fasta), i);
    }

    // ==========================================================================
    // Test storing and retrieving
    // ==========================================================================
    std::cout << "Testing store/retrieve" << '\n';

    CharString store1("ctor_param.dat");
    store(ctor_param, toCString(store1));
    retrieve(from_file, toCString(store1));

    CharString store2("ctor_insta.dat");
    store(ctor_insta, toCString(store2));
    retrieve(from_file, toCString(store2));

    CharString store3("ctor_assig.dat");
    store(ctor_assig, toCString(store3));
    retrieve(from_file, toCString(store3));

    // ==========================================================================
    // Test whichBins()
    // ==========================================================================
    std::cout << "Testing whichBins" << '\n';

    std::vector<uint64_t> ctor_param_set;
    std::vector<uint64_t> ctor_insta_set;
    std::vector<uint64_t> ctor_assig_set;

    std::vector<bool> which = whichBins(ctor_param, DnaString("AAA"), 1);
    (void) whichBins(ctor_param, DnaString("AAA"));
    for (uint64_t i = 0; i < which.size(); ++i)
    {
        if (i % 2)
        {
            assert(which[i]);
            ctor_param_set.push_back(i);
        }
        else
            assert(!which[i]);
    }

    which = whichBins(ctor_insta, DnaString("AAA"), 1);
    (void) whichBins(ctor_insta, DnaString("AAA"));
    for (uint64_t i = 0; i < which.size(); ++i)
    {
        if (!(i % 2))
        {
            assert(which[i]);
            ctor_insta_set.push_back(i);
        }
        else
            assert(!which[i]);
    }

    which = whichBins(ctor_assig, DnaString("AAA"), 1);
    (void) whichBins(ctor_assig, DnaString("AAA"));
    for (uint64_t i = 0; i < which.size(); ++i)
    {
        if ( i==1 || i==7)
        {
            assert(which[i]);
            ctor_assig_set.push_back(i);
        }
        else
            assert(!which[i]);
    }

    // ==========================================================================
    // Test clearBins()
    // ==========================================================================
    std::cout << "Testing clearBins" << '\n';

    // Check if any elements are set in the filters.
    bool ctor_param_any = false;
    for (uint64_t i = 0; i < ctor_param.noOfBlocks * ctor_param.noOfBins; ++i)
    {
        if (ctor_param.filterVector[i])
        {
            ctor_param_any = true;
            break;
        }
    }
    assert(ctor_param_any == true);

    bool ctor_assig_any = false;
    for (uint64_t i = 0; i < ctor_assig.noOfBlocks * ctor_assig.noOfBins; ++i)
    {
        if (ctor_assig.filterVector[i])
        {
            ctor_assig_any = true;
            break;
        }
    }
    assert(ctor_assig_any == true);

    bool ctor_insta_any = false;
    for (uint64_t i = 0; i < ctor_insta.noOfBlocks * ctor_insta.noOfBins; ++i)
    {
        if (ctor_insta.filterVector[i])
        {
            ctor_insta_any = true;
            break;
        }
    }
    assert(ctor_insta_any == true);

    // Reset the filter vectors.
    clearBins(ctor_param, ctor_param_set, threads);
    clearBins(ctor_insta, ctor_insta_set, threads);
    clearBins(ctor_assig, ctor_assig_set, threads);

    // Check if filter Vectors are empty.
    ctor_param_any = false;
    for (uint64_t i = 0; i < ctor_param.noOfBlocks * ctor_param.noOfBins; ++i)
    {
        if (ctor_param.filterVector[i])
        {
            ctor_param_any = true;
            break;
        }
    }
    assert(ctor_param_any == false);

    ctor_assig_any = false;
    for (uint64_t i = 0; i < ctor_assig.noOfBlocks * ctor_assig.noOfBins; ++i)
    {
        if (ctor_assig.filterVector[i])
        {
            ctor_assig_any = true;
            break;
        }
    }
    assert(ctor_assig_any == false);

    ctor_insta_any = false;
    for (uint64_t i = 0; i < ctor_insta.noOfBlocks * ctor_insta.noOfBins; ++i)
    {
        if (ctor_insta.filterVector[i])
        {
            ctor_insta_any = true;
            break;
        }
    }
    assert(ctor_insta_any == false);

    return 0;
}
