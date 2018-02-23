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
    uint64_t t{1};
    std::cout << "Testing ctors" << '\n';
    KmerFilter<Dna, InterleavedBloomFilter> ctor_empty ();
    KmerFilter<Dna, InterleavedBloomFilter> ctor_param (10, 1, 3, 4352);
    KmerFilter<Dna, InterleavedBloomFilter> ctor_insta (ctor_param);
    KmerFilter<Dna, InterleavedBloomFilter> ctor_assig;
    ctor_assig = ctor_insta;

    std::cout << "Testing addKmer" << '\n';
    CharString file0("../test/0.fasta.gz");
    CharString file1("../test/1.fasta.gz");
    CharString file2("../test/2.fasta.gz");
    CharString file3("../test/3.fasta.gz");
    CharString file4("../test/4.fasta.gz");
    CharString file5("../test/5.fasta.gz");
    CharString file6("../test/6.fasta.gz");
    CharString file7("../test/7.fasta.gz");
    CharString file8("../test/8.fasta.gz");
    CharString file9("../test/9.fasta.gz");

    addFastaFile(ctor_param, toCString(file0), 0);
    addFastaFile(ctor_param, toCString(file1), t);
    addFastaFile(ctor_param, toCString(file2), 2);
    addFastaFile(ctor_param, toCString(file3), 3);
    addFastaFile(ctor_insta, toCString(file4), 0);
    addFastaFile(ctor_insta, toCString(file5), 5);
    addFastaFile(ctor_insta, toCString(file6), 6);
    addFastaFile(ctor_assig, toCString(file7), 7);
    addFastaFile(ctor_assig, toCString(file8), 8);
    addFastaFile(ctor_assig, toCString(file9), 9);

    std::cout << "Testing store/retrieve" << '\n';
    KmerFilter<Dna, InterleavedBloomFilter> from_file;

    CharString store1("ctor_param.dat");
    store(ctor_param, toCString(store1));
    retrieve(from_file, toCString(store1));

    CharString store2("ctor_insta.dat");
    store(ctor_insta, toCString(store2));
    retrieve(from_file, toCString(store2));

    CharString store3("ctor_assig.dat");
    store(ctor_assig, toCString(store3));
    retrieve(from_file, toCString(store3));

    std::cout << "Testing clearBins" << '\n';
    // Check if any elements are set in the filters.
    bool ctor_param_any = false;
    for (uint64_t i = 0; i < ctor_param.noOfHashPos; ++i)
    {
        if (ctor_param.filterVector[i])
        {
            ctor_param_any = true;
            break;
        }
    }
    assert(ctor_param_any == true);

    bool ctor_assig_any = false;
    for (uint64_t i = 0; i < ctor_assig.noOfHashPos; ++i)
    {
        if (ctor_assig.filterVector[i])
        {
            ctor_assig_any = true;
            break;
        }
    }
    assert(ctor_assig_any == true);

    bool ctor_insta_any = false;
    for (uint64_t i = 0; i < ctor_insta.noOfHashPos; ++i)
    {
        if (ctor_insta.filterVector[i])
        {
            ctor_insta_any = true;
            break;
        }
    }
    assert(ctor_insta_any == true);

    // Reset the filter vectors.
    std::vector<uint64_t> bins;
    bins.resize(ctor_param.noOfBins, 0);
    for (unsigned i = 0; i < ctor_param.noOfBins; ++i)
        bins[i] = i;
    clearBins(ctor_param, bins, 1);
    clearBins(ctor_insta, bins, t);
    clearBins(ctor_assig, bins, 1);

    // Check if filter Vectors are empty.
    ctor_param_any = false;
    for (uint64_t i = 0; i < ctor_param.noOfHashPos; ++i)
    {
        if (ctor_param.filterVector[i])
        {
            ctor_param_any = true;
            break;
        }
    }
    assert(ctor_param_any == false);

    ctor_assig_any = false;
    for (uint64_t i = 0; i < ctor_assig.noOfHashPos; ++i)
    {
        if (ctor_assig.filterVector[i])
        {
            ctor_assig_any = true;
            break;
        }
    }
    assert(ctor_assig_any == false);

    ctor_insta_any = false;
    for (uint64_t i = 0; i < ctor_insta.noOfHashPos; ++i)
    {
        if (ctor_insta.filterVector[i])
        {
            ctor_insta_any = true;
            break;
        }
    }
    assert(ctor_insta_any == false);

    (void) whichBins(ctor_param, DnaString("AAA"), 1);

    return 0;
}
