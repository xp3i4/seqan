#include <seqan/seq_io.h>

static const uint32_t filterMetadataSize = 256;
static const uint8_t INT_WIDTH = 0x40;

#include "bloom_filter.h"
#include "kdx_filter.h"

using namespace seqan;

int main()
{
    //SeqAnBloomFilter<> filter (10, 1, 12, 16777216);
    SeqAnKDXFilter<> filter (10, 5, 167772170);

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

    filter.addFastaFile(file0, 0);
    filter.addFastaFile(file1, 1);
    filter.addFastaFile(file2, 2);
    filter.addFastaFile(file3, 3);
    filter.addFastaFile(file4, 0);
    filter.addFastaFile(file5, 5);
    filter.addFastaFile(file6, 6);
    filter.addFastaFile(file7, 7);
    filter.addFastaFile(file8, 8);
    filter.addFastaFile(file9, 9);

    return 0;
}
