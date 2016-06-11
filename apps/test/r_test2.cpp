#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <iostream>
#include <math.h>
#include <seqan/basic.h>

using namespace seqan;

int main(int argc, char** argv)
{
    FragmentStore<> fragStore;
    if(!loadContigs(fragStore, argv[1]))
        return 1;
    if(!loadReads(fragStore, argv[2]))
        return 1;
    
    typedef FragmentStore<>::TReadSeqStore TReadSeqStore;
    typedef Index<TReadSeqStore, IndexQGram<UngappedShape<6> > > TIndex;
    TIndex index(fragStore.readSeqStore);
    Shape<Dna5, UngappedShape<6> > myShape;

    for (unsigned k = 0; k < length(fragStore.contigStore); ++k){
        int count = 0;
        int count22 = 0;
        int m = 0;
        hash(myShape, begin(fragStore.contigStore[k].seq));
        unsigned p = length(fragStore.contigStore[k].seq) - length(myShape) + 1;
        for (unsigned i = 1; i < p; ++i)
        {
            count += countOccurrences(index, myShape);
            count22 += countOccurrences(index, myShape);
            if(m++ == 4000)
            {
                std::cout << (m-1)/4000 << " " << count22 << std::endl;
                count22 = 0; 
                m = 0;
            }
            hashNext(myShape, begin(fragStore.contigStore[k].seq) + i);
        }
        
        std::cout << "count = " << count << std::endl;
    }
        
}

