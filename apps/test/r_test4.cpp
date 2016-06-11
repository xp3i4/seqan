#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <iostream>


using namespace seqan;

int estimateReadLength(SeqFileIn &seqFile)
{
    if (atEnd(seqFile))
        return 0;

    typedef String<char, Array<1000> > TBuffer;

    // read chunk into buffer
    TBuffer buffer;
    resize(buffer, capacity(buffer));
    
    size_t len = seqFile.stream.readsome(&buffer[0], length(buffer));
    for (size_t i = 0; i < len; ++i)   
        seqFile.stream.unget();    
    resize(buffer, len);

    // parse record from buffer    
    DirectionIterator<TBuffer, Input>::Type iter = directionIterator(buffer, Input());
    CharString fastaId, seq;
    readRecord(fastaId, seq, iter, seqFile.format);
    std::cout << buffer << std::endl << capacity(buffer) << " " << length(buffer) << " length seq = " << length(seq) << std::endl;
    return length(seq);
}

int main(int argc, char** argv)
{
    FragmentStore<> fragStore;
    if(!loadContigs(fragStore, argv[1]))
        return 1;
    if(!loadReads(fragStore, argv[2]))
        return 1;
    #define N 11 
    typedef FragmentStore<>::TContigStore TContigStore;
    typedef FragmentStore<>::TReadSeqStore TReadSeqStore;
    typedef GetValue<TReadSeqStore>::Type TReadSeq;
    typedef GetValue<TContigStore>::Type TContigSeq;

    typedef Index<TReadSeqStore, IndexQGram<Shape<Dna, UngappedShape<4> >, OpenAddressing> > TIndex11;


    typedef Pattern<TIndex11, Pigeonhole<void> > TPattern11;

    typedef Finder<TReadSeq, Pigeonhole<void> > TFinder;
    
    TIndex11 index11(fragStore.readSeqStore);
 
    TPattern11 pattern11(index11);

    float errorRate = 0.05;
    int count = 0;
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i)
    {
        TFinder finder(fragStore.contigStore[i].seq);
            std::cout << "error rate = " << errorRate << std::endl;
    
            count = 0;
            while (find(finder, pattern11, errorRate))
            {
                count++;
                //std::cout << "shape(11)" << std::endl << position(finder) << " " << position(pattern11) << positionRange(finder) <<  std::endl;
            }
            clear(finder); 
            std::cout << "11 " << count << std::endl;
    }
   
    std::cout << count << std::endl;
       //std::cout << count2 << std::endl;
    //std::cout << count3 << std::endl;
    
}
