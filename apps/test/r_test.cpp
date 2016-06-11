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
    double startTime = sysTime();
    FragmentStore<> fragStore;
    if(!loadContigs(fragStore, argv[1]))
        return 1;
    if(!loadReads(fragStore, argv[2]))
        return 1;
    typedef FragmentStore<>::TContigStore TContigStore;
    typedef FragmentStore<>::TReadSeqStore TReadSeqStore;
    typedef GetValue<TReadSeqStore>::Type TReadSeq;
    typedef GetValue<TContigStore>::Type TContigSeq;

    typedef Index<TReadSeqStore, IndexQGram<Shape<Dna5, UngappedShape<3> >, OpenAddressing> > TIndex3;
    typedef Index<TReadSeqStore, IndexQGram<Shape<Dna5, UngappedShape<4> >, OpenAddressing> > TIndex4;
    typedef Index<TReadSeqStore, IndexQGram<Shape<Dna5, UngappedShape<5> >, OpenAddressing> > TIndex5;
    typedef Index<TReadSeqStore, IndexQGram<Shape<Dna5, UngappedShape<6> >, OpenAddressing> > TIndex6;
    typedef Index<TReadSeqStore, IndexQGram<Shape<Dna5, UngappedShape<7> >, OpenAddressing> > TIndex7;
    typedef Index<TReadSeqStore, IndexQGram<Shape<Dna5, UngappedShape<8> >, OpenAddressing> > TIndex8;
    typedef Index<TReadSeqStore, IndexQGram<Shape<Dna5, UngappedShape<9> >, OpenAddressing> > TIndex9;
    typedef Index<TReadSeqStore, IndexQGram<Shape<Dna5, UngappedShape<11> >, OpenAddressing> > TIndex11;


    typedef Pattern<TIndex3, Swift<SwiftSemiGlobal> > TPattern3;
    typedef Pattern<TIndex4, Swift<SwiftSemiGlobal> > TPattern4;
    typedef Pattern<TIndex5, Swift<SwiftSemiGlobal> > TPattern5;
    typedef Pattern<TIndex6, Swift<SwiftSemiGlobal> > TPattern6;
    typedef Pattern<TIndex7, Swift<SwiftSemiGlobal> > TPattern7;
    typedef Pattern<TIndex8, Swift<SwiftSemiGlobal> > TPattern8;
    typedef Pattern<TIndex9, Swift<SwiftSemiGlobal> > TPattern9;
    typedef Pattern<TIndex11, Swift<SwiftSemiGlobal> > TPattern11;

    typedef Finder<TReadSeq, Swift<SwiftSemiGlobal> > TFinder;
     
    
    TIndex3 index3(fragStore.readSeqStore);
    TIndex4 index4(fragStore.readSeqStore);
    TIndex5 index5(fragStore.readSeqStore);
    TIndex6 index6(fragStore.readSeqStore);
    TIndex7 index7(fragStore.readSeqStore);
    TIndex8 index8(fragStore.readSeqStore);
    TIndex9 index9(fragStore.readSeqStore);
    TIndex11 index11(fragStore.readSeqStore);
 
    TPattern3 pattern3(index3);
    TPattern4 pattern4(index4);
    TPattern5 pattern5(index5);
    TPattern6 pattern6(index6);
    TPattern7 pattern7(index7);
    TPattern8 pattern8(index8);
    TPattern9 pattern9(index9);
    TPattern11 pattern11(index11);
    
    TPattern9 patternX(index9);

    float errorRate = 0.1;
    int count = 0;
    int bucket = 1;
    int tmp = 0;
    int const width = 6000;
    int tp = 0;
    #define N 5000
    int cs[N];
    for (unsigned k = 0; k < N; k++)
        cs[k] = 0;
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i)
    {
        TFinder finder(fragStore.contigStore[i].seq);
        for (unsigned j = 0; j < 1; j++) 
        {
           
            count = 0;
            while (find(finder, patternX, errorRate))
            {
                //std::cout << "po" <<  positionRange(finder).i2<< std::endl;
                //tp = positionRange(finder).i1;
                //if (floor(positionRange(finder).i1 / width ) == tmp)
                //    bucket++;
                //else{
                //    tmp = floor(positionRange(finder).i1 / width );
                //    std::cout << tmp << " " << bucket << std::endl;
                //    bucket = 1;
                //}
                //count++;
                //std::cout << "shape(4)" << std::endl << position(finder) << " " << position(pattern4) << positionRange(finder) << std::endl;
                cs[(int)floor(positionRange(finder).i1) / width]+=1;
            }
            for (unsigned k = 0; k < N; k++)
            {
                if(cs[k] >0)
                {
                    std::cout << k << " " << cs[k] << std::endl;
                    count++;
                }
            }
            //std::cout << tmp << " " << bucket << std::endl;
            //std::cout << "bucketwidth " << finder.hits[0].bucketWidth << std::endl;
            clear(finder); 
            std::cout << "4 " << count << std::endl;

            //count = 0;
            //while (find(finder, pattern5, errorRate))
            //{
            //    count++;
            //    //std::cout << "shape(5)" << std::endl << position(finder) << " " << position(pattern5) << positionRange(finder) << std::endl;
            //}
            //clear(finder);
            //std::cout << "5 " << count << std::endl;

            //count = 0;
            //while (find(finder, pattern9, errorRate))
            //{
            //    count++;
            //    //std::cout << "shape(9)" << std::endl << position(finder) << " " << position(pattern9) << positionRange(finder) << std::endl;
            //}
            //clear(finder); 
            //std::cout << "9 " << count << std::endl;

            //count = 0;
            //while (find(finder, pattern11, errorRate))
            //{
            //    count++;
            //    //std::cout << "shape(11)" << std::endl << position(finder) << " " << position(pattern11) << positionRange(finder) <<  std::endl;
            //}
            //clear(finder); 
            //std::cout << "11 " << count << std::endl;

        errorRate += 0.01;
        }
    }
   
    std::cout << count << "time = " << sysTime() - startTime << std::endl;
    //std::cout << count2 << std::endl;
    //std::cout << count3 << std::endl;
     
}
