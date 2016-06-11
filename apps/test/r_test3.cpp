#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <iostream>
#include <math.h>
#include <seqan/basic.h>

using namespace seqan;

float mean(std::vector<int> frq)
{
    int count = 0;
    float sum = 0;
    //int length = sizeof(frq) / sizeof(*frq);
    for (unsigned k = 0; k < frq.size(); k++)
    {
        if(frq[k]!=0)
        {
            sum += k * frq[k];
            count += frq[k];
        }
    }
    return sum / count;
}

float max(std::vector<int> frq)
{
    int max = 0;
    for(unsigned k = 0; k < frq.size(); k++)
    {
        if(frq[k] !=0 )
            max=k;
    }
    return max;
}

float print(std::vector<int> frq)
{
    for(unsigned k = 0; k < frq.size(); k++)
    {
        if (frq[k] != 0)
            std::cout << k << " " << frq[k] << std::endl;
    } 
}

int main(int argc, char** argv)
{
    FragmentStore<> fragStore;
    if(!loadContigs(fragStore, argv[1]))
        return 1;
    if(!loadReads(fragStore, argv[2]))
        return 1;
    
    typedef FragmentStore<>::TReadSeqStore TReadSeqStore;
    typedef GetValue<TReadSeqStore>::Type TReadSeq;
    //typedef Index<TReadSeqStore, IndexQGram<UngappedShape<8>, OpenAddressing > > TIndex;
    typedef Index<TReadSeqStore, IndexQGram<SimpleShape, OpenAddressing > > TIndex;

    TIndex index(fragStore.readSeqStore);

    typedef Pattern<TIndex, Swift<SwiftLocal> > TPattern;
    typedef Finder<TReadSeq, Swift<SwiftLocal> > TFinder;
 

    TPattern pattern(index);
    float errorRate = atof(argv[4]);
    unsigned minLength = atof(argv[3]);
    unsigned kmerLength = atof(argv[5]);
    unsigned const windowLength = 8000;
    unsigned countWindow[1500];
    resize(indexShape(index), kmerLength);
    int count = 0;
    int count_0 = 0;
    int count_n0 = 0;

    std::vector<int> frq;
    for(unsigned k = 0; k < 1000; k++)    
    {
        frq.push_back(0);
    }
    
    double start=sysTime();
    for (unsigned k = 0; k < length(fragStore.contigStore); ++k){
        TFinder finder(fragStore.contigStore[k].seq);
    //    count = 0;
        for (unsigned j = 0; j < 1500; j++)
        countWindow[j] = 0;

        while(find(finder, pattern, errorRate, minLength)) 
        {  
            countWindow[beginPosition(pattern).i2 / windowLength]++;
            count++;
//            std::cout << 2 * k + 1 << " " << beginPosition(pattern).i2 << std::endl;
        }
        //for (unsigned j = 0; j < 1500; j++ )
        //{ 
        //    if (countWindow[j] != 0)
        //        std::cout << 2*k + 1 << " " << j << " " << countWindow[j] << std::endl;
        //}
        if (count == 0)
            count_0++;
        frq[count]++;
        //std::cout << k << " count " << count << std::endl;
        
    }
    std::cout <<  minLength << " " << kmerLength << " " << count << " " <<sysTime() - start << std::endl;
    //print(frq);
    //std::cout << "mean " << mean(frq) << " max " << max(frq) << " count_0 " << count_0 << std::endl;
}

