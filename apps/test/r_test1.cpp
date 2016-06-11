#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/basic.h>

using namespace seqan;

//int estimateReadLength(SeqFileIn &seqFile)
//{
//    if (atEnd(seqFile))
//        return 0;
//
//    typedef String<char, Array<1000> > TBuffer;
//
//    // read chunk into buffer
//    TBuffer buffer;
//    resize(buffer, capacity(buffer));
//    
//    size_t len = seqFile.stream.readsome(&buffer[0], length(buffer));
//    for (size_t i = 0; i < len; ++i)   
//        seqFile.stream.unget();    
//    resize(buffer, len);
//
//    // parse record from buffer    
//    DirectionIterator<TBuffer, Input>::Type iter = directionIterator(buffer, Input());
//    CharString fastaId, seq;
//    readRecord(fastaId, seq, iter, seqFile.format);
//    std::cout << buffer << std::endl << capacity(buffer) << " " << length(buffer) << " length seq = " << length(seq) << std::endl;
//    return length(seq);
//}

int main(int argc, char** argv)
{
    int const N = 8 * 1024 * 1024 ;
    int *a = new int[N];
    int *b = new int[N];
    int window = atof(argv[1]);
    int sum = 0;
    int sum2 = 0;
    int start = 0;
    
    for (int k = 0; k < N; k++) 
    {
        a[k] = k;
    }
     
    for (int k = 0; k < N; k++)
    {
        b[k] = 1;
    }

    for (int k = 0; k < N; k++)
    {
        sum2 += b[k];
    }


    start = sysTime();
    for (int k = 0; k < window; k++ )
    {
        for (int j = k; j < N; j += window )
        {
            a[j] += 1 -j;
        }
    } 
    std::cout << sysTime() - start << " ";

    for (int k = 0; k < N; k++)
    {
        sum += a[k];
    }
   
    std::cout << window << " " << sum << " " << sum2 << std::endl; 
        
}
