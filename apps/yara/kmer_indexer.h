// ==========================================================================
//                                 kmer_indexer.h
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

std::mutex mtx;

using namespace seqan;
class Semaphore
{
    std::mutex m;
    std::condition_variable cv;
    int count;

public:
    Semaphore(int n) : count{n} {}
    void notify()
    {
        std::unique_lock<std::mutex> l(m);
        ++count;
        cv.notify_one();
    }
    void wait()
    {
        std::unique_lock<std::mutex> l(m);
        cv.wait(l, [this]{ return count!=0; });
        --count;
    }
};

class Critical_section
{
    Semaphore &s;
public:
    Critical_section(Semaphore &ss) : s{ss} { s.wait(); }
    ~Critical_section() { s.notify(); }
};


template <typename TShapeValue, typename TIndex>
inline std::vector <TShapeValue> get_kmers_impl(TIndex & ref_fm_index)
{
    typename Iterator<TIndex, TopDown<ParentLinks< > > >::Type iter(ref_fm_index);

    std::vector <TShapeValue> kmers_vector;

    Shape<Dna, UngappedShape<K_MER_LENGTH> > kmer_shape;
    DnaString kmer;
    resize(kmer, K_MER_LENGTH, 'T');

    do
    {
        while(repLength(iter) < K_MER_LENGTH)
        {

            if ( !goDown(iter) && !goRight(iter))
            {
                while(goUp(iter) && !goRight(iter)) ;
            }
            kmer[repLength(iter)-1] = value(iter).lastChar;
        }

        kmers_vector.push_back(hash(kmer_shape, begin(kmer)));

        do
        {
            if(goRight(iter))
            {
                kmer[repLength(iter)-1] = value(iter).lastChar;
                break;
            }
        } while (goUp(iter));

    }
    while (!isRoot(iter));

    return kmers_vector;
}

template <typename TShapeValue, typename TContigsSize, typename TContigsLen, typename TContigsSum>
inline std::vector <TShapeValue> get_kmers(CharString & file_name)
{

    typedef YaraFMConfig<TContigsSize, TContigsLen, TContigsSum>        TIndexConfig;
    typedef FMIndex<void, TIndexConfig>                             TIndexSpec;
    typedef Index<typename TIndexConfig::Text, TIndexSpec>          TIndex;


    TIndex ref_fm_index;


    if (!open(ref_fm_index, toCString(file_name), OPEN_RDONLY))
        throw "ERROR: Could not open the index.";

    return get_kmers_impl<TShapeValue, TIndex>(ref_fm_index);

}


template <typename TShapeValue, typename TContigsSize, typename TContigsLen>
inline std::vector <TShapeValue> get_kmers(AppOptions & options, CharString & file_name)
{
    if (options.contigsSum <= MaxValue<uint32_t>::VALUE)
        return get_kmers <TShapeValue, TContigsSize, TContigsLen, uint32_t>(file_name);
    else
        return get_kmers <TShapeValue, TContigsSize, TContigsLen, uint64_t>(file_name);
}

template <typename TShapeValue, typename TContigsSize>
inline std::vector <TShapeValue> get_kmers(AppOptions & options, CharString & file_name )
{

    if (options.contigsMaxLength <= MaxValue<uint32_t>::VALUE)
        return get_kmers <TShapeValue, TContigsSize, uint32_t>(options, file_name);
    else
        return get_kmers <TShapeValue, TContigsSize, uint64_t>(options, file_name);
}


template <typename TShapeValue>
inline std::vector <TShapeValue> get_kmers(AppOptions & options, CharString & file_name)
{

    if (options.contigsSize <= MaxValue<uint8_t>::VALUE)
        return get_kmers <TShapeValue, uint8_t>(options, file_name);
    else if (options.contigsSize <= MaxValue<uint16_t>::VALUE)
        return get_kmers <TShapeValue, uint16_t>(options, file_name);
    else
        return get_kmers <TShapeValue, uint32_t>(options, file_name);
}

template <typename TShapeValue, typename TMapValue>
inline bool fill_occurance_table_impl(std::vector <std::vector <TShapeValue>> &kmers_vectors,
                                      std::unordered_map <TShapeValue, TMapValue> &kmer_occurance_table)
{

    std::cout << "Computing the total number of kmers \n";
    std::vector<TShapeValue> all_kmers;
    for(uint32_t file_num=0; file_num < kmers_vectors.size(); file_num++ )
    {
        std::vector<TShapeValue> temp;
        temp.swap(all_kmers);
        std::set_union(temp.begin(), temp.end(), kmers_vectors[file_num].begin(), kmers_vectors[file_num].end(), std::back_inserter(all_kmers));
    }

    std::cout << "Total Number of kmers = " << all_kmers.size() <<"\n";

    // a vector of iterators for all vectors
    std::vector<typename std::vector<TShapeValue>::iterator> iter_vector(NUM_OF_BINS);
    for(uint32_t file_num=0; file_num < kmers_vectors.size(); file_num++ )
        iter_vector[file_num] = kmers_vectors[file_num].begin();

    kmer_occurance_table.reserve(all_kmers.size());

    auto it = kmer_occurance_table.begin();
    for(auto const & kmer: all_kmers)
    {
        TMapValue cur_bitset(0);
        for(uint32_t file_num=0; file_num < NUM_OF_BINS; file_num++ )
        {
            if(*iter_vector[file_num] == kmer)
            {
                cur_bitset.set(file_num);
                ++iter_vector[file_num];
            }
        }
        it = kmer_occurance_table.emplace_hint(it, std::make_pair(kmer, cur_bitset));// [kmer_to_int(kmer)] = occurance_bits;
    }
    return true;
}

template <typename TShapeValue, typename TMapValue>
bool getNext(TShapeValue & kmer, TMapValue & cur_bitset, std::vector<typename std::vector<TShapeValue>::iterator> & iter_vector, std::vector<typename std::vector<TShapeValue>::iterator> & end_iter_vector)
{
    cur_bitset.reset();
    while (iter_vector != end_iter_vector)
    {
        for(uint32_t file_num=0; file_num < NUM_OF_BINS; file_num++ )
        {
            if(*(iter_vector[file_num]) == kmer)
            {
                cur_bitset.set(file_num);
                ++iter_vector[file_num];
            }
        }
        if(cur_bitset.any())
        {
            return true;
        }
        ++kmer;
    }
    return false;
}

template <typename TShapeValue, typename TMapValue>
inline bool fill_occurance_table_impl_new(std::vector <std::vector <TShapeValue>> &kmers_vectors,
                                      std::unordered_map <TShapeValue, TMapValue> &kmer_occurance_table)
{
    // a vector of iterators for all vectors
    std::vector<typename std::vector<TShapeValue>::iterator> iter_vector(NUM_OF_BINS);
    for(uint32_t file_num=0; file_num < kmers_vectors.size(); file_num++ )
        iter_vector[file_num] = kmers_vectors[file_num].begin();

    std::vector<typename std::vector<TShapeValue>::iterator> end_iter_vector(NUM_OF_BINS);
    for(uint32_t file_num=0; file_num < kmers_vectors.size(); file_num++ )
        end_iter_vector[file_num] = kmers_vectors[file_num].end();


    TMapValue cur_bitset(0);
    TShapeValue kmer = 0;
    auto it = kmer_occurance_table.end();
    while(getNext(kmer, cur_bitset, iter_vector, end_iter_vector))
    {
        it = kmer_occurance_table.emplace_hint(it, std::make_pair(kmer, cur_bitset));
        ++kmer;
    }
    return true;
}


template <typename TShapeValue>
inline void index_qGrams(AppOptions const & options, CharString const & file_name, uint32_t const & file_num)
{
    std::mt19937 rng(0xDEADBEEF);

    typedef Dna                                         TAlphabet;
    typedef Alloc<>                                     TSeqSpec;
    typedef Owner<ConcatDirect<Limits<uint64_t> > >     TSeqsSpec;

    typedef String<TAlphabet, TSeqSpec>                 TSeq;
    typedef StringSet<TSeq, TSeqsSpec>                  TSeqs;
    typedef StringSet<CharString, TSeqsSpec>            TNames;

    typedef Shape<TAlphabet, UngappedShape<K_MER_LENGTH> >    TShape;
    typedef IndexQGram<TShape, OpenAddressing>                TIndexSpec;
    typedef Index<TSeqs, TIndexSpec>                           TIndex;


    TNames       ids;
    IupacString  seqs;

    // Get the current file name
    CharString out_file_name = options.output_path;
    append(out_file_name, "_");
    append(out_file_name, std::to_string(file_num));

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(file_name)))
    {
        CharString msg = "Unable to open contigs File: ";
        append (msg, file_name);
        throw toCString(msg);
    }


    readRecords(ids, seqs, seqFileIn);
    std::cout << length(seqs) << "\t" << length(seqs[0]) << std::endl;
    TIndex index(seqs);
    TIndex index2;
    indexRequire(index, QGramDir());

    save2(index, toCString(out_file_name), OPEN_WRONLY | OPEN_CREATE);
    open2(index2, toCString(out_file_name), OPEN_RDONLY);
// ATACGTATAC    1001010011
//    found in bin 0
//    found in bin 1
//    found in bin 4
//    found in bin 6
//    found in bin 9
//    TGTAGGCTCG    1000101001
//    found in bin 0
//    found in bin 3
//    found in bin 5
//    found in bin 9

    DnaString kmer = "TGTAGGCTCG";
    TShape kmer_shape;
    hash(kmer_shape, begin(kmer));
    if (countOccurrences(index2, kmer_shape) > 0)
         std::cout << "found in bin "<< file_num << std::endl;

}



template <typename TShapeValue>
inline void compute_qgram_indices (AppOptions const & options)
{

    Timer<double>                       timer;

    start(timer);

    Semaphore thread_limiter(options.num_threads);
    std::vector<std::future<void>> tasks;

    for (uint32_t file_num = 0; file_num < NUM_OF_BINS; ++file_num )
    {

        tasks.emplace_back(
                           std::async([=, &thread_limiter]
                                      {
                                          Critical_section _(thread_limiter);
                                          mtx.lock();
                                          std::cout << "Indexing bin " << file_num  << " of " << NUM_OF_BINS-1 << "\n";
                                          mtx.unlock();

                                          // Get the current file name
                                          CharString file_name = options.fm_index_path;
                                          append(file_name, std::to_string(file_num));
                                          append(file_name, ".fna");

                                          index_qGrams<TShapeValue>(options, file_name, file_num);

                                          mtx.lock();
                                          std::cout << "Finished indexing bin " << file_num << " of " << NUM_OF_BINS-1 << "\n";
                                          mtx.unlock();
                                      }));
    }

    for (auto &&task : tasks)
    {
        task.get();
    }
    stop(timer);

}


template <typename TShapeValue>
inline std::vector <std::vector <TShapeValue> > get_kmers_from_fm_indices (AppOptions & options)
{

    Timer<double>                       timer;

    start(timer);
    std::vector <std::vector <TShapeValue> > kmers_vectors;
    for (uint32_t file_num = 0; file_num < NUM_OF_BINS; ++file_num )
        kmers_vectors.push_back(std::vector<TShapeValue>());

    Semaphore thread_limiter(options.num_threads);
    std::vector<std::future<void>> tasks;

    for (uint32_t file_num = 0; file_num < NUM_OF_BINS; ++file_num )
    {

        tasks.emplace_back(
                           std::async([=, &thread_limiter, &kmers_vectors]
                                      {
                                          Critical_section _(thread_limiter);
                                          mtx.lock();
                                          std::cout << "Getting kmers from bin " << file_num  << " of " << NUM_OF_BINS-1 << "\n";
                                          mtx.unlock();

                                          AppOptions thread_options = options;

                                          // Get the current file name
                                          CharString file_name = options.fm_index_path;
                                          append(file_name, std::to_string(file_num));

                                          // Get limits

                                          //<TContigsSize, TContigsLen, TContigsSum>
                                          String<uint64_t>    limits;

                                          CharString contigsLimitFile(file_name);
                                          append(contigsLimitFile, ".txt.size");

                                          if (!open(limits, toCString(contigsLimitFile), OPEN_RDONLY))
                                          {
                                              CharString msg = "Unable to open contigsLimitFile: ";
                                              append (msg, contigsLimitFile);
                                              throw toCString(msg);
                                          }

                                          if (length(limits) != 3)
                                              throw "limits should have a length of 3";

                                          thread_options.contigsMaxLength = limits[0];
                                          thread_options.contigsSize = limits[1];
                                          thread_options.contigsSum = limits[2];

                                          kmers_vectors[file_num] = get_kmers<TShapeValue>(thread_options, file_name);
                                          mtx.lock();
                                          std::cout << "Finished getting kmers from bin " << file_num << " of " << NUM_OF_BINS-1 << "\n";
                                          mtx.unlock();
                                      }));
    }

    for (auto &&task : tasks)
    {
        task.get();
    }
    stop(timer);

    std::cout<<"=========================================" << std::endl;
    std::cout<<"finished getting kmers from all bins. in " << getValue(timer) <<" secs." << std::endl;
    std::cout<<"Number of kmers from individual bins:" << std::endl;
    for(uint32_t file_num=0; file_num < kmers_vectors.size(); file_num++ )
    {
        std::cout<< kmers_vectors[file_num].size() << ", ";
    }
    std::cout << std::endl;
    std::cout<<"=========================================" << std::endl;

    return kmers_vectors;
}

template <typename TShapeValue, typename TMapValue>
inline std::unordered_map <TShapeValue, TMapValue>  get_kmer_occurance_table_impl(AppOptions & options)
{
    std::vector <std::vector <TShapeValue> > kmers_vectors;
    kmers_vectors = get_kmers_from_fm_indices<TShapeValue>(options);

    Timer<double>                       timer;
    start(timer);

    std::unordered_map <TShapeValue, TMapValue> kmer_occurance_table;
    kmer_occurance_table.reserve(INTIAL_ELEMENT_SIZE);
    for(uint32_t file_num=0; file_num < NUM_OF_BINS; file_num++ )
    {
        for(auto kmer : kmers_vectors[0])
            (kmer_occurance_table[kmer]).set(file_num);

        kmers_vectors.erase(kmers_vectors.begin());
        kmers_vectors.shrink_to_fit();
    }

    //    fill_occurance_table_impl_old(kmers_vectors, kmer_occurance_table);
    stop(timer);
    std::cout<<"finished adding kmers from all bins in " << getValue(timer) <<" secs." << std::endl;
    std::cout<<"=========================================" << std::endl;

    return kmer_occurance_table;
}


template <typename TShapeValue, typename TMapValue>
inline std::unordered_map <TShapeValue, TMapValue> get_kmer_occurance_table(AppOptions & options)
{
    return get_kmer_occurance_table_impl<TShapeValue, TMapValue>(options);
}
