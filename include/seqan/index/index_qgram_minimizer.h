// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
//
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_MINIMIZER_H
#define SEQAN_HEADER_INDEX_MINIMIZER_H

namespace SEQAN_NAMESPACE_MAIN
{

struct FibreHashSA_;

typedef Tag<FibreHashSA_> const FibreHashSA;
typedef FibreHashSA         QGramHashSA;

struct Minimizer_;
typedef Tag<Minimizer_> Minimizer;

template < typename TText, typename TShapeSpec, typename TSpec > 
struct Fibre< Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreHashSA>
{
    typedef typename Fibre< Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreShape>::Type TShape; 
    typedef typename Value<TShape>::Type THashValue;
    typedef String<THashValue> Type;
};

template < typename TText_, typename TShapeSpec>
class Index<TText_, IndexQGram<TShapeSpec, Minimizer> > {
public:
    typedef typename Member<Index, QGramText>::Type             TTextMember;
    typedef typename Fibre<Index, QGramText>::Type              TText;
    typedef typename Fibre<Index, QGramSA>::Type                TSA;
    typedef typename Fibre<Index, QGramDir>::Type               TDir;
    typedef typename Fibre<Index, QGramCounts>::Type            TCounts;
    typedef typename Fibre<Index, QGramCountsDir>::Type         TCountsDir;
    typedef typename Fibre<Index, QGramShape>::Type             TShape;
    typedef typename Fibre<Index, QGramBucketMap>::Type         TBucketMap;

    typedef typename Fibre<Index, QGramHashSA>::Type            THashSA; 

    typedef typename Cargo<Index>::Type                         TCargo;
    typedef typename Size<Index>::Type                          TSize;

    TTextMember         text;        // underlying text
    TSA                 sa;            // suffix array sorted by the first q chars
    TDir                dir;        // bucket directory
    TCounts             counts;        // counts each q-gram per sequence
    TCountsDir          countsDir;    // directory for count buckets
    TShape              shape;        // underlying shape
    TCargo              cargo;        // user-defined cargo
    TBucketMap          bucketMap;    // bucketMap table (used by open-addressing index)
    THashSA             hashSA;       //hash value table SA
    TSize               stepSize;    // store every <stepSize>'th q-gram in the index

    Index():
    stepSize(1)
    {}

    Index(Index &other):
    text(other.text),
    sa(other.sa),
    dir(other.dir),
    counts(other.counts),
    countsDir(other.countsDir),
    shape(other.shape),
    cargo(other.cargo),
    stepSize(1) {}

    Index(Index const &other):
    text(other.text),
    sa(other.sa),
    dir(other.dir),
    counts(other.counts),
    countsDir(other.countsDir),
    shape(other.shape),
    cargo(other.cargo),
    stepSize(1) {}

    template <typename TText__>
    Index(TText__ &_text):
    text(_text),
    stepSize(1) {}

    template <typename TText__>
    Index(TText__ const &_text):
    text(_text),
    stepSize(1) {}

    template <typename TText__, typename TShape_>
    Index(TText__ &_text, TShape_ const &_shape):
    text(_text),
    shape(_shape),
    stepSize(1) {}

    template <typename TText__, typename TShape_>
    Index(TText__ const &_text, TShape_ const &_shape):
    text(_text),
    shape(_shape),
    stepSize(1) {}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec>, FibreHashSA>::Type &
getFibre(Index<TText, TSpec> &index, FibreHashSA){
    return index.hashSA;
}

template <typename TText, typename TSpec>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec> const, FibreHashSA>::Type &
getFibre(Index<TText, TSpec> const &index, FibreHashSA){
    return index.hashSA;
}

template <typename TText, typename TSpec>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec>, FibreHashSA>::Type & indexHashSA(Index<TText, TSpec> &index) 
{ 
    return getFibre(index, FibreHashSA());
}

template <typename TText, typename TSpec>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec> const, FibreSA>::Type & indexHashSA(Index<TText, TSpec> const &index) 
{ 
    return getFibre(index, FibreHashSA());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSA, typename TText, typename TShape, typename TDir>
inline void
_qgramRefineSuffixArray(TSA & /* sa */, TText const & /* text */, TShape const & /* shape */, TDir const & /* dir */)
{
}

template <typename TSA, typename TText, typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TDir>
inline void
_qgramRefineSuffixArray(TSA & sa, TText const & text, Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > const & shape, TDir const & dir)
{
    _refineQGramIndex(sa, dir, text, 0, length(shape));
}

template <typename TSA, typename THASHSA, typename TText, typename TShape>
inline void
_qgramFillHashSA(TSA & /* sa */, THASHSA & /*hashSA*/, TText const & /* text */, TShape const & /* shape */)
{
}

template <typename TSA, typename THASHSA, typename TText, typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
inline void
_qgramFillHashSA(TSA & sa, THASHSA &hashSA, TText const & text, Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > const & /*shape*/) 
{
    typedef typename Size<TSA>::Type TSize;
    Shape<TValue, UngappedShape<TSPAN> > tmpShape;
    for (TSize k = 0; k < length(sa); k++)
    {
        hashSA[k] = hash(tmpShape, begin(text) + *(begin(sa) + k));
    } 
}

template <typename TSA, typename THASHSA, typename TString, typename TSpec, typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec2>
inline void
_qgramFillHashSA(TSA & sa, THASHSA &hashSA, StringSet<TString, TSpec> const & stringSet , Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec2> > const & /*shape*/) 
{
    typedef typename Size<TSA>::Type TSize;
    //typedef typename RemoveConst<typename CTSize>::Type TSize;
    Shape<TValue, UngappedShape<TSPAN> > tmpShape;
    TSize k;
    for (k = 0; k < length(sa); k++)
    {
        unsigned seqNo = getValueI1(*(begin(sa) + k));
        unsigned pos =  getValueI2(*(begin(sa) + k));
        hashSA[k] = hash(tmpShape, begin(stringSet[seqNo]) + pos);
    } 
}

template < typename TText, unsigned TSPAN, unsigned TWEIGHT, typename TShapeSpec>
void createQGramIndex(Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT, TShapeSpec>, Minimizer> > &index)
{
    SEQAN_CHECKPOINT
    typedef Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT, TShapeSpec>, Minimizer > > TIndex;
    typename Fibre<TIndex, QGramText>::Type const   &text       = indexText(index);
    typename Fibre<TIndex, QGramSA>::Type           &sa         = indexSA(index);
    typename Fibre<TIndex, QGramDir>::Type          &dir        = indexDir(index);
    typename Fibre<TIndex, QGramShape>::Type        &shape      = indexShape(index);
    typename Fibre<TIndex, QGramBucketMap>::Type    &bucketMap  = index.bucketMap;
    typename Fibre<TIndex, QGramHashSA>::Type       &hashSA     = indexHashSA(index);

    // 1. clear counters
    _qgramClearDir(dir, bucketMap);

    // 2. count q-grams
    _qgramCountQGrams(dir, bucketMap, text, shape, getStepSize(index));

    if (_qgramDisableBuckets(index))
    {
        // 3. cumulative sum
        _qgramCummulativeSum(dir, True());

        // 4. fill suffix array
        _qgramFillSuffixArray(sa, text, shape, dir, bucketMap, getStepSize(index), True());

        // 5. correct disabled buckets
        _qgramPostprocessBuckets(dir);
    }
    else
    {
        // 3. cumulative sum
        _qgramCummulativeSum(dir, False());

        // 4. fill suffix array
        _qgramFillSuffixArray(sa, text, shape, dir, bucketMap, getStepSize(index), False());
    }

    // 5. refine suffix array
    _qgramRefineSuffixArray(sa, text, shape, dir);

    //6. fill hash suffix array
    _qgramFillHashSA(sa, hashSA, text, shape); 
    
}

template <typename TText, typename TShapeSpec>
inline bool indexCreate(
                        Index<TText, IndexQGram<TShapeSpec, Minimizer> > &index,
                        FibreSADir,
                        Default const)
{
    resize(indexSA(index), _qgramQGramCount(index), Exact());
    resize(indexDir(index), _fullDirLength(index), Exact());
    resize(indexHashSA(index), _qgramQGramCount(index), Exact());
    createQGramIndex(index);
    resize(indexSA(index), back(indexDir(index)), Exact());     // shrink if some buckets were disabled
    return true;
}

// ----------------------------------------------------------------------------
// Function getOccurrences(); MinimizerShape
// ----------------------------------------------------------------------------

template <typename TText, unsigned TSPAN, unsigned TWEIGHT, typename TShapeSpec, typename TValue>
inline typename Infix<typename Fibre<Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT, TShapeSpec>, Minimizer> >, FibreSA>::Type const>::Type getOccurrences(Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT, TShapeSpec>, Minimizer> > &index, Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > &shape)
{
    typedef Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, Minimizer> >    TIndex;
    typedef typename Fibre<TIndex, FibreHashSA>::Type                               THashSA;
    typedef typename Iterator<THashSA>::Type                                        THashSAIterator; 
   
    unsigned s_begin  = dirAt(shape.hValue, index);
    unsigned s_end = dirAt(shape.hValue + 1, index);
    THashSAIterator s_hsBegin = begin(indexHashSA(index)) + s_begin;

    unsigned m = 0;
    for (; m < s_end - s_begin; m++)
    {
        if(*(s_hsBegin + m) == shape.u_hValue)
        {
            break;
        }
    }

    unsigned n = m;
    for (; n < s_end - s_begin; n++)
    {
        if(*(s_hsBegin + n) != shape.u_hValue)
        {
            break;
        }
    
    } 
    return infix(indexSA(index), s_begin + m, s_begin + n); 
}

/*
template <typename TText, unsigned TSPAN, unsigned TWEIGHT, typename TShapeSpec, typename TValue>
//inline typename Size< typename Fibre< Index< TText, IndexQGram<TShapeSpec, Minimizer> >, FibreSA>::Type const >::Type 
inline uint64_t countOccurrences(Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT, TShapeSpec>, Minimizer> > &index, Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > &shape)
{
    typedef Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, Minimizer> >    TIndex;
    typedef typename Fibre<TIndex, FibreHashSA>::Type                               THashSA;
    typedef typename Iterator<THashSA>::Type                                        THashSAIterator; 
   
    unsigned s_begin  = dirAt(shape.hValue, index); 
    unsigned s_end = dirAt(shape.hValue + 1, index);
    THashSAIterator s_hsBegin = begin(indexHashSA(index)) + s_begin;

    //return 0;//dirAt(shape.hValue + 1, index) ;//- dirAt(shape.hValue, index);
    uint64_t m = 0, n;
    //m = s_begin; 
    //n = s_end;

    if (*(s_hsBegin + m + 4) < shape.u_hValue) 
        for (m = 0; m < 3; m++)
            if (*(s_hsBegin + m ) == shape.u_hValue)
                break;
    else
        for (m = 5; m < 8; m++)
        {
            if(*(s_hsBegin + m) == shape.u_hValue)
            {
                break;
            }
        }

    n = m;

    for (; n < s_end - s_begin; n++)
    {
        if(*(s_hsBegin + n) != shape.u_hValue)
        {
            break;
        }
    
    } 

    return  m; 
}
*/

template <typename TText, unsigned TSPAN, unsigned TWEIGHT, typename TShapeSpec, typename TValue>
//inline typename Size< typename Fibre< Index< TText, IndexQGram<TShapeSpec, Minimizer> >, FibreSA>::Type const >::Type 
inline uint64_t countOccurrences(Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT, TShapeSpec>, Minimizer> > &index, Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > &shape)
{
    typedef Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, Minimizer> >    TIndex;
    typedef typename Fibre<TIndex, FibreHashSA>::Type                               THashSA;
    typedef typename Iterator<THashSA>::Type                                        THashSAIterator; 
   
    unsigned s_begin  = dirAt(shape.hValue, index); 
    
    unsigned s_end = dirAt(shape.hValue + 1, index);
    THashSAIterator s_hsBegin = begin(indexHashSA(index)) + s_begin;

    //return 0;//dirAt(shape.hValue + 1, index) ;//- dirAt(shape.hValue, index);
    unsigned  m = 0, n;
    //m = s_begin; 
    //n = s_end;

    if (*(s_hsBegin + (s_end - s_begin) / 2) < shape.u_hValue) 
        for (m = 0; m < (s_end - s_begin) / 2; m++)
            if (*(s_hsBegin + m ) == shape.u_hValue)
                break;
    else
        for (m = (s_end - s_begin) / 2;  m < s_end; m++)
        {
            if(*(s_hsBegin + m) == shape.u_hValue)
            {
                break;
            }
        }

/*
    n = m;

    for (; n < s_end - s_begin; n++)
    {
        if(*(s_hsBegin + n) != shape.u_hValue)
        {
            break;
        }
    
    } 
*/
    return  m; 
}

template <typename TText, unsigned TSPAN, unsigned TWEIGHT, typename TShapeSpec, typename TValue>
inline uint64_t countOccurrences(Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT, TShapeSpec>, Minimizer> > &index, Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > &shape, unsigned &s_begin, unsigned &s_end, uint64_t &pre_hash, unsigned s_len)
{
    typedef Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, Minimizer> >    TIndex;
    typedef typename Fibre<TIndex, FibreHashSA>::Type                               THashSA;
    typedef typename Iterator<THashSA>::Type                                        THashSAIterator; 
    
    if (shape.hValue != pre_hash) 
    {
        pre_hash = shape.hValue;
        s_begin  = dirAt(shape.hValue, index); 
        s_end = dirAt(shape.hValue + 1, index);
        s_len = (s_end - s_begin) / 2;
    }
    

    THashSAIterator s_hsBegin = begin(indexHashSA(index)) + s_begin;

    //return 0;//dirAt(shape.hValue + 1, index) ;//- dirAt(shape.hValue, index);
    unsigned  m = 0, n;
    //m = s_begin; 
    //n = s_end;
    
        if (*(s_hsBegin + s_len) < shape.u_hValue) 
            for (; m < s_len; m++)
                if (*(s_hsBegin + m ) == shape.u_hValue)
                    break;
        
        else
            for (m = s_len;  m < 2 * s_len;  m++)
                if(*(s_hsBegin + m) == shape.u_hValue)
                    break;


/*
    n = m;

    for (; n < s_end - s_begin; n++)
    {
        if(*(s_hsBegin + n) != shape.u_hValue)
        {
            break;
        }
    
    } 
*/
    return m;
    
}

/*

template <typename TText, unsigned TSPAN, unsigned TWEIGHT, typename TShapeSpec, typename TValue>
//inline typename Size< typename Fibre< Index< TText, IndexQGram<TShapeSpec, Minimizer> >, FibreSA>::Type const >::Type 
inline uint64_t countOccurrences(Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT, TShapeSpec>, Minimizer> > &index, Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > &shape)
{
    typedef Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, Minimizer> >    TIndex;
    typedef typename Fibre<TIndex, FibreHashSA>::Type                               THashSA;
    typedef typename Iterator<THashSA>::Type                                        THashSAIterator; 
    //typedef typename Size< typename Fibre< Index< TText, IndexQGram<TShapeSpec, Minimizer> >, FibreSA>::Type const >::Type TOcc;
   
    unsigned s_lower = dirAt(shape.hValue, index);
    THashSAIterator it = begin(indexHashSA(index));
    
    uint64_t block_left  = dirAt(shape.hValue, index);  
    uint64_t block_right = dirAt(shape.hValue + 1, index) ;
    uint64_t block_p = block_left;
    //uint64_t s_r = s_end - s_begin;        
    //uint64_t block_begin = s_begin / 8;
    //uint64_t block_end = s_end / 8; 
    //THashSAIterator s_hsBegin = hsBegin + block_p;
    //TOcc 

    unsigned s_upper = dirAt(shape.hValue + 1, index);
    
    if(s_upper - s_lower > 32)
    while (s_lower < s_upper)
    {
        if(*(it + (s_upper + s_lower) / 2) == shape.u_hValue) 
            break;
        if(*(it + (s_upper + s_lower) / 2) > shape.u_hValue)
        {
            s_upper = (s_upper + s_lower) / 2;
            continue;
        }
        if(*(it + (s_upper + s_lower) /2) < shape.u_hValue)
        {
            s_lower = (s_upper + s_lower) / 2;
            continue;
        }
    } 
    else
    {
        unsigned m = 0, n;
        for (; m < s_upper - s_lower; m++)
        {
            if(*(it + s_lower + m) == shape.u_hValue)
            {
                break;
            }
        }

        n = m;
        for (; n < s_upper - s_lower; n++)
        {
            if(*(it + s_lower + n) != shape.u_hValue)
            {
                break;
            }
    
        }
    }    

    return s_upper - s_lower;
}
*/
template < typename TText, typename TShapeSpec>
inline bool open(
    Index< TText, IndexQGram<TShapeSpec, Minimizer> > &index,
    const char *fileName,
    int openMode)
{    
    String<char> name;

    name = fileName;    append(name, ".txt");
    if ((!open(getFibre(index, QGramText()), toCString(name), openMode)) &&
        (!open(getFibre(index, QGramText()), fileName, openMode))) return false;

    name = fileName;    append(name, ".sa");
    if (!open(getFibre(index, QGramSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".dir");
    if (!open(getFibre(index, QGramDir()), toCString(name), openMode)) return false;
    
    name = fileName;    append(name, ".hsa");
    if (!open(getFibre(index, QGramHashSA()), toCString(name), openMode)) return false;

    return true;
}

template < typename TText, typename TShapeSpec>
inline bool open(
    Index< TText, IndexQGram<TShapeSpec, Minimizer> > &index,
    const char *fileName)
{
    return open(index, fileName, OPEN_RDONLY);
}


template < typename TText, typename TShapeSpec >
inline bool save(
    Index< TText, IndexQGram<TShapeSpec, Minimizer> > &index,
    const char *fileName,
    int openMode)
{
    String<char> name;
    std::cout << name << std::endl; 
    name = fileName;    append(name, ".txt");
    if ((!save(getFibre(index, QGramText()), toCString(name), openMode)) &&
        (!save(getFibre(index, QGramText()), fileName, openMode))) return false;
    
    name = fileName;    append(name, ".sa");
    if (!save(getFibre(index, QGramSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".dir");
    if (!save(getFibre(index, QGramDir()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".hsa");
    if(!save(getFibre(index, QGramHashSA()), toCString(name), openMode)) return false;

    return true;
}

template < typename TText, typename TShapeSpec>
inline bool save(
    Index< TText, IndexQGram<TShapeSpec, Minimizer> > &index,
    const char *fileName)
{
    return save(index, fileName, OPEN_WRONLY | OPEN_CREATE);
}


}

#endif //#ifndef SEQAN_HEADER_...

