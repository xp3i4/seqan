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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_MINIMIZER_H
#define SEQAN_HEADER_INDEX_MINIMIZER_H

namespace SEQAN_NAMESPACE_MAIN
{

struct FibreHashSA_;

typedef Tag<FibreHashSA_> const FibreHashSA;
typedef FibreHashSA         QGramHashSA;

typedef Tag<Minimizer_> Minimizer;

template < typename TText, typename TShapeSpec, typename TSpec > 
struct Fibre< Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreHashSA>
{
    typedef typename Fibre< Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreShape>::Type TShape; 
    typedef typename Value<TShape>::Type THashValue;
    typedef String<THashValue> Type;
};

template < typename TText_, typename TShapeSpec, typename TSpec >
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

    THashSA             hashSA;
    TCountsSA           countsSA; 
    TText               textSA;    

    TSize               stepSize;    // store every <stepSize>'th q-gram in the index

    Index():
    stepSize(1) {}

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
// Counting sort - Step 4: Fill suffix array
// w/o constraints
template <
typename TSA,
typename TText,
typename TShape,
typename TDir,
typename TBucketMap,
typename TWithConstraints,
typename TStepSize >
inline void
_qgramFillSuffixArray(
                      TSA &sa,
                      TText const &text,
                      TShape shape,
                      TDir &dir,
                      TBucketMap &bucketMap,
                      TStepSize stepSize,
                      TWithConstraints const)
{
    SEQAN_CHECKPOINT
    typedef typename Iterator<TText const, Standard>::Type    TIterator;
    typedef typename Value<TDir>::Type                        TSize;

    if (empty(shape) || length(text) < length(shape)) return;

    TSize num_qgrams = length(text) - length(shape) + 1;
    TIterator itText = begin(text, Standard());

    if (TWithConstraints::VALUE) {
        TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;            // first hash
        if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = 0;                        // if bucket is enabled
    } else
        sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = 0;            // first hash

    if (stepSize == 1)
        for(TSize i = 1; i < num_qgrams; ++i)
        {
            ++itText;
            if (TWithConstraints::VALUE) {
                TSize bktNo = getBucket(bucketMap, hashNext(shape, itText)) + 1;    // next hash
                if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = i;                    // if bucket is enabled
            } else
                sa[dir[getBucket(bucketMap, hashNext(shape, itText)) + 1]++] = i;    // next hash
        }
    else
        for(TSize i = stepSize; i < num_qgrams; i += stepSize)
        {
            itText += stepSize;
            if (TWithConstraints::VALUE) {
                TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;        // next hash (we mustn't use hashNext here)
                if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = i;                    // if bucket is enabled
            } else
                sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = i;        // next hash
        }
}


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

template < typename TIndex >
void createQGramIndex(TIndex &index)
{
    SEQAN_CHECKPOINT
    typename Fibre<TIndex, QGramText>::Type const   &text       = indexText(index);
    typename Fibre<TIndex, QGramSA>::Type           &sa         = indexSA(index);
    typename Fibre<TIndex, QGramDir>::Type          &dir        = indexDir(index);
    typename Fibre<TIndex, QGramShape>::Type        &shape      = indexShape(index);
    typename Fibre<TIndex, QGramBucketMap>::Type    &bucketMap  = index.bucketMap;

    typename Fibre<TIndex, QGramHashSA>::Type       &hashSA     = indexHashSA(index);
    //typename Fibre<TIndex, QGramCountsSA>::Type     &countsSA   = indexCountsSA(index); 

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

// ----------------------------------------------------------------------------
// Function getOccurrences(); MinimizerShape

//linear search
template <typename TText, unsigned TSPAN, unsigned TWEIGHT, typename TShapeSpec, typename TSpec, typename TValue>
inline typename Infix<typename Fibre<Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT, TShapeSpec>, TSpec> >, FibreSA>::Type const>::Type getOccurrences(Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT, TShapeSpec>, TSpec> > &index,
               Shape<TValue, MinimizerShape<TSPAN, TWEIGHT> > &shape)
{
    typedef Index<TText, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, TSpec> >        TIndex;
    typedef typename Fibre<TIndex, FibreHashSA>::Type                               THashSA;
    typedef typename Iterator<THashSA>::Type                                        THashSAIterator; 
   
    THashSAIterator hsBegin = begin(indexHashSA(index));
    unsigned s_begin  = dirAt(shape.hValue, index);
    unsigned s_end = dirAt(shape.hValue + 1, index);
    unsigned s_r = s_end - s_begin;        
    THashSAIterator s_hsBegin = hsBegin + s_begin;

    unsigned m = 0;
    for (; m < s_r; m++)
    {
        if(*(s_hsBegin + m) == shape.u_hValue)
        {
            break;
        }
    }

    unsigned n = m;
    for (; n < s_r; n++)
    {
        if(*(s_hsBegin + n) != shape.u_hValue)
        {
            break;
        }
    
    } 
    return infix(indexSA(index), s_begin + m, s_begin + n); 
}

}

#endif //#ifndef SEQAN_HEADER_...

