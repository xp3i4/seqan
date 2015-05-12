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

#ifndef SEQAN_HEADER_SHAPE_MINIMIZER_H
#define SEQAN_HEADER_SHAPE_MINIMIZER_H

namespace seqan
{

// ----------------------------------------------------------------------------
// Struct MinimizerShape
// ----------------------------------------------------------------------------

template <unsigned TSPAN, unsigned TWEIGHT, typename TSpec = void>
struct MinimizerShape;

// ----------------------------------------------------------------------------
// Class Shape<MinimizerShape>
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
class Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> >
{
public:
    typedef typename Value<Shape>::Type THashValue;

    unsigned span;
    unsigned weight;
    unsigned hPointer;
    THashValue hValue;      //minimizer hash value;
    THashValue u_hValue;    //ungapped hash value;
    THashValue hashBuffer[TSPAN - TWEIGHT + 1];
    static const THashValue leftFactor = Power<ValueSize<TValue>::VALUE, TSPAN - 1>::VALUE;
    static const THashValue m_leftFactor = Power<ValueSize<TValue>::VALUE, TWEIGHT - 1>::VALUE;
    TValue  m_leftChar; 
    TValue  leftChar;

    Shape():
        span(TSPAN),
        weight(TWEIGHT)
    {}
};

// ----------------------------------------------------------------------------
// Metafunction LENGTH
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
struct LENGTH<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >
{
    enum { VALUE = TSPAN };
};

// ----------------------------------------------------------------------------
// Metafunction WEIGHT
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
struct WEIGHT<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >
{
    enum { VALUE = TWEIGHT };
};

// ----------------------------------------------------------------------------
// Function weight()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
inline SEQAN_HOST_DEVICE
typename Size< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
weight(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > const &me)
{
    return me.weight;
}

// ----------------------------------------------------------------------------
// Function _minHash()
// ----------------------------------------------------------------------------
// return lexicographically smaller hash as the minimizer

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type 
_minHash(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const & it)
{
    typedef typename Value<Shape<TValue, UngappedShape<TWEIGHT> > >::Type   THValue;
  
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    SEQAN_ASSERT_GT((unsigned)me.span, (unsigned)me.weight); 

    Shape<TValue, UngappedShape<TSPAN> > u_tmpShape;
    Shape<TValue, UngappedShape<TWEIGHT> > tmpShape;
   
    THValue miniTmp;
    me.hValue = hash(tmpShape, it); 
    me.hashBuffer[0] = tmpShape.hValue;
    for (unsigned k = 1; k < TSPAN - TWEIGHT + 1; k++ )
    {
        me.hashBuffer[k] = hashNext(tmpShape, it + k);
        if (me.hValue > me.hashBuffer[k])
            me.hValue = me.hashBuffer[k];
    }
    me.m_leftChar = *(it + TSPAN - TWEIGHT);
    me.leftChar = *(it + TSPAN - 1);
    me.hPointer = TSPAN - TWEIGHT;
    me.u_hValue = hash(u_tmpShape, it);

    return me.hValue;
}

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
hash(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const & it)
{
     return _minHash(me, it);
}


// ----------------------------------------------------------------------------
// Function hashNext()
// ----------------------------------------------------------------------------


template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline void hashInit(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    Shape<TValue, UngappedShape<TSPAN> > tmpShape;
    Shape<TValue, UngappedShape<TWEIGHT> > u_shape;
    hashInit(u_shape, it); 
    for(unsigned k = 0; k < TSPAN - TWEIGHT + 1; k++) 
    {
        me.hashBuffer[k] = hashNext(u_shape, it + k);
    }
    me.m_leftChar = *(it + TSPAN - TWEIGHT -1);
    me.hPointer = TSPAN - TWEIGHT - 1;
    me.u_hValue = hashInit(tmpShape, it);
    me.leftChar = 0;
}

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type 
hashNext(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    typedef typename Value<Shape<TValue, UngappedShape<TWEIGHT> > >::Type THValue;
    typedef typename Size<Shape<TValue, UngappedShape<TWEIGHT> > >::Type TSize;
    Shape<TValue, UngappedShape<TSPAN> > u_shape;
    unsigned prehPointer = me.hPointer;
    me.hPointer = (me.hPointer + 1) % (TSPAN - TWEIGHT + 1);

    me.hashBuffer[me.hPointer] =  (me.hashBuffer[prehPointer] - ordValue(me.m_leftChar) 
                                * (THValue)me.m_leftFactor) * ValueSize<TValue>::VALUE + ordValue((TValue) * (it + (TSize)TSPAN - 1));

    me.m_leftChar = *(it + TSPAN - TWEIGHT);
    me.hValue = me.hashBuffer[0];
    for(unsigned k = 0; k < TSPAN - TWEIGHT + 1; k++)
    {
        if (me.hValue > me.hashBuffer[k]) 
            me.hValue = me.hashBuffer[k];
    }
    me.u_hValue = (me.u_hValue - ordValue(me.leftChar) * (THValue)me.leftFactor) * ValueSize<TValue>::VALUE 
                  + ordValue((TValue)*(it + ((TSize)TSPAN -1)));
    me.leftChar = *it;
    return me.hValue;
} 

}	// namespace seqan

#endif
