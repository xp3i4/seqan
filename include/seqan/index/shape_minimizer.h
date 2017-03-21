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
const float _boundAlpha = 0.8;
struct ReverseComplement_;
typedef Tag<ReverseComplement_> const   ReverseComplement;

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
    unsigned x;
    int first;
    int bound;
    THashValue hValue;     //hash value 
    THashValue XValue;     //minimizer 
    THashValue YValue;     //Y(h,x)
    static const THashValue leftFactor = Power<ValueSize<TValue>::VALUE, TSPAN - 1>::VALUE;
    static const THashValue m_leftFactor = Power<ValueSize<TValue>::VALUE, TWEIGHT - 1>::VALUE;
    TValue  leftChar;

    Shape():
        span(TSPAN),
        weight(TWEIGHT),
        x(0),
        first(-1),
        bound((unsigned)(TWEIGHT * _boundAlpha)),
        hValue(0),
        XValue(0),
        YValue(0)

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
    enum { VALUE = TWEIGHT - TSPAN};
};

// ----------------------------------------------------------------------------
// Metafunction DELTA=LENGTH-WEIGHT
// ----------------------------------------------------------------------------

template <typename T>
struct WINDOW_BIT_SIZE;

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
struct WINDOW_BIT_SIZE<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >
{
    enum { VALUE = 2 * (TSPAN - TWEIGHT)};
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
inline typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
hash(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    Range<TIter> range(it, it + length(me));
    Shape<TValue, UngappedShape<TSPAN> > u_tmpShape;

    return me.hValue;
}

// ----------------------------------------------------------------------------
// Function hashNext()
// ----------------------------------------------------------------------------

uint64_t hash_key;
uint64_t hash_key1;
uint64_t hash_key2;
uint64_t hash_key3;
template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline void hashInit(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    typedef typename Size< Shape<TValue, SimpleShape> >::Type    TSize;

        SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        me.leftChar = 0;
        me.hValue = ordValue(*it);
        hash_key = ((uint64_t)1 << (TSPAN*2 -2 )) - 1;
        hash_key1 = ((uint64_t)1 << TSPAN * 2) - ((uint64_t)1 << TWEIGHT * 2);
        hash_key2 = ((uint64_t) 1 << (TWEIGHT*2)) - 1;
        hash_key3 = ((uint64_t)1 << (TSPAN*2 )) - 1;

        for(TSize i = 2; i < me.span; ++i) {
            //me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*(it +i-2));
            me.hValue = (me.hValue << 2) + ordValue((TValue)*(it +i-2));
        me.x = 0;

        }
}

/*
    template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
    inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
    hashNext(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
    {
    SEQAN_CHECKPOINT
        // remove first, shift left, and add next character
        //typedef typename Value< Shape<TValue, TSpec> >::Type    THValue;
        typedef typename Size< Shape<TValue, TSpec> >::Type        TSize;
        unsigned t;
        SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        me.hValue=((me.hValue & hash_key)<<2)+ordValue((TValue)*(it + ((TSize)me.span - 1)));

        t = (((me.hValue & hash_key1) & ((me.hValue & hash_key1) << 1)) != 0)?__builtin_ctzll((me.hValue & hash_key1 )& ((me.hValue & hash_key1) << 1)):__builtin_ctzll(me.hValue & hash_key1);
        //t += (WINDOW_BIT_SIZE<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::VALUE
        //  - t) & ((WINDOW_BIT_SIZE<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::VALUE - t)
        //  >> (sizeof(unsigned) * CHAR_BIT - 1));
        //t = ((TWEIGHT << 1) < t)?t:(TWEIGHT << 1);
        //std::cout << t << " " << (32 - TWEIGHT <<1) << " " << hash_key1 << std::endl;
        me.XValue = me.hValue << (64 - t) >> (32 - TWEIGHT <<1);//WINDOW_BIT_SIZE<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::VALUE;
        
        me.YValue = (me.hValue >> (t-1) << (64-t))+(me.hValue <<t>>t)+(t<<(TSPAN - TWEIGHT << 1));
        return me.XValue;
    }

*/

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
hashNext(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{

         // remove first, shift left, and add next character
        //typedef typename Value< Shape<TValue, TSpec> >::Type    THValue;
        typedef typename Size< Shape<TValue, TSpec> >::Type        TSize;
        SEQAN_ASSERT_GT((unsigned)me.span, 0u);
        uint64_t v1, t;// tmp;
        //unsigned pre = 23;
 
        me.hValue=((me.hValue & hash_key)<<2)+ordValue((TValue)*(it + ((TSize)me.span - 1)));
        me.XValue = hash_key3;
        //tmp = me.hValue ^ ((me.hValue << 32)+(me.hValue >> 32));
        //tmp = me.hValue ;
                
        //if (me.x == 0 || me.XValue < (me.hValue & hash_key2))
        //{
        for (unsigned k = 64-(TSPAN << 1) ; k <= 64 - (TWEIGHT << 1); k+=2)
        {
            v1 = me.hValue << k >> (64-(TWEIGHT<<1));
            if(me.XValue > v1)
            //if(v > (me.hValue << (64-((TSPAN - k)<< 1)) >> (64 - (pre << 1))))
            {
                //v = (me.hValue << (64-((TSPAN - k)<< 1)) >> (64 - (pre << 1))); 
                me.XValue=v1;
                t=k;
            }
        } 
        //me.XValue = (me.hValue << t >> (64 - (TWEIGHT << 1)));
        me.YValue = (me.hValue >> (64-t) << (64-t-(TWEIGHT<<1)))+(me.hValue & (((uint64_t)1<<(64-t-(TWEIGHT<<1))) - 1))+(t<<((TSPAN - TWEIGHT) << 1));

        return me.XValue; 
}

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
inline uint64_t h2y(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me,  uint64_t h)
{
    uint64_t x = hash_key3, v1, t;
    for (unsigned k = 64-(TSPAN << 1) ; k <= 64 - (TWEIGHT << 1); k+=2)
    {
        v1 = h << k >> (64-(TWEIGHT<<1));
        if(x > v1)
        { 
            x=v1;
            t=k;
        }
    } 
    return (h>> (64-t) << (64-t-(TWEIGHT<<1)))+(h& (((uint64_t)1<<(64-t-(TWEIGHT<<1))) - 1))+(t<<((TSPAN - TWEIGHT) << 1));

}

}
// namespace seqan



#endif
