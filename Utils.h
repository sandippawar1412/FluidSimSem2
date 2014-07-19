/*
 * @file:   Utils.h
 * @author: Christpher Batty et. al.
 * @brief: code (Selected functions, nicely written)
 *
 */

#ifndef UTILS_H
#define UTILS_H

#include<climits>
using namespace std;

/**
 * Transforms even the sequence 0,1,2,3,... into reasonably good random numbers
 * Challenge: improve on this in speed and "randomness"!
 * This seems to pass several statistical tests, and is a bijective map (of 32-bit unsigned ints)
 */
inline unsigned int randhash(unsigned int seed) {
    unsigned int i = (seed^0xA3C59AC3u)*2654435769u;
    i ^= (i >> 16);
    i *= 2654435769u;
    i ^= (i >> 16);
    i *= 2654435769u;
    return i;
}

/**
 * the inverse of randhash
 * @param h
 * @return
 */
inline unsigned int unhash(unsigned int h) {
    h *= 340573321u;
    h ^= (h >> 16);
    h *= 340573321u;
    h ^= (h >> 16);
    h *= 340573321u;
    h ^= 0xA3C59AC3u;
    return h;
}

/**
 *
 * @param seed
 * @return repeatable stateless pseudo-random number in [0,1] as double
 */
inline double randhashd(unsigned int seed) {
    return randhash(seed) / (double) UINT_MAX;
}

/**
 *
 * @param seed
 * @return repeatable stateless pseudo-random number in [0,1] in float
 */
inline float randhashf(unsigned int seed) {
    return randhash(seed) / (float) UINT_MAX;
}

/**
 *
 * @param seed
 * @param a
 * @param b
 * @return repeatable stateless pseudo-random number in [a,b] as double
 */
inline double randhashd(unsigned int seed, double a, double b) {
    return (b - a)*randhash(seed) / (double) UINT_MAX + a;
}

/**
 *
 * @param seed
 * @param a
 * @param b
 * @return repeatable stateless pseudo-random number in [a,b] as float
 */
inline float randhashf(unsigned int seed, float a, float b) {
    return ( (b - a) * randhash(seed) / (float) UINT_MAX + a);
}

/**
 * This function also handles the boundary cases  
 * at ni and nj indices its takes the value of the cells ni-2 and then its used to
 * average out with ni-1 and ni-2 index cells, (must be done by other function)
 *
 * @param x
 * @param i
 * @param f
 * @param i_low
 * @param i_high
 */
template<class T>
inline void get_barycentric(T x, int& i, T& f, int i_low, int i_high) {
    T s = std::floor(x);
    i = (int) s;
    if (i < i_low) {
        // read the explanation in the if else case below..
        i = i_low;
        f = 0;
    }
    else if (i > i_high - 2) {
        // this thing handles the boundary condition when we try to cross the
        // grid boundaries beyond the ni,nj limit..
        i = i_high - 2;
       
        // f = 1 is important in sense even if it is i_high-2, the value
        // taken would be i_high-1
        // a reasonable choice.. no need to extrapolate, internal values
        // are better.. but then, velocity extrapolation can also been done
        // separately
        f = 1;
    }
    else
        f = (T) (x - s);
        
      //  cout<<"done\n";
}

/**
 *
 * @param value0
 * @param value1
 * @param f
 * @return
 */
template<class S, class T>
inline S lerp(const S& value0, const S& value1, T f) {
    return (1 - f)*value0 + f*value1;
}

/**
 *
 * @param v00
 * @param v10
 * @param v01
 * @param v11
 * @param fx
 * @param fy
 * @return
 */
template<class S, class T>
inline S bilerp(const S& v00, const S& v10,
const S& v01, const S& v11,
T fx, T fy) {
    return lerp(lerp(v00, v10, fx),
            lerp(v01, v11, fx),
            fy);
}

/**
 *
 * @param v000
 * @param v100
 * @param v010
 * @param v110
 * @param v001
 * @param v101
 * @param v011
 * @param v111
 * @param fx
 * @param fy
 * @param fz
 * @return
 */
template<class S, class T>
inline S trilerp(const S& v000, const S& v100,
const S& v010, const S& v110,
const S& v001, const S& v101,
const S& v011, const S& v111,
T fx, T fy, T fz) {
    return lerp(bilerp(v000, v100, v010, v110, fx, fy),
            bilerp(v001, v101, v011, v111, fx, fy),
            fz);
}

/**
 *
 * @param v0000
 * @param v1000
 * @param v0100
 * @param v1100
 * @param v0010
 * @param v1010
 * @param v0110
 * @param v1110
 * @param v0001
 * @param v1001
 * @param v0101
 * @param v1101
 * @param v0011
 * @param v1011
 * @param v0111
 * @param v1111
 * @param fx
 * @param fy
 * @param fz
 * @param ft
 * @return
 */
template<class S, class T>
inline S quadlerp(const S& v0000, const S& v1000,
const S& v0100, const S& v1100,
const S& v0010, const S& v1010,
const S& v0110, const S& v1110,
const S& v0001, const S& v1001,
const S& v0101, const S& v1101,
const S& v0011, const S& v1011,
const S& v0111, const S& v1111,
T fx, T fy, T fz, T ft) {
    return lerp(trilerp(v0000, v1000, v0100, v1100, v0010, v1010, v0110, v1110, fx, fy, fz),
            trilerp(v0001, v1001, v0101, v1101, v0011, v1011, v0111, v1111, fx, fy, fz),
            ft);
}


#endif  /* UTILS_H */

