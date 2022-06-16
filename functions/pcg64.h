/*
 * PCG Random Number Generation for C.
 *
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 * 
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *     http://www.pcg-random.org
 */

/* This 64 bit minimal C adaptation is by Arunkumar Bupathy <arunbupathy@gmail.com> */

#ifndef PCG_VARIANTS_H_INCLUDED
#define PCG_VARIANTS_H_INCLUDED 1

# include <math.h>
# include <inttypes.h>

typedef __uint128_t uint128_t;
    
#define PCG128(high,low) \
            ((((uint128_t) high) << 64) + low)
            
#define PCG64(high,low) \
            ((((uint64_t) high) << 32) + low)

#define PCG_MULTIPLIER_128 \
        PCG128(2549297995355413924ULL,4865540595714422341ULL)
        
#define PCG_INCREMENT_128  \
        PCG128(6364136223846793005ULL,1442695040888963407ULL)

#define pcg64_rng pcg64_random_t

#define seed_rng pcg64_srandom_r

#define rand_int pcg64_random_r

#define rand_real pcg64_random_r_float53

#define rand_bounded pcg64_boundedrand_r

typedef struct {
    uint128_t state;
    uint128_t inc;
} pcg64_random_t;

#if __cplusplus
extern "C" {
#endif
    
////////////////////////////////////////////////////////////////////////////////

inline uint64_t pcg64_random_r(pcg64_random_t* rng)
{
    uint128_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * PCG_MULTIPLIER_128 + rng->inc;
    // Calculate output function (XSL RR), uses old state for max ILP
    uint64_t value = ((uint64_t) (oldstate >> 64u)) ^ ((uint64_t) oldstate);
    uint32_t rot = oldstate >> 122u;
    return (value >> rot) | (value << ((- rot) & 63));
}

////////////////////////////////////////////////////////////////////////////////

inline uint64_t pcg64_boundedrand_r(pcg64_random_t* rng,
                                       uint64_t bound)
{
    uint64_t threshold = -bound % bound;
    for (;;) {
        uint64_t r = pcg64_random_r(rng);
        if (r >= threshold)
            return r % bound;
    }
}

////////////////////////////////////////////////////////////////////////////////

inline double pcg64_random_r_float53(pcg64_random_t* rng)
{
    return (pcg64_random_r(rng) >> 11) * ldexp(1, -53); // this passed the Crush battery (5 out of 5 runs)
    // return (pcg64_random_r(rng) & ((1ULL << 53) - 1ULL)) * ldexp(1, -53);
    // the second seems to fail Crush (1 in 5 runs, one or two of 144 tests)
    // but to be honest, even the Xoshiro256** fails once in a while, with about the same frequency
    // so no big deal
}

////////////////////////////////////////////////////////////////////////////////

inline void pcg64_srandom_r(pcg64_random_t* rng,
                                     uint128_t initstate, uint128_t initseq)
{
    rng->inc = (initseq << 1u) | 1u;
    rng->state = initstate + rng->inc;
    rng->state = rng->state * PCG_MULTIPLIER_128 + rng->inc;
}

////////////////////////////////////////////////////////////////////////////////

int seed_rng_with_entropy(pcg64_random_t *);

////////////////////////////////////////////////////////////////////////////////

int store_rng_state(pcg64_random_t *, char *);

////////////////////////////////////////////////////////////////////////////////

int read_rng_state(pcg64_random_t *, char *);

////////////////////////////////////////////////////////////////////////////////

#if __cplusplus
}
#endif

#endif
