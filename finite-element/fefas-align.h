#ifndef _fefas_align_h
#define _fefas_align_h

#if defined(__AVX__) || defined(__xlc__)  // Assume these compilers support __attribute__((aligned(32)))
#  define _align __attribute__((aligned(32))) /* AVX packed instructions need 32-byte alignment */
#elif defined(__GNUC__)
#  define _align __attribute__((aligned(16))) /* SSE instructions need 16-byte alignment */
#else
#  define _align
#endif

#endif
