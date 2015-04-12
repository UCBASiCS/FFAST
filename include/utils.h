#ifndef UTILS_H
#define UTILS_H

#include <complex>

typedef double ffast_real;
typedef std::complex<double> ffast_complex;

// Perform a mod operation and make sure that the result is positive
inline int positiveMod(int a, int b)
{
  const int result = a % b;

  return result >= 0 ? result : result + b;
}

#endif // UTILS_H
