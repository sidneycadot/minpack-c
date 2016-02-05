
// test_qrfac.c

#include <math.h>
#include <stdio.h>

#include "minpack_c.h"

#define nr 3
#define nc 7

int main()
{
  int m = 3; // number of rows
  int n = 7; // number of columns
  int lda = 5; // leading-dimension increment (must be >= m)

  double a[n * lda];

  printf("\n*** test_qrfac ***\n\n");

  for (int i = 0; i < m; ++i)
    {
      for (int j = 0; j < n; ++j)
	{
	  const int z = 10 * (i + 1) + (j + 1);
	  a[j * lda + i] = 1/sin(z);
	}
    }

  printf("'a' prior to QR factorization:\n\n");

  for (int i = 0; i < m; ++i)
    {
      for (int j = 0; j < n; ++j)
	{
	  printf("%15f", a[j * lda + i]);
	}
      printf("\n");
    }

  printf("\n");

  int ipvt[100];
  int lipvt = 100;

  double rdiag[100];
  double acnorm[100];
  double wa[100];

  qrfac(m, n, a, lda, false, ipvt, lipvt, rdiag, acnorm, wa);

  printf("'a' after QR factorization:\n\n");

  for (int i = 0; i < m; ++i)
    {
      for (int j = 0; j < n; ++j)
	{
	  printf("%15f", a[j * lda + i]);
	}
      printf("\n");
    }

  printf("\n");

  printf("'rdiag' after QR factorization:\n\n");

  for (int i = 0; i < n; ++i)
    {
      printf("%15f", rdiag[i]);
    }

  printf("\n\n");

  printf("'acnorm' after QR factorization:\n\n");

  for (int i = 0; i < n; ++i)
    {
      printf("%15f", acnorm[i]);
    }

  printf("\n\n");

  return 0;
}

