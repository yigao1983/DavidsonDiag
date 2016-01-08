#ifndef LAPACK_HPP
#define LAPACK_HPP

#include "precision.hpp"

extern "C" void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N, const int *K,
                       const dreal *ALPHA, const dreal *A, const int *LDA,
                       const dreal *B, const int *LDB, const dreal *BETA, dreal *C, const int *LDC); 

extern "C" void zgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N, const int *K,
                       const dcmplx *ALPHA, const dcmplx *A, const int *LDA,
                       const dcmplx *B, const int *LDB, const dcmplx *BETA, dcmplx *C, const int *LDC);

extern "C" void dsyev_(const char *JOBZ, const char *UPLO, const int *N, dreal *A, const int *LDA,
                       dreal *W, dreal *WORK, const int *LWORK, int *INFO);

extern "C" void dsyevd_(const char *JOBZ, const char *UPLO, const int *N, dreal *A, const int *LDA,
                        dreal *W, dreal *WORK, const int *LWORK, int *IWORK, const int *LIWORK,
                        int *INFO);

extern "C" void dsyevx_(const char *JOBZ, const char *RANGE, const char *UPLO,
                        const int *N, dreal *A, const int *LDA, const dreal *VL, const dreal *VU,
                        const int *IL, const int *IU, const dreal *ABSTOL, int *M, dreal *W,
                        dreal *Z, const int *LDZ, dreal *WORK, const int *LWORK,
                        int *IWORK, int *IFAIL, int *INFO);

extern "C" void dsyevr_(const char *JOBZ, const char *RANGE, const char *UPLO,
                        const int *N, dreal *A, const int *LDA, const dreal *VL, const dreal *VU,
                        const int *IL, const int *IU, const dreal *ABSTOL, int *M, dreal *W,
                        dreal *Z, const int *LDZ, int *ISUPPZ, dreal *WORK, const int *LWORK,
                        int *IWORK, const int *LIWORK, int *INFO);

extern "C" void zheev_(const char *JOBZ, const char *UPLO, const int *N, dcmplx *A, const int *LDA,
                       dreal *W, dcmplx *WORK, const int *LWORK, dreal *RWORK, int *INFO);

extern "C" void zheevd_(const char *JOBZ, const char *UPLO, const int *N, dcmplx *A, const int *LDA,
                        dreal *W, dcmplx *WORK, const int *LWORK, dreal *RWORK, const int *LRWORK,
                        int *IWORK, const int *LIWORK, int *INFO);

extern "C" void zheevx_(const char *JOBZ, const char *RANGE, const char *UPLO,
                        const int *N, dcmplx *A, const int *LDA, const dreal *VL, const dreal *VU,
                        const int *IL, const int *IU, const dreal *ABSTOL, int *M, dreal *W,
                        dcmplx *Z, const int *LDZ, dcmplx *WORK, const int *LWORK,
                        dreal *RWORK, int *IWORK, int *IFAIL, int *INFO);

extern "C" void zheevr_(const char *JOBZ, const char *RANGE, const char *UPLO,
                        const int *N, dcmplx *A, const int *LDA, const dreal *VL, const dreal *VU,
                        const int *IL, const int *IU, const dreal *ABSTOL, int *M, dreal *W,
                        dcmplx *Z, const int *LDZ, int *ISUPPZ, dcmplx *WORK, const int *LWORK,
                        dreal *RWORK, const int *LRWORK, int *IWORK, const int *LIWORK, int *INFO);

extern "C" void dgetrf_(const int *M, const int *N, dreal *A, const int *LDA, int *IPIV, int *INFO);

extern "C" void dgetri_(const int *N, dreal *A, const int *LDA, const int *IPIV,
                        dreal *WORK, const int *LWORK, int *INFO);

extern "C" void zgetrf_(const int *M, const int *N, dcmplx *A, const int *LDA, int *IPIV, int *INFO);

extern "C" void zgetri_(const int *N, dcmplx *A, const int *LDA, const int *IPIV,
                        dcmplx *WORK, const int *LWORK, int *INFO);

namespace lapack {
  
  void dgemm(int mm, int nn, int kk, char transa, char transb, dreal alpha, const dreal *mata,
             const dreal *matb, dreal beta, dreal *matc);
  
  void zgemm(int mm, int nn, int kk, char transa, char transb, dcmplx alpha, const dcmplx *mata,
             const dcmplx *matb, dcmplx beta, dcmplx *matc);
  
  void dsyev(int nn, dreal *AA, dreal *ww);
  
  void dsyevd(int nn, dreal *AA, dreal *ww);
  
  void dsyevx(int nn, dreal vl, dreal vu, int il, int iu, dreal *AA, dreal *ww, dreal *ZZ);
  
  void dsyevr(int nn, dreal vl, dreal vu, int il, int iu, dreal *AA, dreal *ww, dreal *ZZ);
  
  void zheev(int nn, dcmplx *AA, dreal *ww);
  
  void zheevd(int nn, dcmplx *AA, dreal *ww);
  
  void zheevx(int nn, dreal vl, dreal vr, int il, int iu, dcmplx *AA, dreal *ww, dcmplx *ZZ);
  
  void zheevr(int nn, dreal vl, dreal vr, int il, int iu, dcmplx *AA, dreal *ww, dcmplx *ZZ);
  
  void dgetrf(int mm, int nn, dreal *AA, int *ipiv);
  
  void dgetri(int nn, dreal *AA, int *ipiv);
  
  void dinvs(int nn, dreal *AA);
  
  void zgetrf(int mm, int nn, dcmplx *AA, int *ipiv);
  
  void zgetri(int nn, dcmplx *AA, int *ipiv);
  
  void zinvs(int nn, dcmplx *AA);
}

#endif
