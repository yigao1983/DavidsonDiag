#include <iostream>
#include <cassert>
#include "lapack.hpp"

namespace lapack {
  
  void dgemm(int mm, int nn, int kk, char transa, char transb, dreal alpha, const dreal *mata,
             const dreal *matb, dreal beta, dreal *matc)
  {
    int lda, ldb, ldc;
    
    if(transa == 'N' || transa == 'n') {
      lda = (1 > mm) ? 1 : mm;
    } else {
      lda = (1 > kk) ? 1 : kk;
    }
    
    if(transb == 'N' || transb == 'n') {
      ldb = (1 > kk) ? 1 : kk;
    } else {
      ldb = (1 > nn) ? 1 : nn;
    }
    
    ldc = (1 > mm) ? 1 : mm;
    
    dgemm_(&transa, &transb, &mm, &nn, &kk, &alpha, mata, &lda, matb, &ldb, &beta, matc, &ldc); 
  }
  
  void zgemm(int mm, int nn, int kk, char transa, char transb, dcmplx alpha, const dcmplx *mata,
             const dcmplx *matb, dcmplx beta, dcmplx *matc)
  {
    int lda, ldb, ldc;
    
    if(transa == 'N' || transa == 'n') {
      lda = (1 > mm) ? 1 : mm;
    } else {
      lda = (1 > kk) ? 1 : kk;
    }
    
    if(transb == 'N' || transb == 'n') {
      ldb = (1 > kk) ? 1 : kk;
    } else {
      ldb = (1 > nn) ? 1 : nn;
    }
    
    ldc = (1 > mm) ? 1 : mm;
    
    zgemm_(&transa, &transb, &mm, &nn, &kk, &alpha, mata, &lda, matb, &ldb, &beta, matc, &ldc); 
  }
  
  void dsyev(int nn, dreal *AA, dreal *ww)
  {
    int   lda, lwork, info;
    char  jobz = 'V', uplo = 'U';
    dreal *work = NULL;
    
    lda = (1 > nn) ? 1 : nn;
    lwork = (1 > 3*nn-1) ? 1 : 3*nn-1;
    work = new dreal[lwork];
    
    dsyev_(&jobz, &uplo, &nn, AA, &lda, ww, work, &lwork, &info);
    
    delete [] work;
    
    assert(info == 0);
  }
  
  void dsyevd(int nn, dreal *AA, dreal *ww)
  {
    int   lda, lwork, liwork, info;
    char  jobz = 'V', uplo = 'U';
    int   *iwork = NULL;
    dreal *work = NULL;
    
    lda = (1 > nn) ? 1 : nn;
    lwork = (1 > 1 + 6*nn + 2*nn*nn) ? 1 : 1 + 6*nn + 2*nn*nn;
    liwork = (1 > 3 + 5*nn) ? 1: 3 + 5*nn;
    work = new dreal[lwork];
    iwork = new int[liwork];
    
    dsyevd_(&jobz, &uplo, &nn, AA, &lda, ww, work, &lwork, iwork, &liwork, &info);
    
    delete [] iwork;
    delete [] work;
    
    assert(info == 0);
  }
  
  void dsyevx(int nn, dreal vl, dreal vu, int il, int iu, dreal *AA, dreal *ww, dreal *ZZ)
  {
    const dreal abstol = 1e-8;
    int   lda, ldz, mm, lwork, info;
    char  jobz = 'V', range = 'I', uplo = 'U';
    int   *iwork = NULL, *ifail = NULL;
    dreal *work = NULL;
    
    lda = (1 > nn) ? 1 : nn;
    ldz = (1 > nn) ? 1 : nn;
    lwork = (1 > 8*nn) ? 1 : 8*nn;
    work = new dreal[lwork];
    iwork = new int[(1 > 5*nn) ? 1 : 5*nn];
    ifail = new int[(1 > nn) ? 1 : nn];
    
    dsyevx_(&jobz, &range, &uplo, &nn, AA, &lda, &vl, &vu, &il, &iu, &abstol, &mm, ww, ZZ, &ldz,
            work, &lwork, iwork, ifail, &info);
    
    delete [] ifail;
    delete [] iwork;
    delete [] work;
    
    assert(info == 0);
  }
  
  void dsyevr(int nn, dreal vl, dreal vu, int il, int iu, dreal *AA, dreal *ww, dreal *ZZ)
  {
    const dreal abstol = 1e-8;
    int   lda, ldz, mm, lwork, liwork, info;
    char  jobz = 'V', range = 'I', uplo = 'U';
    int   *iwork = NULL, *isuppz = NULL;
    dreal *work = NULL;
    
    lda = (1 > nn) ? 1 : nn;
    ldz = (1 > nn) ? 1 : nn;
    lwork = (1 > 26*nn) ? 1 : 26*nn;
    liwork = (1 > 10*nn) ? 1 : 10*nn;
    work = new dreal[lwork];
    iwork = new int[liwork];
    isuppz = new int[(2 > 2*nn) ? 2 : 2*nn];
    
    dsyevr_(&jobz, &range, &uplo, &nn, AA, &lda, &vl, &vu, &il, &iu, &abstol, &mm, ww, ZZ, &ldz,
            isuppz, work, &lwork, iwork, &liwork, &info);
    
    delete [] isuppz;
    delete [] iwork;
    delete [] work;
    
    assert(info == 0);
  }
  
  void zheev(int nn, dcmplx *AA, dreal *ww)
  {
    int    lda, lwork, info;
    char   jobz = 'V', uplo = 'U';
    dreal  *rwork = NULL;
    dcmplx *work = NULL;
    
    lda = (1 > nn) ? 1 : nn;
    lwork = (1 > 2*nn-1) ? 1 : 2*nn-1;
    work = new dcmplx[lwork];
    rwork = new dreal[(1 > 3*nn-2) ? 1 : 3*nn-2];
    
    zheev_(&jobz, &uplo, &nn, AA, &lda, ww, work, &lwork, rwork, &info);
    
    delete [] rwork;
    delete [] work;
    
    assert(info == 0);
  }
  
  void zheevd(int nn, dcmplx *AA, dreal *ww)
  {
    int    lda, lwork, lrwork, liwork, info;
    char   jobz = 'V', uplo = 'U';
    int    *iwork = NULL;
    dreal  *rwork = NULL;
    dcmplx *work = NULL;
    
    lda = (1 > nn) ? 1 : nn;
    lwork = (1 > 2*nn + nn*nn) ? 1 : 2*nn + nn*nn;
    lrwork = (1 > 1 + 5*nn + 2*nn*nn) ? 1 : 1 + 5*nn + 2*nn*nn;
    liwork = (1 > 3 + 5*nn) ? 1 : 3 + 5*nn;
    work = new dcmplx[lwork];
    rwork = new dreal[lrwork];
    iwork = new int[liwork];
    
    zheevd_(&jobz, &uplo, &nn, AA, &lda, ww, work, &lwork, rwork, &lrwork,
            iwork, &liwork, &info);
    
    delete [] iwork;
    delete [] rwork;
    delete [] work;
    
    assert(info == 0);
  }
  
  void zheevx(int nn, dreal vl, dreal vr, int il, int iu, dcmplx *AA, dreal *ww, dcmplx *ZZ)
  {
    const dreal abstol = 1e-8;
    int    lda, ldz, mm, lwork, info;
    char   jobz = 'V', range = 'I', uplo = 'U';
    int    *iwork = NULL, *ifail = NULL;
    dreal  *rwork = NULL;
    dcmplx *work = NULL;
    
    lda = (1 > nn) ? 1 : nn;
    ldz = (1 > nn) ? 1 : nn;
    lwork = (1 > 2*nn) ? 1 : 2*nn;
    work = new dcmplx[lwork];
    rwork = new dreal[(1 > 7*nn) ? 1 : 7*nn];
    iwork = new int[(1 > 5*nn) ? 1 : 5*nn];
    ifail = new int[(1 > nn) ? 1 : nn];
    
    zheevx_(&jobz, &range, &uplo, &nn, AA, &lda, &vl, &vr, &il, &iu, &abstol, &mm, ww, ZZ, &ldz,
            work, &lwork, rwork, iwork, ifail, &info);
    
    delete [] ifail;
    delete [] iwork;
    delete [] rwork;
    delete [] work;
    
    assert(info == 0);
  }
  
  void zheevr(int nn, dreal vl, dreal vr, int il, int iu, dcmplx *AA, dreal *ww, dcmplx *ZZ)
  {
    const dreal abstol = 1e-8;
    int    lda, ldz, mm, lwork, lrwork, liwork, info;
    char   jobz = 'V', range = 'I', uplo = 'U';
    int    *iwork = NULL, *isuppz = NULL;
    dreal  *rwork = NULL;
    dcmplx *work = NULL;
    
    lda = (1 > nn) ? 1 : nn;
    ldz = (1 > nn) ? 1 : nn;
    lwork = (1 > 2*nn) ? 1 : 2*nn;
    lrwork = (1 > 24*nn) ? 1 : 24*nn;
    liwork = (1 > 10*nn) ? 1 : 10*nn;
    work = new dcmplx[lwork];
    rwork = new dreal[lrwork];
    iwork = new int[liwork];
    isuppz = new int[(2 > 2*nn) ? 2 : 2*nn];
    
    zheevr_(&jobz, &range, &uplo, &nn, AA, &lda, &vl, &vr, &il, &iu, &abstol, &mm, ww, ZZ, &ldz,
            isuppz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    
    delete [] isuppz;
    delete [] iwork;
    delete [] rwork;
    delete [] work;
    
    assert(info == 0);
  }
  
  void dgetrf(int mm, int nn, dreal *AA, int *ipiv)
  {
    int lda, info;
    
    lda = (1 > mm) ? 1 : mm;
    
    dgetrf_(&mm, &nn, AA, &lda, ipiv, &info);
    
    assert(info == 0);
  }
  
  void dgetri(int nn, dreal *AA, int *ipiv)
  {
    int   lda, lwork, info;
    dreal *work = NULL;
    
    lda = (1 > nn) ? 1 : nn;
    lwork = lda;
    
    work = new dreal[lwork];
    
    dgetri_(&nn, AA, &lda, ipiv, work, &lwork, &info);
    
    delete [] work;
    
    assert(info == 0);
  }
  
  void dinvs(int nn, dreal *AA)
  {
    int *ipiv = NULL;
    
    ipiv = new int[nn];
    
    dgetrf(nn, nn, AA, ipiv);
    dgetri(nn, AA, ipiv);
    
    delete [] ipiv;
  }
  
  void zgetrf(int mm, int nn, dcmplx *AA, int *ipiv)
  {
    int lda, info;
    
    lda = (1 > mm) ? 1 : mm;
    
    zgetrf_(&mm, &nn, AA, &lda, ipiv, &info);
    
    assert(info == 0);
  }
  
  void zgetri(int nn, dcmplx *AA, int *ipiv)
  {
    int    lda, lwork, info;
    dcmplx *work = NULL;
    
    lda = (1 > nn) ? 1 : nn;
    lwork = lda;
    
    work = new dcmplx[lwork];
    
    zgetri_(&nn, AA, &lda, ipiv, work, &lwork, &info);
    
    delete [] work;
    
    assert(info == 0);
  }
  
  void zinvs(int nn, dcmplx *AA)
  {
    int *ipiv = NULL;
    
    ipiv = new int[nn];
    
    zgetrf(nn, nn, AA, ipiv);
    zgetri(nn, AA, ipiv);
    
    delete [] ipiv;
  }
}
