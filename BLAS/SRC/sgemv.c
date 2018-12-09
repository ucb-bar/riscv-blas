/* sgemv.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"
#include "custom-utils.h"

/* > \brief \b SGEMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER INCX,INCY,LDA,M,N */
/*       CHARACTER TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEMV  performs one of the matrix-vector operations */
/* > */
/* >    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are vectors and A is an */
/* > m by n matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry, TRANS specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */
/* > */
/* >              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y. */
/* > */
/* >              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry, M specifies the number of rows of the matrix A. */
/* >           M must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is REAL */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension ( LDA, N ) */
/* >           Before entry, the leading m by n part of the array A must */
/* >           contain the matrix of coefficients. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, m ). */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is REAL array, dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/* >           and at least */
/* >           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
/* >           Before entry, the incremented array X must contain the */
/* >           vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is REAL */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is REAL array, dimension at least */
/* >           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/* >           and at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
/* >           Before entry with BETA non-zero, the incremented array Y */
/* >           must contain the vector y. On exit, Y is overwritten by the */
/* >           updated vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           Y. INCY must not be zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup single_blas_level2 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 2 Blas routine. */
/* >  The vector and matrix arguments are not referenced when N = 0, or M = 0 */
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgemv_(char *trans, integer *m, integer *n, real *alpha, 
	real *a, integer *lda, real *x, integer *incx, real *beta, real *y, 
	integer *incy, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    printf("sgemv %d %d %d %d %d %c\n", *m, *n, *lda, *incx, *incy, *trans);
    /* Local variables */
    static integer i__, j, ix, iy, jx, jy, kx, ky, info;
    static real temp;
    static integer lenx, leny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- Reference BLAS level2 routine (version 3.7.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)
	    ) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("SGEMV ", &info, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0.f && *beta == 1.f) {
	return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

    if (*beta != 1.f) {

        hwacha_init();
        setvcfg(0, 1, 0, 1);
        int vl = 0;
	float* cy = y + ky;
        void* pre = PRELOAD("blas2");
        i__ = 0;
        i__1 = leny;
	if (*incy == 1) {
	    if (*beta == 0.f) {
                while (i__1 - i__ > 0) {
                  vl = setvlen(i__1 - i__);
                  MEMTOUCH(cy, float, vl - 1);
                  asm volatile ("vmca va0, %0" : : "r" (cy));
                  VF("sgemv_zero_loop");
                  cy += vl;
                  i__ += vl;
                }

	    } else {
                asm volatile ("vmcs vs1, %0" : : "r" (*beta));
                while (i__1 - i__ > 0) {
                  vl = setvlen(i__1 - i__);
                  MEMTOUCH(cy, float, vl - 1);
                  asm volatile ("vmca va0, %0" : : "r" (cy));
                  VF("sgemv_beta_loop");
                  cy += vl;
                  i__ += vl;
                }
	    }
	} else {
	    if (*beta == 0.f) {
                asm volatile ("vmca va1, %0" : : "r" (*incy << 2));
                while (i__1 - i__ > 0) {
                  vl = setvlen(max(i__1 - i__,16));
                  MEMTOUCH(cy, float, (vl-1) * *incy);
                  asm volatile ("vmca va0, %0" : : "r" (cy));
                  VF("sgemv_stride_zero_loop");
                  cy += (*incy * vl);
                  i__ += vl;
                }
	    } else {
                asm volatile ("vmcs vs1, %0" : : "r" (*beta));
                asm volatile ("vmca va1, %0" : : "r" (*incy << 2));
                while (i__1 - i__  > 0) {
                  vl = setvlen(max(i__1 - i__,16));
                  MEMTOUCH(cy, float, (vl-1) * *incy);
                  asm volatile ("vmca va0, %0" : : "r" (cy));
                  VF("sgemv_stride_beta_loop");
                  cy += (*incy * vl);
                  i__ += vl;
                }
	    }
	}
    }
    if (*alpha == 0.f) {
      asm volatile("fence");
	return 0;
    }
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  y := alpha*A*x + y. */

        hwacha_init();
        setvcfg(0, 3, 0, 1);
        int vl = 0;
        void* pre = PRELOAD("blas2");
	float* cy = y + ky;
        float* cx = x + kx;
        i__ = 0;
	i__1 = *n;
        i__2 = *m;
	if (*incy == 1) {
            while (i__2 - i__ > 0) {
              vl = setvlen(i__2 - i__);
              jx = kx;
              VF("sgemv_start");
              for (j = 1; j <= i__1; j++) {
                MEMTOUCH(a + i__ + 1 + j * a_dim1, float, vl-1);
                asm volatile ("vmca va0, %0" : : "r" (a + i__ + 1 + j * a_dim1));
                asm volatile ("vmcs vs1, %0" : : "r" (*alpha * x[jx]));
                VF("sgemv_loop");
                jx += *incx;
              }
              MEMTOUCH(cy + i__, float, vl-1);
              asm volatile ("vmca va0, %0" : : "r" (cy + i__));
              VF("sgemv_stop");
              i__ += vl;
            }
	} else {
            asm volatile ("vmca va1, %0" : : "r" (*incy << 2));
            while (i__2 - i__ > 0) {
              vl = setvlen(max(i__2 - i__,16));
              jx = kx;
              VF("sgemv_start");
              for (j = 1; j <= i__1; j++) {
                MEMTOUCH(a + i__ + 1 + j * a_dim1, float, vl-1);
                asm volatile ("vmca va0, %0" : : "r" (a + i__ + 1 + j * a_dim1));
                asm volatile ("vmcs vs1, %0" : : "r" (*alpha * x[jx]));
                VF("sgemv_loop");
                jx += *incx;
              }
              MEMTOUCH(cy + i__ * (*incy), float, (vl-1) * (*incy));
              asm volatile ("vmca va0, %0" : : "r" (cy + i__ * (*incy)));
              VF("sgemv_stride_stop");
              i__ += vl;
            }
	}
    } else {

/*        Form  y := alpha*A**T*x + y. */
        hwacha_init();
        setvcfg(0, 3, 0, 1);
        int vl = 0;
        void* pre = PRELOAD("blas2");
	float* cy = y + ky;
        float* cx = x + kx;
        i__ = 0;
	i__1 = *m;
        i__2 = *n;
	if (*incx == 1) {
            while (i__2 - i__ > 0) {
              vl = setvlen(max(i__2 - i__, 16));
              jx = kx;
              VF("sgemv_start");
              asm volatile ("vmca va1, %0" : : "r" (*incy << 2));
              asm volatile ("vmca va2, %0" : : "r" (a_dim1 << 2));
              for (j = 1; j <= i__1; j++) {
                MEMTOUCH(a + j + (i__ + 1) * a_dim1, float, (vl-1) * a_dim1);
                asm volatile ("vmca va0, %0" : : "r" (a + j + (i__+1) * a_dim1));
                asm volatile ("vmcs vs1, %0" : : "r" (*alpha * x[jx]));
                VF("sgemv_tranpose_loop");
                jx += 1;
              }
              MEMTOUCH(cy + i__ * (*incy), float, (vl-1) * (*incy));
              asm volatile ("vmca va0, %0" : : "r" (cy + i__ * (*incy)));
              VF("sgemv_stride_stop");
              i__ += vl;
            }

	} else {
            while (i__2 - i__ > 0) {
              vl = setvlen(max(i__2 - i__,16));
              jx = kx;
              VF("sgemv_start");
              asm volatile ("vmca va1, %0" : : "r" (*incy << 2));
              asm volatile ("vmca va2, %0" : : "r" (a_dim1 << 2));
              for (j = 1; j <= i__1; j++) {
                MEMTOUCH(a + j + (i__+1)*a_dim1, float, (vl-1) * (a_dim1));
                asm volatile ("vmca va0, %0" : : "r" (a + j + (i__+1) * a_dim1));
                asm volatile ("vmcs vs1, %0" : : "r" (*alpha * x[jx]));
                VF("sgemv_tranpose_loop");
                jx += *incx;
              }
              MEMTOUCH(cy + i__ * (*incy), float, (vl-1) * (*incy));
              asm volatile ("vmca va0, %0" : : "r" (cy + i__ * (*incy)));
              VF("sgemv_stride_stop");
              i__ += vl;
            }
	}
    }
    asm volatile("fence");
    return 0;

/*     End of SGEMV . */

} /* sgemv_ */

