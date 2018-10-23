/* dscal.f -- translated by f2c (version 20160102).
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
/* > \brief \b DSCAL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSCAL(N,DA,DX,INCX) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION DA */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION DX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DSCAL scales a vector by a constant. */
/* >    uses unrolled loops for increment equal to 1. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] DA */
/* > \verbatim */
/* >          DA is DOUBLE PRECISION */
/* >           On entry, DA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DX */
/* > \verbatim */
/* >          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of DX */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup double_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, linpack, 3/11/78. */
/* >     modified 3/93 to return if incx .le. 0. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dscal_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, m, mp1, nincx;


/*  -- Reference BLAS level1 routine (version 3.8.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    hwacha_init();
    setvcfg(1, 0, 0, 1);
    int vl = 0;
    double* dxa = dx + 1;
    void* pre = PRELOAD("blas1");
    asm volatile ("vmcs vs1, %0" : : "r" (*da));
    if (*incx == 1) {

/*        code for increment equal to 1 */
      i__  = 0;
      i__1 = *n;
      while (i__1 - i__ > 0) {
        vl = setvlen(i__1 - i__);
        asm volatile ("vmca va0, %0" : : "r" (dxa));
        VF("dscal_unit");
        dxa += vl;
        i__ += vl;
      }
    } else {
      int ix = 1;
      if (*incx < 0) {
        ix = (-(*n) + 1) * *incx + 1;
      }
      dxa = dx + ix;

      i__  = 0;
      i__1 = *n;
      asm volatile ("vmca va1, %0" : : "r" (*incx << 3));
      while (i__1 - i__ > 0) {
        vl = setvlen(i__1 - i__);
        asm volatile ("vmca va0, %0" : : "r" (dxa));
        VF("dscal_stride");
        dxa += vl * (*incx);
        i__ += vl;
      }
    }
    return 0;
} /* dscal_ */

