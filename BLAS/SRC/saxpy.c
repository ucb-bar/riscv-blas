/* saxpy.f -- translated by f2c (version 20160102).
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

/* > \brief \b SAXPY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY) */

/*       .. Scalar Arguments .. */
/*       REAL SA */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SX(*),SY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SAXPY constant times a vector plus a vector. */
/* >    uses unrolled loops for increments equal to one. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] SA */
/* > \verbatim */
/* >          SA is REAL */
/* >           On entry, SA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] SX */
/* > \verbatim */
/* >          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of SX */
/* > \endverbatim */
/* > */
/* > \param[in,out] SY */
/* > \verbatim */
/* >          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >         storage spacing between elements of SY */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup single_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, linpack, 3/11/78. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int saxpy_(integer *n, real *sa, real *sx, integer *incx, 
	real *sy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;


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
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*sa == 0.f) {
	return 0;
    }


    hwacha_init();
    setvcfg(0, 2, 0, 1);
    int vl = 0;
    float* cx = sx + 1;
    float* cy = sy + 1;
    void* pre = PRELOAD("blas1");

    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

        asm volatile ("vmcs vs1, %0" : : "r" (*sa));
        i__ = 0;
        i__1 = *n;
        while (i__1 - i__ > 0) {
          vl = setvlen(i__1 - i__);
          asm volatile ("vmca va0, %0" : : "r" (cx));
          asm volatile ("vmca va1, %0" : : "r" (cy));
          VF("saxpy_loop");
          cx += vl;
          cy += vl;
          i__ += vl;
        }
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

	ix = 1;
	iy = 1;
	if (*incx < 0) {
	    cx = sx + (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
	    cy = sy + (-(*n) + 1) * *incy + 1;
	}

        asm volatile ("vmcs vs1, %0" : : "r" (*sa));
        asm volatile ("vmca va2, %0" : : "r" (*incx << 2));
        asm volatile ("vmca va3, %0" : : "r" (*incy << 2));
        i__ = 0;
        i__1 = *n;
        while (i__1 - i__ > 0) {
          vl = setvlen(i__1 - i__);
          asm volatile ("vmca va0, %0" : : "r" (cx));
          asm volatile ("vmca va1, %0" : : "r" (cy));
          VF("saxpy_stride_loop");
          cx += (*incx * vl);
          cy += (*incy * vl);
          i__ += vl;
        }
    }
    return 0;
} /* saxpy_ */

