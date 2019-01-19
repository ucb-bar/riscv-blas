/* sdot.f -- translated by f2c (version 20160102).
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

/* > \brief \b SDOT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SDOT(N,SX,INCX,SY,INCY) */

/*       .. Scalar Arguments .. */
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
/* >    SDOT forms the dot product of two vectors. */
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
/* > \param[in] SY */
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
real sdot_(integer *n, real *sx, integer *incx, real *sy, integer *incy)
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;
    static real stemp;


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
    stemp = 0.f;
    ret_val = 0.f;
    if (*n <= 0) {
	return ret_val;
    }
    hwacha_init();
    setvcfg(0, 3, 0, 1);
    PRELOAD("blas1");
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */


/*        clean-up loop */

	m = *n % 32;
	if (m != 0) {
	    i__1 = m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		stemp += sx[i__] * sy[i__];
	    }
	    if (*n < 32) {
		ret_val = stemp;
		return ret_val;
	    }
	}
        float buffer[32];
        memset(buffer, 0, 32 * sizeof(float));
        setvlen(32);
        VF("sdot_pre");
	mp1 = m + 1;
	i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 32) {
          MEMTOUCH(sx + i__, float, 32);
          MEMTOUCH(sy + i__, float, 32);
          asm volatile("vmca va0, %0" : : "r" (sx + i__));
          asm volatile("vmca va1, %0" : : "r" (sy + i__));
          VF("sdot_loop");
        }
        asm volatile ("vmca va2, %0" : : "r" (buffer));
        VF("sdot_post");
        asm volatile ("fence");
        //printf("stemp %.3f\n", stemp);
        for (i__ = 0; i__ < 32; i__++) {
          stemp += buffer[i__];
        }
        //printf("stemp %.3f\n", stemp);
        
	/* for (i__ = mp1; i__ <= i__1; i__ += 5) { */
	/*     stemp = stemp + sx[i__] * sy[i__] + sx[i__ + 1] * sy[i__ + 1] +  */
	/* 	    sx[i__ + 2] * sy[i__ + 2] + sx[i__ + 3] * sy[i__ + 3] +  */
	/* 	    sx[i__ + 4] * sy[i__ + 4]; */
	/* } */
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

	ix = 1;
	iy = 1;
	if (*incx < 0) {
	    ix = (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
	    iy = (-(*n) + 1) * *incy + 1;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    stemp += sx[ix] * sy[iy];
	    ix += *incx;
	    iy += *incy;
	}
    }
    ret_val = stemp;
    return ret_val;
} /* sdot_ */

