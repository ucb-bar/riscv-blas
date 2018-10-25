/* sasum.f -- translated by f2c (version 20160102).
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

/* > \brief \b SASUM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SASUM(N,SX,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SASUM takes the sum of the absolute values. */
/* >    uses unrolled loops for increment equal to one. */
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
/* >     modified 3/93 to return if incx .le. 0. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
real sasum_(integer *n, real *sx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

    /* Local variables */
    static integer i__, m, mp1, nincx;
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
    --sx;

    /* Function Body */
    ret_val = 0.f;
    stemp = 0.f;
    if (*n <= 0 || *incx <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
/*        code for increment equal to 1 */

        hwacha_init();
        setvcfg(0, 2, 0, 2);
        float* cx = sx + 1;
        void* pre = PRELOAD("blas1");
        int vl = setvlen(*n);
        VF("sasum_pre");
        i__ = 0;
        i__1 = *n;
        //multiply accumulate
        while (i__1 - i__ > 0) {
            vl = setvlen(i__1 - i__);
            asm volatile ("vmca va1, %0" : : "r" (cx));
            VF("sasum_loop");
            cx += vl;
            i__ += vl;
        }

        vl = setvlen(*n);
        int vl_pad = vl + vl % 2;
        float* ta = (float*)malloc(vl_pad * sizeof(float));
        ta[vl_pad - 1] = 0.f;
        asm volatile ("vmca va0, %0" : : "r" (ta));
        VF("sasum_post");

        float *ta2;
        i__1 = vl_pad >> 1;
        while (i__1 > 0) {
            vl = setvlen(i__1);
            ta2 = ta + vl;
            asm volatile ("vmca va1, %0" : : "r" (ta2));
            VF("sasum_reduce_loop");
            i__1 = vl >> 1;
        }

        ret_val = *ta;
        free(ta);
        stemp = ret_val;
    } else {

/*        code for increment not equal to 1 */

        hwacha_init();
        setvcfg(0, 2, 0, 2);
        float* cx = sx + 1;
	if (*incx < 0) {
	    cx = sx + (-(*n) + 1) * *incx + 1;
	}
        void* pre = PRELOAD("blas1");
        int vl = setvlen(*n);
        VF("sasum_pre");
        asm volatile ("vmca va2, %0" : : "r" (*incx << 2));
        i__ = 0;
        i__1 = *n;
        //multiply accumulate
        while (i__1 - i__ > 0) {
            vl = setvlen(i__1 - i__);
            asm volatile ("vmca va1, %0" : : "r" (cx));
            VF("sasum_stride_loop");
            cx += vl;
            i__ += vl;
        }

        vl = setvlen(*n);
        int vl_pad = vl + vl % 2;
        float* ta = (float*)malloc(vl_pad * sizeof(float));
        ta[vl_pad - 1] = 0.f;
        asm volatile ("vmca va0, %0" : : "r" (ta));
        VF("sasum_post");

        float *ta2;
        i__1 = vl_pad >> 1;
        while (i__1 > 0) {
            vl = setvlen(i__1);
            ta2 = ta + vl;
            asm volatile ("vmca va1, %0" : : "r" (ta2));
            VF("sasum_reduce_loop");
            i__1 = vl >> 1;
        }

        ret_val = *ta;
        free(ta);
    }
    return ret_val;
} /* sasum_ */

