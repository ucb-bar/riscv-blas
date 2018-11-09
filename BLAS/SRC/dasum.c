/* dasum.f -- translated by f2c (version 20160102).
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

/* > \brief \b DASUM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX) */

/*       .. Scalar Arguments .. */
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
/* >    DASUM takes the sum of the absolute values. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] DX */
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
doublereal dasum_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static integer i__, m, mp1;
    static doublereal dtemp;
    static integer nincx;


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
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0 || *incx <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
/*        code for increment equal to 1 */

        hwacha_init();
        setvcfg(2, 0, 0, 2);
        double* cx = dx + 1;
        void* pre = PRELOAD("blas1");
        int vl = setvlen(*n);
        VF("dasum_pre");
        i__ = 0;
        i__1 = *n;
        //multiply accumulate
        while (i__1 - i__ > 0) {
            vl = setvlen(i__1 - i__);
            asm volatile ("vmca va1, %0" : : "r" (cx));
            VF("dasum_loop");
            cx += vl;
            i__ += vl;
        }

        vl = setvlen(*n);
        int vl_pad = vl + vl % 2;
        double* ta = (double*)malloc(vl_pad * sizeof(double));
        ta[vl_pad - 1] = 0.f;
        asm volatile ("vmca va0, %0" : : "r" (ta));
        VF("dasum_post");

        double *ta2;
        i__1 = vl_pad >> 1;
        while (i__1 > 0) {
            vl = setvlen(i__1);
            ta2 = ta + vl;
            asm volatile ("vmca va1, %0" : : "r" (ta2));
            VF("dasum_reduce_loop");
            i__1 = vl >> 1;
        }
        asm volatile("fence");
        ret_val = *ta;
        free(ta);
    } else {

/*        code for increment not equal to 1 */

        hwacha_init();
        setvcfg(0, 2, 0, 2);
        double* cx = dx + 1;
	if (*incx < 0) {
	    cx = dx + (-(*n) + 1) * *incx + 1;
	}
        void* pre = PRELOAD("blas1");
        int vl = setvlen(*n);
        VF("dasum_pre");
        asm volatile ("vmca va2, %0" : : "r" (*incx << 3));
        i__ = 0;
        i__1 = *n;
        //multiply accumulate
        while (i__1 - i__ > 0) {
            vl = setvlen(i__1 - i__);
            asm volatile ("vmca va1, %0" : : "r" (cx));
            VF("dasum_stride_loop");
            cx += vl;
            i__ += vl;
        }

        vl = setvlen(*n);
        int vl_pad = vl + vl % 2;
        double* ta = (double*)malloc(vl_pad * sizeof(double));
        ta[vl_pad - 1] = 0.f;
        asm volatile ("vmca va0, %0" : : "r" (ta));
        VF("dasum_post");

        double *ta2;
        i__1 = vl_pad >> 1;
        while (i__1 > 0) {
            vl = setvlen(i__1);
            ta2 = ta + vl;
            asm volatile ("vmca va1, %0" : : "r" (ta2));
            VF("dasum_reduce_loop");
            i__1 = vl >> 1;
        }
        asm volatile("fence");
        ret_val = *ta;
        free(ta);
    }
    return ret_val;
} /* dasum_ */

