/* ddot.f -- translated by f2c (version 20160102).
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

/* > \brief \b DDOT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION DX(*),DY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DDOT forms the dot product of two vectors. */
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
/* > */
/* > \param[in] DY */
/* > \verbatim */
/* >          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >         storage spacing between elements of DY */
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
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;
    static doublereal dtemp;


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
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */
        hwacha_init();
        setvcfg(3, 0, 0, 1);
        int vl = 0;
        double* cy = dy+1;
        double* cx = dx+1;
        void* pre = PRELOAD("blas1");
        vl = setvlen(*n);
        VF("ddot_pre");
        i__ = 0;
        i__1 = *n;
        //multiply accumulate
        while (i__1 - i__ > 0) {
            vl = setvlen(i__1 - i__);
            asm volatile ("vmca va0, %0" : : "r" (cx));
            asm volatile ("vmca va1, %0" : : "r" (cy));
            VF("ddot_loop");
            cx += vl;
            cy += vl;
            i__ += vl;
        }

        //reduce 
        vl = setvlen(*n);
        int vl_pad = vl + vl % 2;
        double* ta = (double*)malloc(vl_pad * sizeof(double));
        ta[vl_pad - 1] = 0.f;
        asm volatile ("vmca va2, %0" : : "r" (ta));
        VF("ddot_post");

        double *ta2;
        i__1 = vl_pad >> 1;
        while (i__1 > 0) {
            vl = setvlen(i__1);
            ta2 = ta + vl;
            asm volatile ("vmca va1, %0" : : "r" (ta2));
            VF("ddot_reduce_loop");
            i__1 = vl >> 1;
        }

        ret_val = *ta;
        free(ta);
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */
        hwacha_init();
        setvcfg(3, 0, 0, 1);
        int vl = 0;
        double* cy = dy+1;
        double* cx = dx+1;
	if (*incx < 0) {
	    cx = dx + (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
	    cy = dy + (-(*n) + 1) * *incy + 1;
	}
        void* pre = PRELOAD("blas1");
        vl = setvlen(*n);
        VF("ddot_pre");
        asm volatile ("vmca va3, %0" : : "r" (*incx << 3));
        asm volatile ("vmca va4, %0" : : "r" (*incy << 3));
        i__ = 0;
        i__1 = *n;
        //multiply accumulate
        while (i__1 - i__ > 0) {
            vl = setvlen(i__1 - i__);
            asm volatile ("vmca va0, %0" : : "r" (cx));
            asm volatile ("vmca va1, %0" : : "r" (cy));
            VF("ddot_stride_loop");
            cx +=  vl;
            cy +=  vl;
            i__ += vl;
        }

        //reduce 
        vl = setvlen(*n);
        int vl_pad = vl + vl % 2;
        double* ta = (double*)malloc(vl_pad * sizeof(double));
        ta[vl_pad - 1] = 0.f;
        asm volatile ("vmca va2, %0" : : "r" (ta));
        VF("ddot_post");

        double *ta2;
        i__1 = vl_pad >> 1;
        while (i__1 > 0) {
            vl = setvlen(i__1);
            ta2 = ta + vl;
            asm volatile ("vmca va1, %0" : : "r" (ta2));
            VF("ddot_reduce_loop");
            i__1 = vl >> 1;
        }

        ret_val = *ta;
        free(ta);
    }
    return ret_val;
} /* ddot_ */

