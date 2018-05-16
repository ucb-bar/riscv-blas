/* drotm.f -- translated by f2c (version 20160102).
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

/* > \brief \b DROTM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DROTM(N,DX,INCX,DY,INCY,DPARAM) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION DPARAM(5),DX(*),DY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX */
/* > */
/* >    (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN */
/* >    (DY**T) */
/* > */
/* >    DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE */
/* >    LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY. */
/* >    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */
/* > */
/* >    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0 */
/* > */
/* >      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0) */
/* >    H=(          )    (          )    (          )    (          ) */
/* >      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0). */
/* >    SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
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
/* > */
/* > \param[in,out] DY */
/* > \verbatim */
/* >          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >         storage spacing between elements of DY */
/* > \endverbatim */
/* > */
/* > \param[in] DPARAM */
/* > \verbatim */
/* >          DPARAM is DOUBLE PRECISION array, dimension (5) */
/* >     DPARAM(1)=DFLAG */
/* >     DPARAM(2)=DH11 */
/* >     DPARAM(3)=DH21 */
/* >     DPARAM(4)=DH12 */
/* >     DPARAM(5)=DH22 */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup double_blas_level1 */

/*  ===================================================================== */
/* Subroutine */ int drotm_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *dparam)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal two = 2.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal w, z__;
    static integer kx, ky;
    static doublereal dh11, dh12, dh21, dh22, dflag;
    static integer nsteps;


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
/*     .. Data statements .. */
    /* Parameter adjustments */
    --dparam;
    --dy;
    --dx;

    /* Function Body */
/*     .. */

    dflag = dparam[1];
    if (*n <= 0 || dflag + two == zero) {
	return 0;
    }
    if (*incx == *incy && *incx > 0) {

	nsteps = *n * *incx;
	if (dflag < zero) {
	    dh11 = dparam[2];
	    dh12 = dparam[4];
	    dh21 = dparam[3];
	    dh22 = dparam[5];
	    i__1 = nsteps;
	    i__2 = *incx;
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
		w = dx[i__];
		z__ = dy[i__];
		dx[i__] = w * dh11 + z__ * dh12;
		dy[i__] = w * dh21 + z__ * dh22;
	    }
	} else if (dflag == zero) {
	    dh12 = dparam[4];
	    dh21 = dparam[3];
	    i__2 = nsteps;
	    i__1 = *incx;
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
		w = dx[i__];
		z__ = dy[i__];
		dx[i__] = w + z__ * dh12;
		dy[i__] = w * dh21 + z__;
	    }
	} else {
	    dh11 = dparam[2];
	    dh22 = dparam[5];
	    i__1 = nsteps;
	    i__2 = *incx;
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
		w = dx[i__];
		z__ = dy[i__];
		dx[i__] = w * dh11 + z__;
		dy[i__] = -w + dh22 * z__;
	    }
	}
    } else {
	kx = 1;
	ky = 1;
	if (*incx < 0) {
	    kx = (1 - *n) * *incx + 1;
	}
	if (*incy < 0) {
	    ky = (1 - *n) * *incy + 1;
	}

	if (dflag < zero) {
	    dh11 = dparam[2];
	    dh12 = dparam[4];
	    dh21 = dparam[3];
	    dh22 = dparam[5];
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		w = dx[kx];
		z__ = dy[ky];
		dx[kx] = w * dh11 + z__ * dh12;
		dy[ky] = w * dh21 + z__ * dh22;
		kx += *incx;
		ky += *incy;
	    }
	} else if (dflag == zero) {
	    dh12 = dparam[4];
	    dh21 = dparam[3];
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		w = dx[kx];
		z__ = dy[ky];
		dx[kx] = w + z__ * dh12;
		dy[ky] = w * dh21 + z__;
		kx += *incx;
		ky += *incy;
	    }
	} else {
	    dh11 = dparam[2];
	    dh22 = dparam[5];
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		w = dx[kx];
		z__ = dy[ky];
		dx[kx] = w * dh11 + z__;
		dy[ky] = -w + dh22 * z__;
		kx += *incx;
		ky += *incy;
	    }
	}
    }
    return 0;
} /* drotm_ */

