/* zdrot.f -- translated by f2c (version 20160102).
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

/* > \brief \b ZDROT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZDROT( N, CX, INCX, CY, INCY, C, S ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, INCY, N */
/*       DOUBLE PRECISION   C, S */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         CX( * ), CY( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Applies a plane rotation, where the cos and sin (c and s) are real */
/* > and the vectors cx and cy are complex. */
/* > jack dongarra, linpack, 3/11/78. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the vectors cx and cy. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CX */
/* > \verbatim */
/* >          CX is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( N - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array CX must contain the n */
/* >           element vector cx. On exit, CX is overwritten by the updated */
/* >           vector cx. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           CX. INCX must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CY */
/* > \verbatim */
/* >          CY is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( N - 1 )*abs( INCY ) ). */
/* >           Before entry, the incremented array CY must contain the n */
/* >           element vector cy. On exit, CY is overwritten by the updated */
/* >           vector cy. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           CY. INCY must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION */
/* >           On entry, C specifies the cosine, cos. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION */
/* >           On entry, S specifies the sine, sin. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16_blas_level1 */

/*  ===================================================================== */
/* Subroutine */ int zdrot_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer i__, ix, iy;
    static doublecomplex ctemp;


/*  -- Reference BLAS level1 routine (version 3.7.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = i__;
	    z__2.r = *c__ * cx[i__2].r, z__2.i = *c__ * cx[i__2].i;
	    i__3 = i__;
	    z__3.r = *s * cy[i__3].r, z__3.i = *s * cy[i__3].i;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
	    i__2 = i__;
	    i__3 = i__;
	    z__2.r = *c__ * cy[i__3].r, z__2.i = *c__ * cy[i__3].i;
	    i__4 = i__;
	    z__3.r = *s * cx[i__4].r, z__3.i = *s * cx[i__4].i;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
	    i__2 = i__;
	    cx[i__2].r = ctemp.r, cx[i__2].i = ctemp.i;
	}
    } else {

/*        code for unequal increments or equal increments not equal */
/*          to 1 */

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
	    i__2 = ix;
	    z__2.r = *c__ * cx[i__2].r, z__2.i = *c__ * cx[i__2].i;
	    i__3 = iy;
	    z__3.r = *s * cy[i__3].r, z__3.i = *s * cy[i__3].i;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
	    i__2 = iy;
	    i__3 = iy;
	    z__2.r = *c__ * cy[i__3].r, z__2.i = *c__ * cy[i__3].i;
	    i__4 = ix;
	    z__3.r = *s * cx[i__4].r, z__3.i = *s * cx[i__4].i;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
	    i__2 = ix;
	    cx[i__2].r = ctemp.r, cx[i__2].i = ctemp.i;
	    ix += *incx;
	    iy += *incy;
	}
    }
    return 0;
} /* zdrot_ */

