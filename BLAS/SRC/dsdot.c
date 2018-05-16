/* dsdot.f -- translated by f2c (version 20160102).
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

/* > \brief \b DSDOT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DSDOT(N,SX,INCX,SY,INCY) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SX(*),SY(*) */
/*       .. */

/*    AUTHORS */
/*    ======= */
/*    Lawson, C. L., (JPL), Hanson, R. J., (SNLA), */
/*    Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Compute the inner product of two vectors with extended */
/* > precision accumulation and result. */
/* > */
/* > Returns D.P. dot product accumulated in D.P., for S.P. SX and SY */
/* > DSDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY), */
/* > where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
/* > defined in a similar way using INCY. */
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
/* >          SX is REAL array, dimension(N) */
/* >         single precision vector with N elements */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          storage spacing between elements of SX */
/* > \endverbatim */
/* > */
/* > \param[in] SY */
/* > \verbatim */
/* >          SY is REAL array, dimension(N) */
/* >         single precision vector with N elements */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >         storage spacing between elements of SY */
/* > \endverbatim */
/* > */
/* > \result DSDOT */
/* > \verbatim */
/* >          DSDOT is DOUBLE PRECISION */
/* >         DSDOT  double precision dot product (zero if N.LE.0) */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup double_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > \endverbatim */

/* > \par References: */
/*  ================ */
/* > */
/* > \verbatim */
/* > */
/* > */
/* >  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/* >  Krogh, Basic linear algebra subprograms for Fortran */
/* >  usage, Algorithm No. 539, Transactions on Mathematical */
/* >  Software 5, 3 (September 1979), pp. 308-323. */
/* > */
/* >  REVISION HISTORY  (YYMMDD) */
/* > */
/* >  791001  DATE WRITTEN */
/* >  890831  Modified array declarations.  (WRB) */
/* >  890831  REVISION DATE from Version 3.2 */
/* >  891214  Prologue converted to Version 4.0 format.  (BAB) */
/* >  920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/* >  920501  Reformatted the REFERENCES section.  (WRB) */
/* >  070118  Reformat to LAPACK style (JL) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
doublereal dsdot_(integer *n, real *sx, integer *incx, real *sy, integer *
	incy)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    static integer i__, ns, kx, ky;


/*  -- Reference BLAS level1 routine (version 3.7.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Authors: */
/*  ======== */
/*  Lawson, C. L., (JPL), Hanson, R. J., (SNLA), */
/*  Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL) */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    ret_val = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == *incy && *incx > 0) {

/*     Code for equal, positive, non-unit increments. */

	ns = *n * *incx;
	i__1 = ns;
	i__2 = *incx;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    ret_val += (doublereal) sx[i__] * (doublereal) sy[i__];
	}
    } else {

/*     Code for unequal or nonpositive increments. */

	kx = 1;
	ky = 1;
	if (*incx < 0) {
	    kx = (1 - *n) * *incx + 1;
	}
	if (*incy < 0) {
	    ky = (1 - *n) * *incy + 1;
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ret_val += (doublereal) sx[kx] * (doublereal) sy[ky];
	    kx += *incx;
	    ky += *incy;
	}
    }
    return ret_val;
} /* dsdot_ */

