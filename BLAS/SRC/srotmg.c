/* srotmg.f -- translated by f2c (version 20160102).
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

/* > \brief \b SROTMG */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SROTMG(SD1,SD2,SX1,SY1,SPARAM) */

/*       .. Scalar Arguments .. */
/*       REAL SD1,SD2,SX1,SY1 */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SPARAM(5) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS */
/* >    THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(SD1)*SX1,SQRT(SD2)*>    SY2)**T. */
/* >    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */
/* > */
/* >    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0 */
/* > */
/* >      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0) */
/* >    H=(          )    (          )    (          )    (          ) */
/* >      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0). */
/* >    LOCATIONS 2-4 OF SPARAM CONTAIN SH11,SH21,SH12, AND SH22 */
/* >    RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE */
/* >    VALUE OF SPARAM(1) ARE NOT STORED IN SPARAM.) */
/* > */
/* >    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE */
/* >    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE */
/* >    OF SD1 AND SD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in,out] SD1 */
/* > \verbatim */
/* >          SD1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in,out] SD2 */
/* > \verbatim */
/* >          SD2 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in,out] SX1 */
/* > \verbatim */
/* >          SX1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] SY1 */
/* > \verbatim */
/* >          SY1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] SPARAM */
/* > \verbatim */
/* >          SPARAM is REAL array, dimension (5) */
/* >     SPARAM(1)=SFLAG */
/* >     SPARAM(2)=SH11 */
/* >     SPARAM(3)=SH21 */
/* >     SPARAM(4)=SH12 */
/* >     SPARAM(5)=SH22 */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup single_blas_level1 */

/*  ===================================================================== */
/* Subroutine */ int srotmg_(real *sd1, real *sd2, real *sx1, real *sy1, real 
	*sparam)
{
    /* Initialized data */

    static real zero = 0.f;
    static real one = 1.f;
    static real two = 2.f;
    static real gam = 4096.f;
    static real gamsq = 16777200.f;
    static real rgamsq = 5.96046e-8f;

    /* System generated locals */
    real r__1;

    /* Local variables */
    static real su, sp1, sp2, sq1, sq2, sh11, sh12, sh21, sh22, sflag, stemp;


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
/*     .. Data statements .. */

    /* Parameter adjustments */
    --sparam;

    /* Function Body */
/*     .. */
    if (*sd1 < zero) {
/*        GO ZERO-H-D-AND-SX1.. */
	sflag = -one;
	sh11 = zero;
	sh12 = zero;
	sh21 = zero;
	sh22 = zero;

	*sd1 = zero;
	*sd2 = zero;
	*sx1 = zero;
    } else {
/*        CASE-SD1-NONNEGATIVE */
	sp2 = *sd2 * *sy1;
	if (sp2 == zero) {
	    sflag = -two;
	    sparam[1] = sflag;
	    return 0;
	}
/*        REGULAR-CASE.. */
	sp1 = *sd1 * *sx1;
	sq2 = sp2 * *sy1;
	sq1 = sp1 * *sx1;

	if (dabs(sq1) > dabs(sq2)) {
	    sh21 = -(*sy1) / *sx1;
	    sh12 = sp2 / sp1;

	    su = one - sh12 * sh21;

	    if (su > zero) {
		sflag = zero;
		*sd1 /= su;
		*sd2 /= su;
		*sx1 *= su;
	    }
	} else {
	    if (sq2 < zero) {
/*              GO ZERO-H-D-AND-SX1.. */
		sflag = -one;
		sh11 = zero;
		sh12 = zero;
		sh21 = zero;
		sh22 = zero;

		*sd1 = zero;
		*sd2 = zero;
		*sx1 = zero;
	    } else {
		sflag = one;
		sh11 = sp1 / sp2;
		sh22 = *sx1 / *sy1;
		su = one + sh11 * sh22;
		stemp = *sd2 / su;
		*sd2 = *sd1 / su;
		*sd1 = stemp;
		*sx1 = *sy1 * su;
	    }
	}
/*     PROCESURE..SCALE-CHECK */
	if (*sd1 != zero) {
	    while(*sd1 <= rgamsq || *sd1 >= gamsq) {
		if (sflag == zero) {
		    sh11 = one;
		    sh22 = one;
		    sflag = -one;
		} else {
		    sh21 = -one;
		    sh12 = one;
		    sflag = -one;
		}
		if (*sd1 <= rgamsq) {
/* Computing 2nd power */
		    r__1 = gam;
		    *sd1 *= r__1 * r__1;
		    *sx1 /= gam;
		    sh11 /= gam;
		    sh12 /= gam;
		} else {
/* Computing 2nd power */
		    r__1 = gam;
		    *sd1 /= r__1 * r__1;
		    *sx1 *= gam;
		    sh11 *= gam;
		    sh12 *= gam;
		}
	    }
	}
	if (*sd2 != zero) {
	    while(dabs(*sd2) <= rgamsq || dabs(*sd2) >= gamsq) {
		if (sflag == zero) {
		    sh11 = one;
		    sh22 = one;
		    sflag = -one;
		} else {
		    sh21 = -one;
		    sh12 = one;
		    sflag = -one;
		}
		if (dabs(*sd2) <= rgamsq) {
/* Computing 2nd power */
		    r__1 = gam;
		    *sd2 *= r__1 * r__1;
		    sh21 /= gam;
		    sh22 /= gam;
		} else {
/* Computing 2nd power */
		    r__1 = gam;
		    *sd2 /= r__1 * r__1;
		    sh21 *= gam;
		    sh22 *= gam;
		}
	    }
	}
    }
    if (sflag < zero) {
	sparam[2] = sh11;
	sparam[3] = sh21;
	sparam[4] = sh12;
	sparam[5] = sh22;
    } else if (sflag == zero) {
	sparam[3] = sh21;
	sparam[4] = sh12;
    } else {
	sparam[2] = sh11;
	sparam[5] = sh22;
    }
    sparam[1] = sflag;
    return 0;
} /* srotmg_ */

