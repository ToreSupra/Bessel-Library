#pragma once

// Bessel library: A C++ library with routines to evaluate Bessel functions of real or complex arguments (https://github.com/jodesarro/bessel-library)
// MIT License Copyright (c) 2021 Jhonas Olivati de Sarro (https://github.com/jodesarro/bessel-library/blob/main/LICENSE)

// See Refs. [[1–3](#references)] for more information concerning Bessel functions and their computation.

//-v------------------------------------------------------------------------
// CHANGELOG
// 0.1.3 Jun 09, 2024 (current version)
//  - Inclusion of Hankel functions cyl_h1 and cyl_h2 for integer or real orders and real or complex arguments.
//  - Inclusion of Airy functions airy_ai and airy_bi for real or complex arguments.
//  - Functions mod_i and mod_k were consistently renamed to cyl_i and cyl_k.
//  - The flags now print the number of components set to zero due to underflow.
//  - Routines zairy_, zbesh_, zbesj_, zbesy_, zbesi_, zbesk_ and zbiry_, based on version 930101 of D. E. Amos routines (https://doi.org/10.1145/7921.214331, ACM domain)
//      were changed (reverted) to be based on slatec ([3], public domain) versions to avoid copyright conflicts between ACM and MIT licenses and permissions. The
//      versions 0.1.1 and 0.1.2 of this code, and github commits related to them, shall be deleted and must be disconsidered and discarded by all users.
//  - Revision and reorganization of all slatec functions.
//  - Creation of functions d1mach and i1mach to make easier to compare with original slatec versions. 
// 0.1.2 Jun 06, 2024
//  - Inclusion of modified Bessel functions mod_i and mod_k for integer or real orders and real or complex arguments.
// 0.1.1 May 27, 2024
//  - Routines zairy_, zbesh_, zbesj_, and zbesy_, updated to the version 930101 of D. E. Amos routines (https://doi.org/10.1145/7921.214331).
//  - Inclusion of routines zbesi_, zbesk_, zbiry_ accordingly to version 930101 of D. E. Amos routines (https://doi.org/10.1145/7921.214331).
//  - Inclusion of C++ callable functions to overload cyl_j and cyl_y for real arguments.
//  - Static declarations removed for thread safety.
// 0.1.0 May 26, 2024
//  - Routines for cyl_j based on Ref. [2] were replaced by D. E. Amos Fortran 77 routines of SLATEC library [3].
//  - D. E. Amos routines zairy_.f, zbesh_.f, zbesj_.f, zbesy_.f, and all their dependencies, were converted to C using f2c (Availabe at: https://www.netlib.org/f2c/. Accessed: May 25, 2024).
//  - Replacement of all functions d1mach amd i1mach by C macros of float.h.
//  - Corrections of the translated f2c version and elimination of external dependencies.
//  - Reorganization of the whole code to be easily callable from C++.
//  - Inclusion of cylindrical Bessel functions of the second kind (or Neumann functions) cyl_y.
//  - Calculation of negative orders for cyl_j and cyl_y through Eq. (5.5.4) of Ref. [2].
//  - Now, cyl, Bessel functions of the first and second kinds, cyl_j and cyl_y, are available also for real (positive or negative) orders.
//  - Inclusion of cyl_j and cyl_y that returns an array of an int sequence of orders.
//  - Inclusion of parameters to print flag messages, and to return scaled versions of cyl_j and cyl_y.
//  - Inclusion of namespace bessel::slatec to call all slatec routines.
// 0.0.0 until May 12, 2024 
//  - Routines for cylindrical Bessel functions of the first kind and int order written based on Ref. [2].
// CHANGELOG
//-^------------------------------------------------------------------------


//-v------------------------------------------------------------------------
// REFERENCES
// [1] M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions With Formulas,
//      Graphs, and Mathematical Tables. Washington, D. C.: National Bureau of Standards, 1972.
// [2] S. Zhang and J. Jin, Computation of Special Functions. New York: Wiley, 1996.
// [3] SLATEC Common Mathematical Library, Version 4.1, July 1993. Comprehensive software library containing
//      over 1400 general purpose mathematical and statistical routines written in Fortran 77. Available
//      at https://www.netlib.org/slatec/ (Accessed: May 25, 2024).
// REFERENCES
//-^------------------------------------------------------------------------


//-v------------------------------------------------------------------------
// FORTRAN TRANSLATED TO C CODE

//-v------------------------------------------------------------------------
// C LIBRARIES
#define _USE_MATH_DEFINES // For M_PI, and other constant macros.
#include <math.h>   // For functions such as sin(), cos(), abs(), ... .
#include <float.h>  // For macro constants such as DBL_MIN, DBL_MAX_EXP, ... .
// C LIBRARIES
//-^------------------------------------------------------------------------


namespace bessel::slatec
{

//-v------------------------------------------------------------------------
// MAIN D. E. AMOS (SLATEC) ROUTINES DECLARATIONS
int zbesj_(double *zr, double *zi, double *fnu, int *kode, int *n, double *cyr, double *cyi, int *nz, int *ierr);
int zbesy_(double *zr, double *zi, double *fnu, int *kode, int *n, double *cyr, double *cyi, int *nz, double *cwrkr, double *cwrki, int *ierr);
int zbesh_(double *zr, double *zi, double *fnu, int *kode, int *m, int *n, double *cyr, double *cyi, int *nz, int *ierr);
int zbesi_(double *zr, double *zi, double *fnu, int *kode, int *n, double *cyr, double *cyi, int *nz, int *ierr);
int zbesk_(double *zr, double *zi, double *fnu, int *kode, int *n, double *cyr, double *cyi, int *nz, int *ierr);
int zairy_(double *zr, double *zi, int *id, int *kode, double *air, double *aii, int *nz, int *ierr);
int zbiry_(double *zr, double *zi, int *id, int *kode, double *bir, double *bii, int *ierr);
// MAIN D. E. AMOS (SLATEC) ROUTINES DECLARATIONS
//-^------------------------------------------------------------------------


//-v------------------------------------------------------------------------
// DEPENDENCY AMOS/SLATEC ROUTINES DECLARATIONS
double zabs_(double *zr, double *zi);
int zexp_(double *ar, double *ai, double *br, double *bi);
int zdiv_(double *ar, double *ai, double *br, double *bi, double *cr, double *ci);
int zsqrt_(double *ar, double *ai, double *br, double *bi);
int zlog_(double *ar, double *ai, double *br, double *bi, int *ierr);
int zs1s2_(double *zrr, double *zri, double *s1r, double *s1i, double *s2r, double *s2i, int *nz, double *ascle, double *alim, int *iuf);
int zasyi_(double *zr, double *zi, double *fnu, int *kode, int *n, double *yr, double *yi, int *nz, double *rl, double *tol, double *elim, double *alim);
int zacai_(double *zr, double *zi, double *fnu,	int *kode, int *mr, int *n, double *yr, double *yi, int *nz, double *rl, double *tol, double *elim,	double *alim);
int zuni1_(double *zr, double *zi, double *fnu, int *kode, int *n, double *yr, double *yi, int *nz, int *nlast, double *fnul, double *tol, double *elim, double *alim);
int zuni2_(double *zr, double *zi, double *fnu, int *kode, int *n, double *yr, double *yi, int *nz, int *nlast, double *fnul, double *tol, double *elim, double *alim);
int zbuni_(double *zr, double *zi, double *fnu, int *kode, int *n, double *yr, double *yi, int *nz, int *nui, int *nlast, double *fnul, double *tol, double *elim, double *alim);
int zmlri_(double *zr, double *zi, double *fnu, int *kode, int *n, double *yr, double *yi, int *nz, double *tol);
int zmlt_(double *ar, double *ai, double *br, double *bi, double *cr, double *ci);
double dgamln_(double *z__, int *ierr);
int zseri_(double *zr, double *zi, double *fnu, int *kode, int *n, double *yr, double *yi, int *nz, double *tol, double *elim, double *alim);
int zunik_(double *zrr, double *zri, double *fnu, int *ikflg, int *ipmtr, double *tol, int *init, double *phir, double *phii, double *zeta1r, double *zeta1i, double *zeta2r, double *zeta2i, double *sumr, double *sumi, double *cwrkr, double *cwrki);
int zunhj_(double *zr, double *zi, double *fnu, int *ipmtr, double *tol, double *phir, double *phii, double *argr, double *argi, double *zeta1r, double *zeta1i, double *zeta2r, double *zeta2i, double *asumr, double *asumi, double *bsumr, double *bsumi);
int zuchk_(double *yr, double *yi, int *nz, double *ascle, double *tol);
int zuoik_(double *zr, double *zi, double *fnu, int *kode, int *ikflg, int *n, double *yr, double *yi, int *nuf, double *tol, double *elim, double *alim);
int zbknu_(double *zr, double *zi, double *fnu, int *kode, int *n, double *yr, double *yi, int *nz, double *tol, double *elim, double *alim);
int zrati_(double *zr, double *zi, double *fnu, int *n, double *cyr, double *cyi, double *tol);
int zwrsk_(double *zrr, double *zri, double *fnu, int *kode, int *n, double *yr, double *yi, int *nz, double *cwr, double *cwi, double *tol, double *elim, double *alim);
int zbinu_(double *zr, double *zi, double *fnu, int *kode, int *n, double *cyr, double *cyi, int *nz, double *rl, double *fnul, double *tol, double *elim, double *alim);
int zshch_(double *zr, double *zi, double *cshr,	double *cshi, double *cchr, double *cchi);
int zkscl_(double *zrr, double *zri, double *fnu, int *n, double *yr, double *yi, int *nz, double *rzr, double *rzi, double *ascle, double *tol, double *elim);
int zacon_(double *zr, double *zi, double *fnu, int *kode, int *mr, int *n, double *yr, double *yi, int *nz, double *rl, double *fnul, double *tol, double *elim, double *alim);
int zbunk_(double *zr, double *zi, double *fnu, int *kode, int *mr, int *n, double *yr, double *yi, int *nz, double *tol, double *elim, double *alim);
int zunk1_(double *zr, double *zi, double *fnu, int *kode, int *mr, int *n, double *yr, double *yi, int *nz, double *tol, double *elim, double *alim);
int zunk2_(double *zr, double *zi, double *fnu, int *kode, int *mr, int *n, double *yr, double *yi, int *nz, double *tol, double *elim, double *alim);
// DEPENDENCY AMOS/SLATEC ROUTINES DECLARATIONS
//-^------------------------------------------------------------------------


//-v------------------------------------------------------------------------
// TABLE OF GLOBAL CONSTANT VALUES
int c__0 = 0;
int c__1 = 1;
int c__2 = 2;
int c__4 = 4;
int c__5 = 5;
int c__9 = 9;
int c__14 = 14;
int c__15 = 15;
int c__16 = 16;
double c_b10 = .5;
double c_b11 = 0.;
// TABLE OF GLOBAL CONSTANT VALUES
//-^------------------------------------------------------------------------


//-v------------------------------------------------------------------------
// FUNCTIONS OF GLOBAL CONSTANT VALUES
double d1mach_(int *c__n);

int i1mach_(int *c__n);
// FUNCTIONS OF GLOBAL CONSTANT VALUES
//-^------------------------------------------------------------------------


//-v------------------------------------------------------------------------
// DEPENDENCY C ROUTINES
inline double max(double x, double y) { return((x) > (y) ? x : y); }
inline double min(double x, double y) { return((x) < (y) ? x : y); }
inline double d_sign(double *x, double *y) { return ((*y >= 0.) ? abs(*x) : -abs(*x)); }
inline double pow_dd(double *x, double *y) { return pow(*x,*y); }
// DEPENDENCY C ROUTINES
//-^------------------------------------------------------------------------


//-v------------------------------------------------------------------------
// DEPENDENCY AMOS/SLATEC ROUTINES ORIGINALLY TRANSLATED WITH F2C
/* zabs.f, zexp.f, zdiv.f, zsqrt.f, zlog.f,
    zs1s2.f, zasyi.f, zacai.f, zuni1.f,
    zuni2.f, zbuni.f, zmlri.f, zmlt.f, dgamln.f,
    zseri.f, zunik.f, zunhj.f, zuchk.f, zuoik.f,
    zbknu.f, zrati.f, zwrsk.f, zbinu.f, zshch.f,
    zkscl.f, zacon.f, zbunk.f, zunk1.f, zunk2.f
    -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

/* DECK ZABS */
/* Subroutine */
double zabs_(double *zr, double *zi);

/* DECK ZEXP */
/* Subroutine */
int zexp_(double *ar, double *ai, double *br, 
	double *bi);

/* DECK ZDIV */
/* Subroutine */
int zdiv_(double *ar, double *ai, double *br, 
	double *bi, double *cr, double *ci);

/* DECK ZSQRT */
/* Subroutine */
int zsqrt_(double *ar, double *ai, double *br, 
	double *bi);

/* DECK ZLOG */
/* Subroutine */
int zlog_(double *ar, double *ai, double *br, 
	double *bi, int *ierr);

/* DECK ZS1S2 */
/* Subroutine */
int zs1s2_(double *zrr, double *zri, double *s1r,
	 double *s1i, double *s2r, double *s2i, int *nz, 
	double *ascle, double *alim, int *iuf);

/* DECK ZASYI */
/* Subroutine */
int zasyi_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *yr, double *yi, int *
	nz, double *rl, double *tol, double *elim, double *
	alim);

/* DECK ZACAI */
/* Subroutine */
int zacai_(double *zr, double *zi, double *fnu, 
	int *kode, int *mr, int *n, double *yr, double *
	yi, int *nz, double *rl, double *tol, double *elim, 
	double *alim);

/* DECK ZUNI1 */
/* Subroutine */
int zuni1_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *yr, double *yi, int *
	nz, int *nlast, double *fnul, double *tol, double *
	elim, double *alim);

/* DECK ZUNI2 */
/* Subroutine */
int zuni2_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *yr, double *yi, int *
	nz, int *nlast, double *fnul, double *tol, double *
	elim, double *alim);

/* DECK ZBUNI */
/* Subroutine */
int zbuni_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *yr, double *yi, int *
	nz, int *nui, int *nlast, double *fnul, double *tol, 
	double *elim, double *alim);

/* DECK ZMLRI */
/* Subroutine */
int zmlri_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *yr, double *yi, int *
	nz, double *tol);

/* DECK ZMLT */
/* Subroutine */
int zmlt_(double *ar, double *ai, double *br, 
	double *bi, double *cr, double *ci);

/* DECK DGAMLN */
/* Subroutine */
double dgamln_(double *z__, int *ierr);

/* DECK ZSERI */
/* Subroutine */
int zseri_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *yr, double *yi, int *
	nz, double *tol, double *elim, double *alim);

/* DECK ZUNIK */
/* Subroutine */
int zunik_(double *zrr, double *zri, double *fnu,
	 int *ikflg, int *ipmtr, double *tol, int *init, 
	double *phir, double *phii, double *zeta1r, double *
	zeta1i, double *zeta2r, double *zeta2i, double *sumr, 
	double *sumi, double *cwrkr, double *cwrki);

/* DECK ZUNHJ */
/* Subroutine */
int zunhj_(double *zr, double *zi, double *fnu, 
	int *ipmtr, double *tol, double *phir, double *phii, 
	double *argr, double *argi, double *zeta1r, double *
	zeta1i, double *zeta2r, double *zeta2i, double *asumr, 
	double *asumi, double *bsumr, double *bsumi);

/* DECK ZUCHK */
/* Subroutine */
int zuchk_(double *yr, double *yi, int *nz, 
	double *ascle, double *tol);

/* DECK ZUOIK */
/* Subroutine */
int zuoik_(double *zr, double *zi, double *fnu, 
	int *kode, int *ikflg, int *n, double *yr, double 
	*yi, int *nuf, double *tol, double *elim, double *
	alim);

/* DECK ZBKNU */
/* Subroutine */
int zbknu_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *yr, double *yi, int *
	nz, double *tol, double *elim, double *alim);

/* DECK ZRATI */
/* Subroutine */
int zrati_(double *zr, double *zi, double *fnu, 
	int *n, double *cyr, double *cyi, double *tol);

/* DECK ZWRSK */
/* Subroutine */
int zwrsk_(double *zrr, double *zri, double *fnu,
	 int *kode, int *n, double *yr, double *yi, int *
	nz, double *cwr, double *cwi, double *tol, double *
	elim, double *alim);

/* DECK ZBINU */
/* Subroutine */
int zbinu_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *cyr, double *cyi, int *
	nz, double *rl, double *fnul, double *tol, double *
	elim, double *alim);

/* DECK ZSHCH */
/* Subroutine */
int zshch_(double *zr, double *zi, double *cshr, 
	double *cshi, double *cchr, double *cchi);

/* DECK ZKSCL */
/* Subroutine */
int zkscl_(double *zrr, double *zri, double *fnu,
	 int *n, double *yr, double *yi, int *nz, double *
	rzr, double *rzi, double *ascle, double *tol, double *
	elim);

/* DECK ZACON */
/* Subroutine */
int zacon_(double *zr, double *zi, double *fnu, 
	int *kode, int *mr, int *n, double *yr, double *
	yi, int *nz, double *rl, double *fnul, double *tol, 
	double *elim, double *alim);

/* DECK ZBUNK */
/* Subroutine */
int zbunk_(double *zr, double *zi, double *fnu, 
	int *kode, int *mr, int *n, double *yr, double *
	yi, int *nz, double *tol, double *elim, double *alim);

/* DECK ZUNK1 */
/* Subroutine */
int zunk1_(double *zr, double *zi, double *fnu, 
	int *kode, int *mr, int *n, double *yr, double *
	yi, int *nz, double *tol, double *elim, double *alim);

/* DECK ZUNK2 */
/* Subroutine */
int zunk2_(double *zr, double *zi, double *fnu, 
	int *kode, int *mr, int *n, double *yr, double *
	yi, int *nz, double *tol, double *elim, double *alim);

// DEPENDENCY AMOS/SLATEC ROUTINES ORIGINALLY TRANSLATED WITH F2C
//-^------------------------------------------------------------------------


//-v------------------------------------------------------------------------
// MAIN D. E. AMOS (SLATEC) ROUTINES ORIGINALLY TRANSLATED WITH F2C
/* zbesj.f, zbesy.f, zbesh.f, zbesi.f, zbesk.f,
   zairy.f, zbiry.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

/* DECK ZBESJ */
/* Subroutine */
int zbesj_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *cyr, double *cyi, int *
	nz, int *ierr);

/* DECK ZBESY */
/* Subroutine */
int zbesy_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *cyr, double *cyi, int *
	nz, double *cwrkr, double *cwrki, int *ierr);

/* DECK ZBESH */
/* Subroutine */
int zbesh_(double *zr, double *zi, double *fnu, 
	int *kode, int *m, int *n, double *cyr, double *
	cyi, int *nz, int *ierr);

/* DECK ZBESI */
/* Subroutine */
int zbesi_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *cyr, double *cyi, int *
	nz, int *ierr);

/* DECK ZBESK */
/* Subroutine */
int zbesk_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *cyr, double *cyi, int *
	nz, int *ierr);

/* DECK ZAIRY */
/* Subroutine */
int zairy_(double *zr, double *zi, int *id, 
	int *kode, double *air, double *aii, int *nz, int 
	*ierr);

/* DECK ZBIRY */
/* Subroutine */
int zbiry_(double *zr, double *zi, int *id, 
	int *kode, double *bir, double *bii, int *ierr);

// MAIN D. E. AMOS (SLATEC) ROUTINES ORIGINALLY TRANSLATED WITH F2C
//-^------------------------------------------------------------------------


} //namespace bessel::slatec

// FORTRAN TRANSLATED TO C CODE
//-^------------------------------------------------------------------------


//-v------------------------------------------------------------------------
// C++ CODE

//-v------------------------------------------------------------------------
// C++ LIBRARIES
#include <complex>  // For complex variables. 
#include <iostream> // For std::cerr, std::endl, ... .
// C++ LIBRARIES
//-^------------------------------------------------------------------------


namespace bessel
{

//-v------------------------------------------------------------------------
// C++ AMOS/SLATEC FLAGS
void _flag_zbesj_(const int _ierr, const int _nz);
void _flag_zbesy_(const int _ierr, const int _nz);
void _flag_zbesh_(const int _ierr, const int _nz);
void _flag_zbesi_(const int _ierr, const int _nz);
void _flag_zbesk_(const int _ierr, const int _nz);
void _flag_zairy_(const int _ierr, const int _nz);
void _flag_zbiry_(const int _ierr);
// C++ AMOS/SLATEC FLAGS
//-^------------------------------------------------------------------------


//-v------------------------------------------------------------------------
// C++ CALLABLE BESSEL ROUTINES

// Cylindrical Bessel function of the first kind J_nu(z) 
template<typename T1, typename T2>
std::complex<T2> cyl_j( const T1 _nu, const std::complex<T2> _z, bool _scaled = false, bool _flags = false );

template<typename T1, typename T2>
T2 cyl_j( const T1 _nu, const T2 _z, bool _scaled = false, bool _flags = false );

// Array of cylindrical Bessel functions of the first kind J_nu(z)
template<typename T1, typename T2>
void cyl_j( const T1 _nu, const int _n, const std::complex<T2> _z, std::complex<T2> * _cyl_j, bool _scaled = false, bool _flags = false );

template<typename T1, typename T2>
void cyl_j( const T1 _nu, const int _n, T2 _z, T2 * _cyl_j, bool _scaled = false, bool _flags = false );

// Cylindrical Bessel function of the second kind Y_nu(z)
template<typename T1, typename T2>
std::complex<T2> cyl_y( const T1 _nu, const std::complex<T2> _z, bool _scaled = false, bool _flags = false );

template<typename T1, typename T2>
T2 cyl_y( const T1 _nu, const T2 _z, bool _scaled = false, bool _flags = false );

// Array of cylindrical Bessel functions of the second kind Y_nu(z)
template<typename T1, typename T2>
void cyl_y( const T1 _nu, const int _n, const std::complex<T2> _z, std::complex<T2> * _cyl_y, bool _scaled = false, bool _flags = false );

template<typename T1, typename T2>
void cyl_y( const T1 _nu, const int _n, T2 _z, T2 * _cyl_y, bool _scaled = false, bool _flags = false );

// Hankel function of the first kind H1_nu(z)
template<typename T1, typename T2>
std::complex<T2> cyl_h1( const T1 _nu, const std::complex<T2> _z, bool _scaled = false, bool _flags = false );

template<typename T1, typename T2>
T2 cyl_h1( const T1 _nu, const T2 _z, bool _scaled = false, bool _flags = false );

// Array of Hankel functiond of the first kind H1_nu(z)
template<typename T1, typename T2>
void cyl_h1( const T1 _nu, const int _n, const std::complex<T2> _z, std::complex<T2> * _cyl_h1, bool _scaled = false, bool _flags = false );

template<typename T1, typename T2>
void cyl_h1( const T1 _nu, const int _n, T2 _z, T2 * _cyl_h1, bool _scaled = false, bool _flags = false );

// Hankel function of the second kind H2_nu(z)
template<typename T1, typename T2>
std::complex<T2> cyl_h2( const T1 _nu, const std::complex<T2> _z, bool _scaled = false, bool _flags = false );

template<typename T1, typename T2>
T2 cyl_h2( const T1 _nu, const T2 _z, bool _scaled = false, bool _flags = false );

// Array of Hankel functions of the second kind H2_nu(z)
template<typename T1, typename T2>
void cyl_h2( const T1 _nu, const int _n, const std::complex<T2> _z, std::complex<T2> * _cyl_h2, bool _scaled = false, bool _flags = false );

template<typename T1, typename T2>
void cyl_h2( const T1 _nu, const int _n, T2 _z, T2 * _cyl_h2, bool _scaled = false, bool _flags = false );

// Modified cylindrical Bessel function of the first kind I_nu(z)
template<typename T1, typename T2>
std::complex<T2> cyl_i( const T1 _nu, const std::complex<T2> _z, bool _scaled = false, bool _flags = false );

template<typename T1, typename T2>
T2 cyl_i( const T1 _nu, const T2 _z, bool _scaled = false, bool _flags = false );

// Array of modified cylindrical Bessel functions of the first kind I_nu(z)
template<typename T1, typename T2>
void cyl_i( const T1 _nu, const int _n, const std::complex<T2> _z, std::complex<T2> * _cyl_i, bool _scaled = false, bool _flags = false );

template<typename T1, typename T2>
void cyl_i( const T1 _nu, const int _n, T2 _z, T2 * _cyl_i, bool _scaled = false, bool _flags = false );

// Modified cylindrical Bessel function of the second kind K_nu(z) [Basset or MacDonald function]
template<typename T1, typename T2>
std::complex<T2> cyl_k( const T1 _nu, const std::complex<T2> _z, bool _scaled = false, bool _flags = false );

template<typename T1, typename T2>
T2 cyl_k( const T1 _nu, const T2 _z, bool _scaled = false, bool _flags = false );

// Array of modified cylindrical Bessel functions of the second kind K_nu(z) [Basset or MacDonald functions]
template<typename T1, typename T2>
void cyl_k( const T1 _nu, const int _n, const std::complex<T2> _z, std::complex<T2> * _cyl_k, bool _scaled = false, bool _flags = false );

template<typename T1, typename T2>
void cyl_k( const T1 _nu, const int _n, T2 _z, T2 * _cyl_k, bool _scaled = false, bool _flags = false );

// Airy Ai function
template<typename T1>
std::complex<T1> airy_ai( const std::complex<T1> _z, bool _scaled = false, bool _flags = false );

template<typename T1>
T1 airy_ai( const T1 _z, bool _scaled = false, bool _flags = false );

// Airy Bi function
template<typename T1>
std::complex<T1> airy_bi( const std::complex<T1> _z, bool _scaled = false, bool _flags = false );

template<typename T1>
T1 airy_bi( const T1 _z, bool _scaled = false, bool _flags = false );

// C++ CALLABLE BESSEL ROUTINES
//-^------------------------------------------------------------------------


} // namespace bessel

// C++ CODE
//-^------------------------------------------------------------------------