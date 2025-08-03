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
	// FUNCTIONS OF GLOBAL CONSTANT VALUES
	double d1mach_(int *c__n);

	int i1mach_(int *c__n);
	// FUNCTIONS OF GLOBAL CONSTANT VALUES
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
    std::complex<T2> cyl_j( const T1 _nu, const std::complex<T2> _z, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _nu:     real order.
        //  _z:      complex argument.
        //  _scaled: if true, returns the scaled version J_nu(z)*exp(-|Im(z)|).
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int nz;
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        std::complex<double> z = std::complex<double> (x, y);
        int n = 1;
        double nu = double(_nu);
        double fnu = abs(nu);

        if ( nu >= 0. )
        {
            // Positive order
            double cjr, cji;
            slatec::zbesj_(&x, &y, &fnu, &kode, &n, &cjr, &cji, &nz, &ierr);
            if (_flags) _flag_zbesj_(ierr, nz);
            return std::complex<T2>(cjr, cji);
        }
        else
        {
            // Negative order
            if ( double(floor(fnu)) == fnu )
            {
                // int negative order
                double cjr, cji;
                slatec::zbesj_(&x, &y, &fnu, &kode, &n, &cjr, &cji, &nz, &ierr);
                if (_flags) _flag_zbesj_(ierr, nz);
                return pow(-1., fnu) * std::complex<T2>(cjr, cji);
            }
            else if ( double(floor(fnu)) != fnu && abs(z) == 0. )
            {
                // Non-int negative order with abs(z)=0,
                //  J_nu = complex infinity since Y_fnu = complex infinity.
                double dinf = std::numeric_limits<double>::infinity();
                return std::complex<double>(dinf, dinf);
            }
            else
            {
                // Non-int negative order with abs(z)!=0.
                double cjr, cji;
                slatec::zbesj_(&x, &y, &fnu, &kode, &n, &cjr, &cji, &nz, &ierr);
                if (_flags) _flag_zbesj_(ierr, nz);

                double cyr, cyi;
                double * cwrkr = new double [1];
                double * cwrki = new double [1];
                slatec::zbesy_(&x, &y, &fnu, &kode, &n, &cyr, &cyi, &nz, cwrkr, cwrki, &ierr);
                if (_flags) _flag_zbesy_(ierr, nz);

                //  Eq. (5.5.4) of Zhang and Jin (1996) [2].
                double tmp = fnu*M_PI;
                return ( std::complex<T2>(cjr, cji) )*cos(tmp) - ( std::complex<T2>(cyr, cyi) )*sin(tmp);
            }
        }
    }

    template<typename T1, typename T2>
    T2 cyl_j( const T1 _nu, const T2 _z, bool _scaled = false, bool _flags = false )
    {
        return real( cyl_j( _nu, std::complex<T2>(_z,0.), _scaled, _flags ) ) ;
    }

	// Array of cylindrical Bessel functions of the first kind J_nu(z)
    template<typename T1, typename T2>
    void cyl_j( const T1 _nu, const int _n, const std::complex<T2> _z, std::complex<T2> * _cyl_j, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _nu:     real first order to be calculated.
        //  _n:      int to return an array with the orders _nu, _nu+1, ..., _nu+_n.
        //  _z:      complex argument.
        //  _cyl_j:  array of size _n to output J_nu(z) for the orders _nu, _nu+1, ..., _nu+_n.
        //  _scaled: if true, returns an array of size _n with the scaled version J_nu(z)*exp(-|Im(z)|) for the orders _nu, _nu+1, ..., _nu+_n.
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int nz;
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        std::complex<double> z = std::complex<double> (x, y);
        int n = _n;
        double nu = double(_nu);
        double fnu = abs(nu);
        double nu_m = nu + double(n-1);

        if ( n <= 0 )
        {
            if (_flags) std::cerr << "[bessel::cyl_j] Flag -> Input error: No computation, _n must be an integer greater than zero." << std::endl;
        }
        else if ( nu >= 0. )
        {
            // Array of only positive orders, including the order 0
            double * cjr = new double [n];
            double * cji = new double [n];
            slatec::zbesj_(&x, &y, &fnu, &kode, &n, cjr, cji, &nz, &ierr);
            if (_flags) _flag_zbesj_(ierr, nz);

            for ( int i=0; i<n; i++ )
            {
                _cyl_j[i] = std::complex<T2>(cjr[i], cji[i]);
            }
        }
        else if ( nu_m <= 0. )
        {
            // Array of only negative orders, including the order 0
            if ( double(floor(fnu)) == fnu )
            {
                // int negative orders, or order 0
                double * cjr = new double [n];
                double * cji = new double [n];
                double fnu_m = abs( nu_m );
                slatec::zbesj_(&x, &y, &fnu_m, &kode, &n, cjr, cji, &nz, &ierr);
                if (_flags) _flag_zbesj_(ierr, nz);

                for (int i=0; i<n; i++)
                {
                    int tmp = n-1-i;
                    _cyl_j[i] = pow(-1., double(fnu_m+tmp)) * std::complex<T2>(cjr[tmp], cji[tmp]);
                }
            }
            else if ( double(floor(nu_m)) != nu_m && abs(z) == 0. )
            {
                // Non-int negative orders with abs(z)=0,
                //  Jnu = complex infinity since Yfnu = complex infinity.
                double dinf = std::numeric_limits<double>::infinity();
                for (int i=0; i<n; i++)
                {
                    _cyl_j[i] = std::complex<T2>(dinf, dinf);
                }
            }
            else
            {
                // Non-int negative orders with abs(z)!=0,
                double * cjr = new double [n];
                double * cji = new double [n];
                double fnu_m = abs( nu_m );
                slatec::zbesj_(&x, &y, &fnu_m, &kode, &n, cjr, cji, &nz, &ierr);
                if (_flags) _flag_zbesj_(ierr, nz);

                double * cyr = new double [n];
                double * cyi = new double [n];
                double * cwrkr = new double [n];
                double * cwrki = new double [n];
                slatec::zbesy_(&x, &y, &fnu_m, &kode, &n, cyr, cyi, &nz, cwrkr, cwrki, &ierr);
                if (_flags) _flag_zbesy_(ierr, nz);

                for (int i=0; i<n; i++)
                {
                    //  Eq. (5.5.4) of Zhang and Jin (1996) [2].
                    int tmp1 = n-1-i;
                    double tmp2 = ( fnu_m + double(tmp1) )*M_PI;
                    _cyl_j[i] = ( std::complex<T2>(cjr[tmp1], cji[tmp1]) )*cos(tmp2) - ( std::complex<T2>(cyr[tmp1], cyi[tmp1]) )*sin(tmp2);
                }
            }
        }
        else
        {
            // Array of negative and positive orders
            int n_m = floor( abs(nu) ) + 1;
            int n_p = n - n_m;
            double fnu_m = abs( nu + double(n_m-1) );

            // The negative orders
            if ( double(floor(fnu)) == fnu )
            {
                // int negative orders, or order 0
                double * cjr_m = new double [n_m];
                double * cji_m = new double [n_m];
                slatec::zbesj_(&x, &y, &fnu_m, &kode, &n_m, cjr_m, cji_m, &nz, &ierr);
                if (_flags) _flag_zbesj_(ierr, nz);

                for (int i=0; i<n_m; i++)
                {
                    int tmp = n_m-1-i;
                    _cyl_j[i] = pow(-1., double(fnu_m+tmp)) * std::complex<T2>(cjr_m[tmp], cji_m[tmp]);
                }
            }
            else if ( double(floor(nu_m)) != nu_m && abs(z) == 0. )
            {
                // Non-int negative orders with abs(z)=0,
                //  Jnu = complex infinity since Yfnu = complex infinity.
                double dinf = std::numeric_limits<double>::infinity();
                for (int i=0; i<n_m; i++)
                {
                    // Complex infinity
                    _cyl_j[i] = std::complex<T2>(dinf, dinf);
                }
            }
            else
            {
                // Non-int negative orders with abs(z)!=0,
                double * cjr_m = new double [n_m];
                double * cji_m = new double [n_m];
                slatec::zbesj_(&x, &y, &fnu_m, &kode, &n_m, cjr_m, cji_m, &nz, &ierr);
                if (_flags) _flag_zbesj_(ierr, nz);

                double * cyr_m = new double [n_m];
                double * cyi_m = new double [n_m];
                double * cwrkr = new double [n_m];
                double * cwrki = new double [n_m];
                slatec::zbesy_(&x, &y, &fnu_m, &kode, &n_m, cyr_m, cyi_m, &nz, cwrkr, cwrki, &ierr);
                if (_flags) _flag_zbesy_(ierr, nz);

                for (int i=0; i<n_m; i++)
                {
                    //  Eq. (5.5.4) of Zhang and Jin (1996) [2].
                    int tmp1 = n_m-1-i;
                    double tmp2 = ( fnu_m + double(tmp1) )*M_PI;
                    _cyl_j[i] = ( std::complex<T2>(cjr_m[tmp1], cji_m[tmp1]) )*cos(tmp2) - ( std::complex<T2>(cyr_m[tmp1], cyi_m[tmp1]) )*sin(tmp2);
                }
            }

            // The positive orders
            double * cjr_p = new double [n_p];
            double * cji_p = new double [n_p];
            double fnu_p = nu + double(n_m);
            slatec::zbesj_(&x, &y, &fnu_p, &kode, &n_p, cjr_p, cji_p, &nz, &ierr);
            if (_flags) _flag_zbesj_(ierr, nz);

            for (int i=n_m; i<n; i++)
            {
                int tmp1 = i-n_m;
                _cyl_j[i] = std::complex<T2>(cjr_p[tmp1], cji_p[tmp1]);
            }
        }
    }

    template<typename T1, typename T2>
    void cyl_j( const T1 _nu, const int _n, T2 _z, T2 * _cyl_j, bool _scaled = false, bool _flags = false )
    {
        std::complex<T2> * tmp_cylj = new std::complex<T2> [_n];
        cyl_j( _nu, _n, std::complex<T2>(_z,0.), tmp_cylj, _scaled, _flags );
        for (int i=0; i<_n; i++)
        {
            _cyl_j[i] = real(tmp_cylj[i]);
        }
    }

	// Cylindrical Bessel function of the second kind Y_nu(z)
    template<typename T1, typename T2>
    std::complex<T2> cyl_y( const T1 _nu, const std::complex<T2> _z, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _nu:     real order.
        //  _z:      complex argument.
        //  _scaled: if true, returns the scaled version Y_nu(z)*exp(-|Im(z)|).
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int nz;
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        std::complex<double> z = std::complex<double> (x, y);
        int n = 1;
        double nu = double(_nu);
        double fnu = abs(nu);

        if ( abs(z) == 0. )
        {
            double dinf = std::numeric_limits<double>::infinity();
            if ( nu == 0. )
            {
                // Infinity
                return std::complex<double>(-dinf, 0.);
            }
            else
            {
                // Complex infinity
                return std::complex<double>(dinf, dinf);
            }
        }
        else if ( nu > 0. )
        {
            // Positive order
            double cyr, cyi;
            double * cwrkr = new double [1];
            double * cwrki = new double [1];
            slatec::zbesy_(&x, &y, &fnu, &kode, &n, &cyr, &cyi, &nz, cwrkr, cwrki, &ierr);
            if (_flags) _flag_zbesy_(ierr, nz);
            return std::complex<T2>(cyr, cyi);
        }
        else
        {
            // Negative order
            if ( double(floor(fnu)) == fnu )
            {
                // int negative order
                double cyr, cyi;
                double * cwrkr = new double [1];
                double * cwrki = new double [1];
                slatec::zbesy_(&x, &y, &fnu, &kode, &n, &cyr, &cyi, &nz, cwrkr, cwrki, &ierr);
                if (_flags) _flag_zbesy_(ierr, nz);
                return pow(-1., fnu) * std::complex<T2>(cyr, cyi);
            }
            else
            {
                // Non-int negative order
                double cyr, cyi;
                double * cwrkr = new double [1];
                double * cwrki = new double [1];
                slatec::zbesy_(&x, &y, &fnu, &kode, &n, &cyr, &cyi, &nz, cwrkr, cwrki, &ierr);
                if (_flags) _flag_zbesy_(ierr, nz);

                double cjr, cji;
                slatec::zbesj_(&x, &y, &fnu, &kode, &n, &cjr, &cji, &nz, &ierr);
                if (_flags) _flag_zbesj_(ierr, nz);

                //  Eq. (5.5.4) of Zhang and Jin (1996) [2].
                double tmp = fnu*M_PI;
                return ( std::complex<T2>(cjr, cji) )*sin(tmp) + ( std::complex<T2>(cyr, cyi) )*cos(tmp);
            }
        }
    }

	template<typename T1, typename T2>
    T2 cyl_y( const T1 _nu, const T2 _z, bool _scaled = false, bool _flags = false )
    {
        return real( cyl_y( _nu, std::complex<T2>(_z,0.), _scaled, _flags ) ) ;
    }

	// Array of cylindrical Bessel functions of the second kind Y_nu(z)
    template<typename T1, typename T2>
    void cyl_y( const T1 _nu, const int _n, const std::complex<T2> _z, std::complex<T2> * _cyl_y, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _nu:     real first order to be calculated.
        //  _n:      int to return an array with the orders _nu, _nu+1, ..., _nu+_n.
        //  _z:      complex argument.
        //  _cyl_y:  array of size _n to output Y_nu(z) for the orders _nu, _nu+1, ..., _nu+_n.
        //  _scaled: if true, returns an array of size _n with the scaled version Y_nu(z)*exp(-|Im(z)|) for the orders _nu, _nu+1, ..., _nu+_n.
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int nz;
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        std::complex<double> z = std::complex<double> (x, y);
        int n = _n;
        double nu = double(_nu);
        double fnu = abs(nu);
        double nu_m = nu + double(n-1);

        if ( n <= 0 )
        {
            if (_flags) std::cerr << "[bessel::cyl_y] Flag -> Input error: No computation, _n must be an integer greater than zero." << std::endl;
        }
        else if ( abs(z) == 0. )
        {
            double dinf = std::numeric_limits<double>::infinity();
            for (int i=0; i<n; i++)
            {
                // Complex infinity for all orders
                _cyl_y[i] = std::complex<T2>(dinf, dinf);
            }
            if ( nu_m == floor(nu_m) && n > fnu )
            {
                // Infinity for the int order 0
                _cyl_y[ int(floor(fnu)) ] = std::complex<double>(-dinf, 0.);
            }
        }
        else if ( nu >= 0. )
        {
            // Array of only positive orders, including the order 0, with abs(z)!=0
            double * cyr = new double [n];
            double * cyi = new double [n];
            double * cwrkr = new double [n];
            double * cwrki = new double [n];
            slatec::zbesy_(&x, &y, &fnu, &kode, &n, cyr, cyi, &nz, cwrkr, cwrki, &ierr);
            if (_flags) _flag_zbesy_(ierr, nz);

            for ( int i=0; i<n; i++ )
            {
                _cyl_y[i] = std::complex<T2>(cyr[i], cyi[i]);
            }
        }
        else if ( nu_m <= 0. )
        {
            // Array of only negative orders, including the order 0, with abs(z)!=0
            if ( double(floor(fnu)) == fnu )
            {
                // int negative orders, or order 0, with abs(z)!=0
                double * cyr = new double [n];
                double * cyi = new double [n];
                double * cwrkr = new double [n];
                double * cwrki = new double [n];
                double fnu_m = abs( nu_m );
                slatec::zbesy_(&x, &y, &fnu_m, &kode, &n, cyr, cyi, &nz, cwrkr, cwrki, &ierr);
                if (_flags) _flag_zbesy_(ierr, nz);

                for (int i=0; i<n; i++)
                {
                    int tmp = n-1-i;
                    _cyl_y[i] = pow(-1., double(fnu_m+tmp)) * std::complex<T2>(cyr[tmp], cyi[tmp]);
                }
            }
            else
            {
                // Non-int negative orders with abs(z)!=0
                double fnu_m = abs( nu_m );

                double * cyr = new double [n];
                double * cyi = new double [n];
                double * cwrkr = new double [n];
                double * cwrki = new double [n];
                slatec::zbesy_(&x, &y, &fnu_m, &kode, &n, cyr, cyi, &nz, cwrkr, cwrki, &ierr);
                if (_flags) _flag_zbesy_(ierr, nz);

                double * cjr = new double [n];
                double * cji = new double [n];
                slatec::zbesj_(&x, &y, &fnu_m, &kode, &n, cjr, cji, &nz, &ierr);
                if (_flags) _flag_zbesj_(ierr, nz);

                for (int i=0; i<n; i++)
                {
                    //  Eq. (5.5.4) of Zhang and Jin (1996) [2].
                    int tmp1 = n-1-i;
                    double tmp2 = ( fnu_m + double(tmp1) )*M_PI;
                    _cyl_y[i] = ( std::complex<T2>(cjr[tmp1], cji[tmp1]) )*sin(tmp2) + ( std::complex<T2>(cyr[tmp1], cyi[tmp1]) )*cos(tmp2);
                }
            }
        }
        else
        {
            // Array of negative and positive orders with abs(z)!=0
            int n_m = floor( abs(nu) ) + 1;
            int n_p = n - n_m;
            double fnu_m = abs( nu + double(n_m-1) );

            // The negative orders
            if ( double(floor(fnu)) == fnu )
            {
                // int negative orders, or order 0, with abs(z)!=0
                double * cyr_m = new double [n_m];
                double * cyi_m = new double [n_m];
                double * cwrkr = new double [n_m];
                double * cwrki = new double [n_m];
                slatec::zbesy_(&x, &y, &fnu_m, &kode, &n_m, cyr_m, cyi_m, &nz, cwrkr, cwrki, &ierr);
                if (_flags) _flag_zbesy_(ierr, nz);

                for (int i=0; i<n_m; i++)
                {
                    int tmp = n_m-1-i;
                    _cyl_y[i] = pow(-1., double(fnu_m+tmp)) * std::complex<T2>(cyr_m[tmp], cyi_m[tmp]);
                }
            }
            else
            {
                // Non-int negative orders with abs(z)!=0,
                double * cyr_m = new double [n_m];
                double * cyi_m = new double [n_m];
                double * cwrkr = new double [n_m];
                double * cwrki = new double [n_m];
                slatec::zbesy_(&x, &y, &fnu_m, &kode, &n_m, cyr_m, cyi_m, &nz, cwrkr, cwrki, &ierr);
                if (_flags) _flag_zbesy_(ierr, nz);

                double * cjr_m = new double [n_m];
                double * cji_m = new double [n_m];
                slatec::zbesj_(&x, &y, &fnu_m, &kode, &n_m, cjr_m, cji_m, &nz, &ierr);
                if (_flags) _flag_zbesj_(ierr, nz);

                for (int i=0; i<n_m; i++)
                {
                    //  Eq. (5.5.4) of Zhang and Jin (1996) [2].
                    int tmp1 = n_m-1-i;
                    double tmp2 = ( fnu_m + double(tmp1) )*M_PI;
                    _cyl_y[i] = ( std::complex<T2>(cjr_m[tmp1], cji_m[tmp1]) )*sin(tmp2) + ( std::complex<T2>(cyr_m[tmp1], cyi_m[tmp1]) )*cos(tmp2);
                }
            }

            // The positive orders
            double * cyr_p = new double [n_p];
            double * cyi_p = new double [n_p];
            double * cwrkr = new double [n_p];
            double * cwrki = new double [n_p];
            double fnu_p = nu + double(n_m);
            slatec::zbesy_(&x, &y, &fnu_p, &kode, &n_p, cyr_p, cyi_p, &nz, cwrkr, cwrki, &ierr);
            if (_flags) _flag_zbesy_(ierr, nz);

            for (int i=n_m; i<n; i++)
            {
                int tmp1 = i-n_m;
                _cyl_y[i] = std::complex<T2>(cyr_p[tmp1], cyi_p[tmp1]);
            }
        }
    }

    template<typename T1, typename T2>
    void cyl_y( const T1 _nu, const int _n, T2 _z, T2 * _cyl_y, bool _scaled = false, bool _flags = false )
    {
        std::complex<T2> * tmp_cyly = new std::complex<T2> [_n];
        cyl_y( _nu, _n, std::complex<T2>(_z,0.), tmp_cyly, _scaled, _flags );
        for (int i=0; i<_n; i++)
        {
            _cyl_y[i] = real(tmp_cyly[i]);
        }
    }

	// Hankel function of the first kind H1_nu(z)
    template<typename T1, typename T2>
    std::complex<T2> cyl_h1( const T1 _nu, const std::complex<T2> _z, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _nu:     real order.
        //  _z:      complex argument.
        //  _scaled: if true, returns the scaled version H1_nu(z)*exp(-iz).
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int nz;
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        std::complex<double> z = std::complex<double> (x, y);
        int n = 1;
        double nu = double(_nu);
        double fnu = abs(nu);
        int m = 1;

        if ( abs(z) == 0. )
        {
            // Complex infinity
            double dinf = std::numeric_limits<double>::infinity();
            return std::complex<double>(dinf, dinf);
        }
        else
        {
            double ch1r, ch1i;
            slatec::zbesh_(&x, &y, &fnu, &kode, &m, &n, &ch1r, &ch1i, &nz, &ierr);
            if (_flags) _flag_zbesh_(ierr, nz);
            if ( nu >= 0. )
            {
                // Positive order
                return std::complex<T2>(ch1r, ch1i);
            }
            else
            {
                // Negative order, see Eq. (9.1.6) of Ref. [1].
                std::complex<T2> exp_of_Ipifnu = exp( std::complex<T2> (0.,M_PI*fnu) );
                return std::complex<T2>(ch1r, ch1i)*exp_of_Ipifnu;
            }
        }
    }

    template<typename T1, typename T2>
    T2 cyl_h1( const T1 _nu, const T2 _z, bool _scaled = false, bool _flags = false )
    {
        return real( cyl_h1( _nu, std::complex<T2>(_z,0.), _scaled, _flags ) ) ;
    }

	// Array of Hankel functiond of the first kind H1_nu(z)
    template<typename T1, typename T2>
    void cyl_h1( const T1 _nu, const int _n, const std::complex<T2> _z, std::complex<T2> * _cyl_h1, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _nu:     real first order to be calculated.
        //  _n:      int to return an array with the orders _nu, _nu+1, ..., _nu+_n.
        //  _z:      complex argument.
        //  _cyl_h1: array of size _n to output H1_nu(z) for the orders _nu, _nu+1, ..., _nu+_n.
        //  _scaled: if true, returns an array of size _n with the scaled version H1_nu(z)*exp(-iz). for the orders _nu, _nu+1, ..., _nu+_n.
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int nz;
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        std::complex<double> z = std::complex<double> (x, y);
        int n = _n;
        double nu = double(_nu);
        double fnu = abs(nu);
        double nu_m = nu + double(n-1);
        int m = 1;

        if ( n <= 0 )
        {
            if (_flags) std::cerr << "[bessel::cyl_h1] Flag -> Input error: No computation, _n must be an integer greater than zero." << std::endl;
        }
        else if ( abs(z) == 0. )
        {
            // Complex infinity for all orders
            double dinf = std::numeric_limits<double>::infinity();
            for (int i=0; i<n; i++)
            {
                _cyl_h1[i] = std::complex<T2>(dinf, dinf);
            }
        }
        else if ( nu >= 0. )
        {
            // Array of only positive orders, including the order 0, with abs(z)!=0
            double * ch1r = new double [n];
            double * ch1i = new double [n];
            slatec::zbesh_(&x, &y, &fnu, &kode, &m, &n, ch1r, ch1i, &nz, &ierr);
            if (_flags) _flag_zbesh_(ierr, nz);
            for ( int i=0; i<n; i++ )
            {
                _cyl_h1[i] = std::complex<T2>(ch1r[i], ch1i[i]);
            }
        }
        else if ( nu_m <= 0. )
        {
            // Array of only negative orders, including the order 0, with abs(z)!=0
            double * ch1r_m = new double [n];
            double * ch1i_m = new double [n];
            double fnu_m = abs( nu_m );
            slatec::zbesh_(&x, &y, &fnu_m, &kode, &m, &n, ch1r_m, ch1i_m, &nz, &ierr);
            if (_flags) _flag_zbesh_(ierr, nz);
            for (int i=0; i<n; i++)
            {
                // Eq. (9.1.6) of Ref. [1].
                int tmp1 = n-1-i;
                double tmp2 = ( fnu_m + double(tmp1) )*M_PI;
                std::complex<T2> exp_of_Ipifnu = exp( std::complex<T2> (0.,tmp2) );
                _cyl_h1[i] = std::complex<T2>(ch1r_m[tmp1], ch1i_m[tmp1])*exp_of_Ipifnu;
            }
        }
        else
        {
            // Array of negative and positive orders with abs(z)!=0
            int n_m = floor( abs(nu) ) + 1;
            int n_p = n - n_m;
            double fnu_m = abs( nu + double(n_m-1) );

            // The negative orders
            // int negative orders, or order 0, with abs(z)!=0
            double * ch1r_m = new double [n_m];
            double * ch1i_m = new double [n_m];
            slatec::zbesh_(&x, &y, &fnu_m, &kode, &m, &n_m, ch1r_m, ch1i_m, &nz, &ierr);
            if (_flags) _flag_zbesh_(ierr, nz);
            for (int i=0; i<n_m; i++)
            {
                // Eq. (9.1.6) of Ref. [1].
                int tmp1 = n_m-1-i;
                double tmp2 = ( fnu_m + double(tmp1) )*M_PI;
                std::complex<T2> exp_of_Ipifnu = exp( std::complex<T2> (0.,tmp2) );
                _cyl_h1[i] = std::complex<T2>(ch1r_m[tmp1], ch1i_m[tmp1])*exp_of_Ipifnu;
            }

            // The positive orders
            double * ch1r_p = new double [n_p];
            double * ch1i_p = new double [n_p];
            double fnu_p = nu + double(n_m);
            slatec::zbesh_(&x, &y, &fnu_p, &kode, &m, &n_p, ch1r_p, ch1i_p, &nz, &ierr);
            if (_flags) _flag_zbesh_(ierr, nz);
            for (int i=n_m; i<n; i++)
            {
                int tmp1 = i-n_m;
                _cyl_h1[i] = std::complex<T2>(ch1r_p[tmp1], ch1i_p[tmp1]);
            }
        }
    }

    template<typename T1, typename T2>
    void cyl_h1( const T1 _nu, const int _n, T2 _z, T2 * _cyl_h1, bool _scaled = false, bool _flags = false )
    {
        std::complex<T2> * tmp_cylh1 = new std::complex<T2> [_n];
        cyl_h1( _nu, _n, std::complex<T2>(_z,0.), tmp_cylh1, _scaled, _flags );
        for (int i=0; i<_n; i++)
        {
            _cyl_h1[i] = real(tmp_cylh1[i]);
        }
    }

	// Hankel function of the second kind H2_nu(z)
    template<typename T1, typename T2>
    std::complex<T2> cyl_h2( const T1 _nu, const std::complex<T2> _z, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _nu:     real order.
        //  _z:      complex argument.
        //  _scaled: if true, returns the scaled version H2_nu(z)*exp(+iz).
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int nz;
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        std::complex<double> z = std::complex<double> (x, y);
        int n = 1;
        double nu = double(_nu);
        double fnu = abs(nu);
        int m = 2;

        if ( abs(z) == 0. )
        {
            // Complex infinity
            double dinf = std::numeric_limits<double>::infinity();
            return std::complex<double>(dinf, dinf);
        }
        else
        {
            double ch2r, ch2i;
            slatec::zbesh_(&x, &y, &fnu, &kode, &m, &n, &ch2r, &ch2i, &nz, &ierr);
            if (_flags) _flag_zbesh_(ierr, nz);
            if ( nu >= 0. )
            {
                // Positive order
                return std::complex<T2>(ch2r, ch2i);
            }
            else
            {
                // Negative order, see Eq. (9.1.6) of Ref. [1].
                std::complex<T2> exp_of_mIpifnu = exp( std::complex<T2> (0.,-M_PI*fnu) );
                return std::complex<T2>(ch2r, ch2i)*exp_of_mIpifnu;
            }
        }
    }

    template<typename T1, typename T2>
    T2 cyl_h2( const T1 _nu, const T2 _z, bool _scaled = false, bool _flags = false )
    {
        return real( cyl_h2( _nu, std::complex<T2>(_z,0.), _scaled, _flags ) ) ;
    }

	// Array of Hankel functions of the second kind H2_nu(z)
    template<typename T1, typename T2>
    void cyl_h2( const T1 _nu, const int _n, const std::complex<T2> _z, std::complex<T2> * _cyl_h2, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _nu:     real first order to be calculated.
        //  _n:      int to return an array with the orders _nu, _nu+1, ..., _nu+_n.
        //  _z:      complex argument.
        //  _cyl_h2: array of size _n to output H2_nu(z) for the orders _nu, _nu+1, ..., _nu+_n.
        //  _scaled: if true, returns an array of size _n with the scaled version H2_nu(z)*exp(+iz). for the orders _nu, _nu+1, ..., _nu+_n.
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int nz;
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        std::complex<double> z = std::complex<double> (x, y);
        int n = _n;
        double nu = double(_nu);
        double fnu = abs(nu);
        double nu_m = nu + double(n-1);
        int m = 2;

        if ( n <= 0 )
        {
            if (_flags) std::cerr << "[bessel::cyl_h2] Flag -> Input error: No computation, _n must be an integer greater than zero." << std::endl;
        }
        else if ( abs(z) == 0. )
        {
            // Complex infinity for all orders
            double dinf = std::numeric_limits<double>::infinity();
            for (int i=0; i<n; i++)
            {
                _cyl_h2[i] = std::complex<T2>(dinf, dinf);
            }
        }
        else if ( nu >= 0. )
        {
            // Array of only positive orders, including the order 0, with abs(z)!=0
            double * ch2r = new double [n];
            double * ch2i = new double [n];
            slatec::zbesh_(&x, &y, &fnu, &kode, &m, &n, ch2r, ch2i, &nz, &ierr);
            if (_flags) _flag_zbesh_(ierr, nz);
            for ( int i=0; i<n; i++ )
            {
                _cyl_h2[i] = std::complex<T2>(ch2r[i], ch2i[i]);
            }
        }
        else if ( nu_m <= 0. )
        {
            // Array of only negative orders, including the order 0, with abs(z)!=0
            double * ch2r_m = new double [n];
            double * ch2i_m = new double [n];
            double fnu_m = abs( nu_m );
            slatec::zbesh_(&x, &y, &fnu_m, &kode, &m, &n, ch2r_m, ch2i_m, &nz, &ierr);
            if (_flags) _flag_zbesh_(ierr, nz);
            for (int i=0; i<n; i++)
            {
                // Eq. (9.1.6) of Ref. [1].
                int tmp1 = n-1-i;
                double tmp2 = ( fnu_m + double(tmp1) )*M_PI;
                std::complex<T2> exp_of_mIpifnu = exp( std::complex<T2> (0.,-tmp2) );
                _cyl_h2[i] = std::complex<T2>(ch2r_m[tmp1], ch2i_m[tmp1])*exp_of_mIpifnu;
            }
        }
        else
        {
            // Array of negative and positive orders with abs(z)!=0
            int n_m = floor( abs(nu) ) + 1;
            int n_p = n - n_m;
            double fnu_m = abs( nu + double(n_m-1) );

            // The negative orders
            // int negative orders, or order 0, with abs(z)!=0
            double * ch2r_m = new double [n_m];
            double * ch2i_m = new double [n_m];
            slatec::zbesh_(&x, &y, &fnu_m, &kode, &m, &n_m, ch2r_m, ch2i_m, &nz, &ierr);
            if (_flags) _flag_zbesh_(ierr, nz);
            for (int i=0; i<n_m; i++)
            {
                // Eq. (9.1.6) of Ref. [1].
                int tmp1 = n_m-1-i;
                double tmp2 = ( fnu_m + double(tmp1) )*M_PI;
                std::complex<T2> exp_of_mIpifnu = exp( std::complex<T2> (0.,-tmp2) );
                _cyl_h2[i] = std::complex<T2>(ch2r_m[tmp1], ch2i_m[tmp1])*exp_of_mIpifnu;
            }

            // The positive orders
            double * ch2r_p = new double [n_p];
            double * ch2i_p = new double [n_p];
            double fnu_p = nu + double(n_m);
            slatec::zbesh_(&x, &y, &fnu_p, &kode, &m, &n_p, ch2r_p, ch2i_p, &nz, &ierr);
            if (_flags) _flag_zbesh_(ierr, nz);
            for (int i=n_m; i<n; i++)
            {
                int tmp1 = i-n_m;
                _cyl_h2[i] = std::complex<T2>(ch2r_p[tmp1], ch2i_p[tmp1]);
            }
        }
    }

    template<typename T1, typename T2>
    void cyl_h2( const T1 _nu, const int _n, T2 _z, T2 * _cyl_h2, bool _scaled = false, bool _flags = false )
    {
        std::complex<T2> * tmp_cylh2 = new std::complex<T2> [_n];
        cyl_h2( _nu, _n, std::complex<T2>(_z,0.), tmp_cylh2, _scaled, _flags );
        for (int i=0; i<_n; i++)
        {
            _cyl_h2[i] = real(tmp_cylh2[i]);
        }
    }

	// Modified cylindrical Bessel function of the first kind I_nu(z)
    template<typename T1, typename T2>
    std::complex<T2> cyl_i( const T1 _nu, const std::complex<T2> _z, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _nu:     real order.
        //  _z:      complex argument.
        //  _scaled: if true, returns the scaled version I_nu(z)*exp(-|Re(z)|).
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int nz;
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        std::complex<double> z = std::complex<double> (x, y);
        int n = 1;
        double nu = double(_nu);
        double fnu = abs(nu);

        if ( nu >= 0. )
        {
            // Positive order
            double cir, cii;
            slatec::zbesi_(&x, &y, &fnu, &kode, &n, &cir, &cii, &nz, &ierr);
            if (_flags) _flag_zbesi_(ierr, nz);
            return std::complex<T2>(cir, cii);
        }
        else
        {
            // Negative order
            if ( double(floor(fnu)) == fnu )
            {
                // int negative order, see Eq. (6.1.5) of [2].
                double cir, cii;
                slatec::zbesi_(&x, &y, &fnu, &kode, &n, &cir, &cii, &nz, &ierr);
                if (_flags) _flag_zbesi_(ierr, nz);
                return std::complex<T2>(cir, cii);
            }
            else if ( double(floor(fnu)) != fnu && abs(z) == 0. )
            {
                // Non-int negative order with abs(z)=0,
                //  I_nu = complex infinity since K_fnu = complex infinity.
                double dinf = std::numeric_limits<double>::infinity();
                return std::complex<double>(dinf, dinf);
            }
            else
            {
                // Non-int negative order with abs(z)!=0, see Eq. (6.5.4) of [2].
                double cir, cii;
                slatec::zbesi_(&x, &y, &fnu, &kode, &n, &cir, &cii, &nz, &ierr);
                if (_flags) _flag_zbesi_(ierr, nz);

                double ckr, cki;
                slatec::zbesk_(&x, &y, &fnu, &kode, &n, &ckr, &cki, &nz, &ierr);
                if (_flags) _flag_zbesk_(ierr, nz);

                //  Eq. (5.5.4) of Zhang and Jin (1996) [2].
                double tmp = fnu*M_PI;
                return std::complex<T2>(cir, cii) + std::complex<T2>(ckr, cki)*sin(tmp)*M_2_PI;
            }
        }
    }

    template<typename T1, typename T2>
    T2 cyl_i( const T1 _nu, const T2 _z, bool _scaled = false, bool _flags = false )
    {
        return real( cyl_i( _nu, std::complex<T2>(_z,0.), _scaled, _flags ) ) ;
    }

	// Array of modified cylindrical Bessel functions of the first kind I_nu(z)
    template<typename T1, typename T2>
    void cyl_i( const T1 _nu, const int _n, const std::complex<T2> _z, std::complex<T2> * _cyl_i, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _nu:     real first order to be calculated.
        //  _n:      int to return an array with the orders _nu, _nu+1, ..., _nu+_n.
        //  _z:      complex argument.
        //  _cyl_i:  array of size _n to output I_nu(z) for the orders _nu, _nu+1, ..., _nu+_n.
        //  _scaled: if true, returns an array of size _n with the scaled version I_nu(z)*exp(-|Re(z)|) for the orders _nu, _nu+1, ..., _nu+_n.
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int nz;
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        std::complex<double> z = std::complex<double> (x, y);
        int n = _n;
        double nu = double(_nu);
        double fnu = abs(nu);
        double nu_m = nu + double(n-1);

        if ( n <= 0 )
        {
            if (_flags) std::cerr << "[bessel::cyl_i] Flag -> Input error: No computation, _n must be an integer greater than zero." << std::endl;
        }
        else if ( nu >= 0. )
        {
            // Array of only positive orders, including the order 0
            double * cir = new double [n];
            double * cii = new double [n];
            slatec::zbesi_(&x, &y, &fnu, &kode, &n, cir, cii, &nz, &ierr);
            if (_flags) _flag_zbesi_(ierr, nz);

            for ( int i=0; i<n; i++ )
            {
                _cyl_i[i] = std::complex<T2>(cir[i], cii[i]);
            }
        }
        else if ( nu_m <= 0. )
        {
            // Array of only negative orders, including the order 0
            if ( double(floor(fnu)) == fnu )
            {
                // int negative orders, or order 0, see Eq. (6.1.5) of [2].
                double * cir = new double [n];
                double * cii = new double [n];
                double fnu_m = abs( nu_m );
                slatec::zbesi_(&x, &y, &fnu_m, &kode, &n, cir, cii, &nz, &ierr);
                if (_flags) _flag_zbesi_(ierr, nz);

                for (int i=0; i<n; i++)
                {
                    int tmp = n-1-i;
                    _cyl_i[i] = std::complex<T2>(cir[tmp], cii[tmp]);
                }
            }
            else if ( double(floor(nu_m)) != nu_m && abs(z) == 0. )
            {
                // Non-int negative order with abs(z)=0,
                //  I_nu = complex infinity since K_fnu = complex infinity.
                double dinf = std::numeric_limits<double>::infinity();
                for (int i=0; i<n; i++)
                {
                    _cyl_i[i] = std::complex<T2>(dinf, dinf);
                }
            }
            else
            {
                // Non-int negative orders with abs(z)!=0, see Eq. (6.5.4) of [2].
                double * cir = new double [n];
                double * cii = new double [n];
                double fnu_m = abs( nu_m );
                slatec::zbesi_(&x, &y, &fnu_m, &kode, &n, cir, cii, &nz, &ierr);
                if (_flags) _flag_zbesi_(ierr, nz);

                double * ckr = new double [n];
                double * cki = new double [n];
                slatec::zbesk_(&x, &y, &fnu_m, &kode, &n, ckr, cki, &nz, &ierr);
                if (_flags) _flag_zbesk_(ierr, nz);

                for (int i=0; i<n; i++)
                {
                    //  Eq. (6.5.4) of Zhang and Jin (1996) [2].
                    int tmp1 = n-1-i;
                    double tmp2 = ( fnu_m + double(tmp1) )*M_PI;
                    _cyl_i[i] = ( std::complex<T2>(cir[tmp1], cii[tmp1]) ) + ( std::complex<T2>(ckr[tmp1], cki[tmp1]) )*sin(tmp2)*M_2_PI;
                }
            }
        }
        else
        {
            // Array of negative and positive orders
            int n_m = floor( abs(nu) ) + 1;
            int n_p = n - n_m;
            double fnu_m = abs( nu + double(n_m-1) );

            // The negative orders
            if ( double(floor(fnu)) == fnu )
            {
                // int negative orders, or order 0, see Eq. (6.1.5) of [2].
                double * cir_m = new double [n_m];
                double * cii_m = new double [n_m];
                slatec::zbesi_(&x, &y, &fnu_m, &kode, &n_m, cir_m, cii_m, &nz, &ierr);
                if (_flags) _flag_zbesi_(ierr, nz);

                for (int i=0; i<n_m; i++)
                {
                    int tmp = n_m-1-i;
                    _cyl_i[i] = std::complex<T2>(cir_m[tmp], cii_m[tmp]);
                }
            }
            else if ( double(floor(nu_m)) != nu_m && abs(z) == 0. )
            {
                // Non-int negative orders with abs(z)=0,
                //  I_nu = complex infinity since K_fnu = complex infinity.
                double dinf = std::numeric_limits<double>::infinity();
                for (int i=0; i<n_m; i++)
                {
                    // Complex infinity
                    _cyl_i[i] = std::complex<T2>(dinf, dinf);
                }
            }
            else
            {
                // Non-int negative orders with abs(z)!=0,
                double * cir_m = new double [n_m];
                double * cii_m = new double [n_m];
                slatec::zbesi_(&x, &y, &fnu_m, &kode, &n_m, cir_m, cii_m, &nz, &ierr);
                if (_flags) _flag_zbesi_(ierr, nz);

                double * ckr_m = new double [n_m];
                double * cki_m = new double [n_m];
                slatec::zbesk_(&x, &y, &fnu_m, &kode, &n_m, ckr_m, cki_m, &nz, &ierr);
                if (_flags) _flag_zbesk_(ierr, nz);

                for (int i=0; i<n_m; i++)
                {
                    //  Eq. (6.5.4) of Zhang and Jin (1996) [2].
                    int tmp1 = n_m-1-i;
                    double tmp2 = ( fnu_m + double(tmp1) )*M_PI;
                    _cyl_i[i] = ( std::complex<T2>(cir_m[tmp1], cii_m[tmp1]) ) + ( std::complex<T2>(ckr_m[tmp1], cki_m[tmp1]) )*sin(tmp2)*M_2_PI;
                }
            }

            // The positive orders
            double * cir_p = new double [n_p];
            double * cii_p = new double [n_p];
            double fnu_p = nu + double(n_m);
            slatec::zbesi_(&x, &y, &fnu_p, &kode, &n_p, cir_p, cii_p, &nz, &ierr);
            if (_flags) _flag_zbesi_(ierr, nz);

            for (int i=n_m; i<n; i++)
            {
                int tmp1 = i-n_m;
                _cyl_i[i] = std::complex<T2>(cir_p[tmp1], cii_p[tmp1]);
            }
        }
    }

    template<typename T1, typename T2>
    void cyl_i( const T1 _nu, const int _n, T2 _z, T2 * _cyl_i, bool _scaled = false, bool _flags = false )
    {
        std::complex<T2> * tmp_modi = new std::complex<T2> [_n];
        cyl_i( _nu, _n, std::complex<T2>(_z,0.), tmp_modi, _scaled, _flags );
        for (int i=0; i<_n; i++)
        {
            _cyl_i[i] = real(tmp_modi[i]);
        }
    }

	// Modified cylindrical Bessel function of the second kind K_nu(z) [Basset or MacDonald function]
    template<typename T1, typename T2>
    std::complex<T2> cyl_k( const T1 _nu, const std::complex<T2> _z, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _nu:     real order.
        //  _z:      complex argument.
        //  _scaled: if true, returns the scaled version K_nu(z)*exp(z).
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int nz;
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        std::complex<double> z = std::complex<double> (x, y);
        int n = 1;
        double nu = double(_nu);
        double fnu = abs(nu);

        if ( abs(z) == 0. )
        {
            double dinf = std::numeric_limits<double>::infinity();
            if ( nu == 0. )
            {
                // Infinity
                return std::complex<double>(dinf, 0.);
            }
            else
            {
                // Complex infinity
                return std::complex<double>(dinf, dinf);
            }
        }
        else
        {
            // Positive or negative order, see Eq. (6.5.5) of Ref. [2].
            double ckr, cki;
            slatec::zbesk_(&x, &y, &fnu, &kode, &n, &ckr, &cki, &nz, &ierr);
            if (_flags) _flag_zbesk_(ierr, nz);
            return std::complex<T2>(ckr, cki);
        }
    }

    template<typename T1, typename T2>
    T2 cyl_k( const T1 _nu, const T2 _z, bool _scaled = false, bool _flags = false )
    {
        return real( cyl_k( _nu, std::complex<T2>(_z,0.), _scaled, _flags ) ) ;
    }

	// Array of modified cylindrical Bessel functions of the second kind K_nu(z) [Basset or MacDonald functions]
    template<typename T1, typename T2>
    void cyl_k( const T1 _nu, const int _n, const std::complex<T2> _z, std::complex<T2> * _cyl_k, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _nu:     real first order to be calculated.
        //  _n:      int to return an array with the orders _nu, _nu+1, ..., _nu+_n.
        //  _z:      complex argument.
        //  _cyl_k:  array of size _n to output K_nu(z) for the orders _nu, _nu+1, ..., _nu+_n.
        //  _scaled: if true, returns an array of size _n with the scaled version K_nu(z)*exp(z) for the orders _nu, _nu+1, ..., _nu+_n.
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int nz;
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        std::complex<double> z = std::complex<double> (x, y);
        int n = _n;
        double nu = double(_nu);
        double fnu = abs(nu);
        double nu_m = nu + double(n-1);

        if ( n <= 0 )
        {
            if (_flags) std::cerr << "[bessel::cyl_k] Flag -> Input error: No computation, _n must be an integer greater than zero." << std::endl;
        }
        else if ( abs(z) == 0. )
        {
            double dinf = std::numeric_limits<double>::infinity();
            for (int i=0; i<n; i++)
            {
                // Complex infinity for all orders
                _cyl_k[i] = std::complex<T2>(dinf, dinf);
            }
            if ( nu_m == floor(nu_m) && n > fnu )
            {
                // Infinity for the int order 0
                _cyl_k[ int(floor(fnu)) ] = std::complex<double>(dinf, 0.);
            }
        }
        else if ( nu >= 0. )
        {
            // Array of only positive orders, including the order 0, with abs(z)!=0
            double * ckr = new double [n];
            double * cki = new double [n];
            slatec::zbesk_(&x, &y, &fnu, &kode, &n, ckr, cki, &nz, &ierr);
            if (_flags) _flag_zbesk_(ierr, nz);

            for ( int i=0; i<n; i++ )
            {
                _cyl_k[i] = std::complex<T2>(ckr[i], cki[i]);
            }
        }
        else if ( nu_m <= 0. )
        {
            // Array of only negative orders, including the order 0, with abs(z)!=0
            // See Eq. (6.5.5) of Ref. [2].
            double * ckr = new double [n];
            double * cki = new double [n];
            double fnu_m = abs( nu_m );
            slatec::zbesk_(&x, &y, &fnu_m, &kode, &n, ckr, cki, &nz, &ierr);
            if (_flags) _flag_zbesk_(ierr, nz);

            for (int i=0; i<n; i++)
            {
                int tmp = n-1-i;
                _cyl_k[i] = std::complex<T2>(ckr[tmp], cki[tmp]);
            }
        }
        else
        {
            // Array of negative and positive orders with abs(z)!=0
            int n_m = floor( abs(nu) ) + 1;
            int n_p = n - n_m;
            double fnu_m = abs( nu + double(n_m-1) );

            // The negative orders
            // int negative orders, or order 0, with abs(z)!=0
            // See Eq. (6.5.5) of Ref. [2].
            double * ckr_m = new double [n_m];
            double * cki_m = new double [n_m];
            slatec::zbesk_(&x, &y, &fnu_m, &kode, &n_m, ckr_m, cki_m, &nz, &ierr);
            if (_flags) _flag_zbesk_(ierr, nz);

            for (int i=0; i<n_m; i++)
            {
                int tmp = n_m-1-i;
                _cyl_k[i] = std::complex<T2>(ckr_m[tmp], cki_m[tmp]);
            }

            // The positive orders
            double * ckr_p = new double [n_p];
            double * cki_p = new double [n_p];
            double fnu_p = nu + double(n_m);
            slatec::zbesk_(&x, &y, &fnu_p, &kode, &n_p, ckr_p, cki_p, &nz, &ierr);
            if (_flags) _flag_zbesk_(ierr, nz);

            for (int i=n_m; i<n; i++)
            {
                int tmp1 = i-n_m;
                _cyl_k[i] = std::complex<T2>(ckr_p[tmp1], cki_p[tmp1]);
            }
        }
    }

    template<typename T1, typename T2>
    void cyl_k( const T1 _nu, const int _n, T2 _z, T2 * _cyl_k, bool _scaled = false, bool _flags = false )
    {
        std::complex<T2> * tmp_modk = new std::complex<T2> [_n];
        cyl_k( _nu, _n, std::complex<T2>(_z,0.), tmp_modk, _scaled, _flags );
        for (int i=0; i<_n; i++)
        {
            _cyl_k[i] = real(tmp_modk[i]);
        }
    }

	// Airy Ai function
    template<typename T1>
    std::complex<T1> airy_ai( const std::complex<T1> _z, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _z:      complex argument.
        //  _scaled: if true, returns the scaled version exp[(2/3)*z*sqrt(z)]*Ai(z).
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int nz;
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        int id = 0;
        double air, aii;
        slatec::zairy_(&x, &y, &id, &kode, &air, &aii, &nz, &ierr);
        if (_flags) _flag_zairy_(ierr, nz);
        return std::complex<T1>(air, aii);
    }

    template<typename T1>
    T1 airy_ai( const T1 _z, bool _scaled = false, bool _flags = false )
    {
        return real( airy_ai( std::complex<T1>(_z,0.), _scaled, _flags ) ) ;
    }

	// Airy Bi function
    template<typename T1>
    std::complex<T1> airy_bi( const std::complex<T1> _z, bool _scaled = false, bool _flags = false )
    {
        // Parameters
        //  _z:      complex argument.
        //  _scaled: if true, returns the scaled version exp{-|Re[(2/3)*z*sqrt(z)]|}*Bi(z).
        //  _flags:  if true, print error messages.

        int kode = ( _scaled ? 2 : 1 );
        int ierr;
        double x = double( real(_z) );
        double y = double( imag(_z) );
        int id = 0;
        double bir, bii;
        slatec::zbiry_(&x, &y, &id, &kode, &bir, &bii, &ierr);
        if (_flags) _flag_zbiry_(ierr);
        return std::complex<T1>(bir, bii);
    }

    template<typename T1>
    T1 airy_bi( const T1 _z, bool _scaled = false, bool _flags = false )
    {
        return real( airy_bi( std::complex<T1>(_z,0.), _scaled, _flags ) ) ;
    }

	// C++ CALLABLE BESSEL ROUTINES
	//-^------------------------------------------------------------------------
} // namespace bessel

// C++ CODE
//-^------------------------------------------------------------------------
