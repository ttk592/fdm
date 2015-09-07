/*
 * ---------------------------------------------------------------------
 * Copyright (C) 2002, 2009, 2010, 2015 Tino Kluge (ttk448 at gmail.com)
 * Copyright (C) 2009, 2010 MathFinance AG (info@mathfinance.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see
 *  <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 */


#ifndef MF_FDM_FUNCTORS_H
#define MF_FDM_FUNCTORS_H


// to get access to M_PI and other constants in Visual Studio we need
// to define _USE_MATH_DEFINES before inclusion of <cmath>
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>

// for any compiler not C99 compliant certain math functions
// are unavailable (eg erf, asinh), include boost to provide those
#ifndef __GNUC__
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/asinh.hpp>
using boost::math::erf;
using boost::math::asinh;
#endif



namespace fdm
{

// base class for one dimensional function
class FunctionREAL
{
protected:
public:
    virtual ~FunctionREAL() { };
    virtual double operator () (double x)  const = 0;

    // integral of 1/f(x)
    double integrate_simpson_inv(double a, double b, int n_points=10000)
    {
        assert( (b>a) && (n_points>2) );
        n_points += n_points%2;		// we need an even number of points

        double	sum=0.0;
        double	dx=(b-a)/n_points;
        double	x=a+0.5*dx;
        while(x<b) {
            sum += 1.0/this->operator()(x);
            sum += 4.0/this->operator()(x+dx);
            sum += 1.0/this->operator()(x+2.0*dx);
            x   += 2.0*dx;
        }
        return sum*2.0*dx/6.0;
    }
    double integrate_simple_inv(double a, double b, int n_points=10000)
    {
        assert( (b>a) && (n_points>2) );
        double	sum=0.0;
        double	dx=(b-a)/n_points;
        double	x=a+0.5*dx;
        while(x<b) {
            sum += 1.0/this->operator()(x);
            x   += dx;
        }
        return sum*dx;
    }

};


// numerical solution of an ode x' = g(x) in the range [a,b] with x(a)=x_0
std::vector<double> solve_ode_runge(const FunctionREAL& g,
                                    double a, double b, int n_points, double x_0);

// numerical solution of an ode x' = g(x) in the range [a,b] with x(a)=x_0
std::vector<double> solve_ode_euler(const FunctionREAL& g,
                                    double a, double b, int n_points, double x_0);


// function class, which gives a messure for the distance between two
// adjacent grid points: g(x) says how dense grid points are compared
// to a uniform net (g(x)=1.0 means same as uniform, smaller = denser)
// more abstract
// g:[a,b] --> R+  which satisfies  int_a^b 1/g(x) dx = b-a
class FunctionDistance : public FunctionREAL
{
protected:
    // the intention is that operator() contains wScaleFactor as a parameter
    // so that it can be set appropriately to satisfy
    // int_a^b 1/operator()(x) dx = b-a
    double	wScaleFactor;
    double	wA, wB;

public:
    FunctionDistance()
    {
        wScaleFactor=1.0;
    };
    FunctionDistance(double a, double b)
    {
        assert(a<b);
        wScaleFactor=1.0;
        wA=a;
        wB=b;
    };
    void set_bounds(double a, double b, double scale=-1.0)
    {
        wA=a;
        wB=b;
        if(scale>0.0) {
            wScaleFactor=scale;
        }
    }
    virtual ~FunctionDistance() { };

    // if wScaleFactor is linear in operator() then we can easily
    // asjust it to satisfy the integral constraint
    bool adapt_scale_factor(double a, double b, int n_points=10000)
    {
        wScaleFactor=1.0;
        double	value=this->integrate_simpson_inv(a,b,n_points);
        wScaleFactor=1.0/value;
        assert( fabs(this->integrate_simpson_inv(a,b,n_points)-(b-a))
                < 1e-15*(b-a) );
        return true;
    }
    // if wScaleFactor is non-linear in operator() then we have to use
    // an iterative method, very experimental, might fail!
    bool adapt_scale_factor_nonlin(double a, double b, int n_points=10000);

};


// apparently optimal mesh for u_t=wDiffusion u_xx
// tk diploma thesis version
class FunctionDistanceDiffusionOptim : public FunctionDistance
{
public:
    double	wDiffusion;	// strength of diffusion
    double	wTime;		// time to maturity
    double	wCentre;	// coordinate of the centre

    double operator () (double x) const
    {
        double	dist = fabs(x-wCentre);
        double	c_1  = 1.0/sqrt(4.0*M_PI*wDiffusion);
        double	c_2  = dist*dist/(4.0*wDiffusion);
        double	value_centre=1.0/sqrt(2.0*c_1*(sqrt(wTime)));
        double	value=2.0*c_1*(
                          sqrt(wTime)*exp(-c_2/wTime) +
                          sqrt(c_2*M_PI)*erf(sqrt(c_2/wTime)) - sqrt(c_2*M_PI)
                      );
        if(value<=0.0) value=1E3*value_centre;
        else value=1.0/sqrt(value);

        return wScaleFactor*std::min(value,1E3*value_centre);
    }
};

// apparently optimal mesh for u_t= wDiffusion (x u_xx +1/2 u_x)
// tk diploma thesis version
class FunctionDistanceLinDiffusionOptim : public FunctionDistance
{
public:
    double	wDiffusion;		// strength of diffusion
    double	wTime;			// time to maturity
    double	wCentre;		// coordinate of the centre

    double operator () (double x) const
    {
        if(x<=0.0) x=1E-5;
        double	dist = fabs(sqrt(x)-sqrt(wCentre));
        double	c_1  = 1/sqrt(4*M_PI*wDiffusion);
        double	c_2  = dist*dist/wDiffusion;
        double	value_centre=2*c_1*( sqrt(wTime));
        value_centre=wScaleFactor/sqrt(value_centre);
        double	value=2*c_1*(
                          sqrt(wTime)*exp(-c_2/wTime) + sqrt(c_2*M_PI)*erf(sqrt(c_2/wTime)) - sqrt(
                              c_2*M_PI)
                      );
        if(value==0) value=1E3*value_centre;
        else value=wScaleFactor/sqrt(value);
        return wScaleFactor*std::min(value,1E3*value_centre);
    }
};


// points are concentrated arround a centre
// do not use in production environment, use the explicit solution
// to this, see FunctionMapping1
class FunctionDistance1 : public FunctionDistance
{
public:
    double	Beta;		// how dense are points at the centre
    double	Centre;		// coordinate of the centre
    double operator () (double x) const
    {
        return sqrt(Beta*Beta+wScaleFactor*(x-Centre)*(x-Centre));
    }
};



// function class for a monotone bijective mapping
// f:[0,1]-->[a,b]
class FunctionMapping : public FunctionREAL
{
protected:
    double	wA;
    double	wB;
    // optional parameter in case f(1)=b is not automatically satisfied
    double	wParam;
public:
    FunctionMapping()
    {
        wParam=1.0;
    }
    FunctionMapping(double a, double b)
    {
        wParam=1.0;
        wA=a;
        wB=b;
    }
    double lower() const
    {
        return wA;
    }
    double upper() const
    {
        return wB;
    }
    // uses newton iteration to modify wParam to satisfy f(b)=b
    bool adapt_parameter();
};


// explicit solution to FunctionDistance1
class FunctionMapping1 : public FunctionMapping
{
protected:
    double	wBeta;		// how dense are points at the centre
    double	wCentre;	// coordinate of the centre in [0,1]
public:
    void set_param(double a, double b, double centre, double beta)
    {
        bool	status;
        assert( (a<b) && (beta>0.0) );
        wA=a;
        wB=b;
        wBeta=beta;
        wCentre=(centre-a)/(b-a);
        status=this->adapt_parameter();
        if(status==false) {
            printf("FunctionMapping1 failed: f(1)=%f, b=%f, error=%e\n",
                   this->operator()(1.0), wB, this->operator()(1.0)-wB);
            abort();		// replace this with some error handling?
        }
        assert( fabs(this->operator()(1.0)-wB) < fabs(wB)*1e-14 );
    }
    double operator() (double x) const
    {
        return wA + (wB-wA) * (
                   wCentre + wBeta/wParam * sinh(wParam*x+asinh(-wParam/wBeta*wCentre))
               );
    }
};


} // namespace fdm

#ifndef USING_C99_CMATH
#undef _USE_MATH_DEFINES
#endif


// implementation part (could be moved to cpp file)
#include <fdm/functors.inl>


#endif // MF_FDM_FUNCTORS_H
