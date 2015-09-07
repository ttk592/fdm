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

#include "pde1d.h"

#include <cstdio>
#ifdef USE_FENV
#  include <fenv.h>
#endif

#include <fdm/common.hpp>
#include <fdm/grid.hpp>
#include <fdm/functors.hpp>
#include <fdm/solve_pde.hpp>


// this programme solves 1d const coefficient parabolic pde's
//  u_t = a u_xx + b u_x + c u ,   x in [-inf, inf]
// with the initial condition of
//  u(x,0)=(e^x-K)^+
// we know the solution explicitly and can compare against the
// numerical results
// this is also an ideal test place for comparing different
// non-uniform grids


// Fundamental solution to u_t = a u_xx + b u_x + c u
double	G(double x, double y, double t, double a, double b, double c)
{
    return exp(c*t - (b*t+x-y)*(b*t+x-y)/(4.0*a*t)) / sqrt(4.0*M_PI*a*t);
}

// solution of u_t = a u_xx + b u_x + c u
// with u(x,0)=(e^x-K)^+
double U(double x, double t, double a, double b, double c, double K)
{
    double	value;
    value = K*exp(c*t) * (erf( (log(K)-x-b*t)/(2.0*sqrt(a*t)) ) - 1.0);
    value-= exp(x+(a+b+c)*t) *
            (erf( (log(K)-x-(2.0*a+b)*t)/(2.0*sqrt(a*t)) ) - 1.0);
    return 0.5*value;
}




int main(int argc, char** argv)
{

    if(argc<1) {
        printf("usage: %s <>\n", argv[0]);
        exit(EXIT_FAILURE);
    }


#ifdef USE_FENV
    // programme will terminate if any of the following errors occur
    // so we catch nan, inf, etc
    //feenableexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
    feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
#endif

    boost::timer	time_all;

    // reading pde parameters,  u_t = diffusion u_xx - convection u_x + c u
    double	T=1.0;
    double	diffusion=1.0;
    double	convection=0.5;
    double	c=0.0;
    double	K=1.0;			// strike for initial condition

    // defining lower and upper boundaries
    fdm::SVector<double,1>	x;	// point of interest
    x[0]=0.0;
    double	lower=-5.0;
    double	upper=5.0;

    // reading calculation parameters
    int	n_time=50;
    int	n_x=200;


    // generating grids
    fdm::Grid<1>		grid_time;
    fdm::GridC2<1>	grid_space;

    grid_time.set_uniform(0,0.0,T,n_time+1);
    grid_space.set_uniform(0,lower,upper,n_x);
    grid_space.calculate();


    // defining the pde
    fdm::PDEConst1D		pde;

    pde.a = diffusion;
    pde.b = -convection;
    pde.c = c;

    // defining the boundary condition
    // make sure we prescribe fixed values at boundaries with incoming
    // convection, and free boundary condition for out flowing boundaries
    fdm::BoundaryConst<1>	boundary;

    boundary.TypeLower[0]=fdm::BoundaryType::Value;
    boundary.ValueLower[0]=0.0;
    boundary.TypeUpper[0]=fdm::BoundaryType::Free;
    boundary.ValueUpper[0]=0.0;


    // setting initial conditions
    fdm::GridFunction<1>    u(grid_space);
    u.iterator.begin();
    double	z;
    while(u.iterator.isend()==false) {
        z=u.coord()[0];
        u.value()=std::max(exp(z)-K,0.0);
        u.iterator++;
    }


    // solving the pde
    double	t,exact,value;
    for(size_t k=1; k<grid_time.size(0); k++) {
        fdm::FDMSolveStep(pde,boundary,u,grid_time,x,k,fdm::type::cn_1d,false,false);
        // check error
        t=grid_time.coord(0,k);
        exact=U(x[0],t, pde.a, pde.b, pde.c, K);
        value=u(x);
        printf("t=%.3f: num = %f, exact = %f   (error=%f %%)\n",
               t,value,exact, (value-exact)/exact*100.0 );
        /*
           for(int i=0;i<n_x;i++){
              z=u.coord(i);
              u.value(i);
           }
        */
    }
    value=u(x);

    printf("-------------------- result --------------------\n");

    printf("Numerical value:\t %12.10f\n",value);
    printf("elapsed time: %.1f s:\n", time_all.elapsed());


    return 0;		// success
}
