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

#include "heat.h"

#include <fdm/grid.hpp>
#include <fdm/solve_pde.hpp>
#include "png_io.h"


int main(int argc, char** argv)
{

    if(argc<1) {
        printf("usage: %s <>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    fdm::enableexcept();    // for testing, raise fpu exceptions


    boost::timer	time_all;

    // reading pde parameters
    double	T=20.0;
    double	diffusion=0.5;
    fdm::SVector<double,2>	convection;
    fdm::SVector<double,2>	x;
    convection[0]=1.0;
    convection[1]=0.5;

    x[0]=0.0;
    x[1]=0.0;


    // reading calculation parameters
    int	n_time=(int) (T*25);
    int	n_x=720;
    int	n_y=480;

    // n_x=360;
    // n_y=240;


    // generating grids
    fdm::Grid<1>	grid_time;
    fdm::GridC2<2>	grid_space;

    grid_time.set_uniform(0,0.0,T,n_time+1);
    grid_space.set_uniform(0,-7.2,7.2,n_x);
    grid_space.set_uniform(1,0.-4.8,4.8,n_y);
    grid_space.calculate();


    // defining the pde
    fdm::PDEHeat2D		pde;

    pde.Diffusion[0][0] = diffusion;
    pde.Diffusion[1][1] = diffusion;
    pde.Diffusion[0][1] = 0.0;
    pde.Diffusion[1][0] = 0.0;

    pde.Convection[0] = convection[0];
    pde.Convection[1] = convection[1];

    pde.Const = 0.0;

    // defining the boundary condition
    // make sure we prescribe fixed values at boundaries with incoming
    // convection, and free boundary condition for out flowing boundaries
    fdm::BoundaryConst<2>	boundary;

    boundary.TypeLower[0]=fdm::BoundaryType::Value;
    boundary.ValueLower[0]=0.0;
    boundary.TypeUpper[0]=fdm::BoundaryType::Free;
    boundary.ValueUpper[0]=0.0;

    boundary.TypeLower[1]=fdm::BoundaryType::Value;
    boundary.ValueLower[1]=0.0;
    boundary.TypeUpper[1]=fdm::BoundaryType::Free;
    boundary.ValueUpper[1]=0.0;


    // setting initial conditions
    fdm::GridFunction<2>    u(grid_space);
    u.iterator.begin();
    while(u.iterator.isend()==false) {
        if( u.coord().norm() <3.0 ) {
            u.value()=1.0;
        } else {
            u.value()=0.0;
        }
        u.iterator++;
    }


    // solving the pde
    boost::multi_array<double, 2>	pixels(boost::extents[n_y][n_x]);
    fdm::SVector<int,2>	idx;
    char	filename[30];
    for(size_t k=1; k<grid_time.size(0); k++) {
        fdm::FDMSolveStep(pde,boundary,u,grid_time,x,k,fdm::type::pred_corr,true,false);
        // saving the result in a 2d array
        for(int i=0; i<n_y; i++) {
            for(int j=0; j<n_x; j++) {
                idx[0]=j;
                idx[1]=i;
                pixels[i][j]=u.value(idx);
            }
        }
        // output of the result in a png in form of a heat map
        sprintf(filename,"out/heatmap-%.6lu.png",k);
        png::write_png(pixels,filename);
    }
    double	value=u(x);

    printf("-------------------- result --------------------\n");

    printf("Numerical value:\t %12.10f\n",value);
    printf("elapsed time: %.1f s:\n", time_all.elapsed());

    printf("\n");
    printf("generate an animation by calling\n\n");
    printf("mencoder \"mf://out/*.png\" -mf fps=25 -o out.avi -ovc lavc -lavcopts vcodec=mpeg4\n\n");

    return 0;		// success
}
