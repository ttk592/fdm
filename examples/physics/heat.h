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


#ifndef MF_FDM_HEAT_H
#define MF_FDM_HEAT_H


#include <fdm/solve_pde.hpp>

namespace fdm
{

// convection diffusion equation (constant coefficients)
// u_t = div(A grad u) - div(u b) + c u + f
// u_t = A u_xx - b u_x + c u + f (because of const coeff)
class PDEHeat2D : public PDE<2>
{
public:
    boost::array<boost::array<double, 2>, 2>	Diffusion;	// A
    boost::array<double, 2>			Convection;	// b
    double					Const;		// c

    inline double deriv2(int dim1,int dim2) const
    {
        assert( (0<=dim1) && (dim1<this->dim()) );
        assert( (0<=dim2) && (dim2<this->dim()) );
        return this->Diffusion[dim1][dim2];
    }
    inline double deriv1(int dim) const
    {
        assert( (0<=dim) && (dim<this->dim()) );
        return - this->Convection[dim];
    }
    inline double deriv0() const
    {
        return this->Const;
    }
    inline double f() const
    {
        return 0.0;
    }
    inline bool is_time_indep() const
    {
        return true;
    }
};

} // namespace fdm

#endif // MF_FDM_HEAT_H

