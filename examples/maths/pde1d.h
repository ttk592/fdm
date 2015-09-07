/*
 * ---------------------------------------------------------------------
 * Copyright (C) 2009, 2010, 2015 Tino Kluge (ttk448 at gmail.com)
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


#ifndef MF_FDM_PDE1D_H
#define MF_FDM_PDE1D_H

#include <fdm/common.hpp>
#include <fdm/solve_pde.hpp>

namespace fdm
{

// constant coefficient pde in one dimension
// u_t = a u_xx + b u_x + c u + f
class PDEConst1D : public PDE<1>
{
public:
    double	a,b,c;		// coeff of the pde

    inline double deriv2(int UNUSED(dim1), int UNUSED(dim2) ) const
    {
        assert( (dim1==0) && (dim2==0) );
        return this->a;
    }
    inline double deriv1(int UNUSED(dim) ) const
    {
        assert( dim==0 );
        return this->b;
    }
    inline double deriv0() const
    {
        return this->c;
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

#endif // MF_FDM_PDE1D_H
