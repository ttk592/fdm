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


#ifndef MF_FDM_PDE_HPP
#define MF_FDM_PDE_HPP

#include <cassert>
#include <fdm/svector.hpp>


namespace fdm
{


// Represenation of the coefficients of parabolic PDEs
//
// u_t = sum_ij deriv2(i,j) d^2/(dx_i dx_j) u
//     + sum_i  deriv1(i)   d/dx_i u
//     + deriv0() u + f
//
// by convention we only use deriv2(i,j) with j>=i, e.g. in 2d
// deriv2(1,0) u_{yx} will be zero and all the coefficients go into
// deriv2(0,1) u_{xy}
//
// allows for auxilliary variables, i.e. coefficients can depend on other
// variables but PDE itself has no derivatives with respect to that variable,
// e.g. u_t = a(y) u_xx
//
// disadvantage of this representation: we'll have to use if(dim==0) etc
//      statements in the derived classes to obtain particular coefficients,
//	due to the virtual functions calls, and the use via base class
//	reference, FDMSolveStep(PDE<>&, ...), the compiler won't simplify
//	(not even with -funroll-loops -finline-functions)
//	TODO: bench to see whether this is a performance issue
//
template <std::size_t GridDim, std::size_t AuxDim=0>
class PDE
{
protected:
    SVector<double, GridDim>	m_x;	// space coordinate x
    SVector<double, AuxDim>	m_y;	// auxilliary coordinate y
    double			m_t;	// time t
public:
    int dim() const
    {
        return GridDim;
    }
    void set_space(const SVector<double, GridDim>& x)
    {
        m_x=x;
    }
    void set_time(double t)
    {
        m_t=t;
    }
    void set_aux(const SVector<double, AuxDim>& y)
    {
        m_y=y;
    }
    virtual double deriv2(int dim1, int dim2) const = 0;
    virtual double deriv1(int dim) const = 0;
    virtual double deriv0() const = 0;
    virtual double f() const = 0;       // TODO all solvers currently ignore this
    virtual bool is_time_indep() const = 0;	// if all coeff's time indep
};





// not sure - could that be more efficient?
// advantage: no need for: if(dim==0) etc statements to get to the coefficients
// disadvantage: have to calculate all the coeff, eventhough we might
//	only need a few
/*
template <std::size_t GridDim>
class PDE {
 protected:
   boost::array< boost::array<double, GridDim>, GridDim >	wDeriv2;
   boost::array< double, GridDim >				wDeriv1;
   double							wDeriv0;
   double							wF;
 public:
   const int dim() const { return GridDim; }
   // sets all coefficients wDeriv based on x and t
   virtual void set(const SVector<double, GridDim>& x, double t) = 0;
   virtual bool is_time_indep() const = 0;
   inline double deriv2(int dim1, int dim2) const {
	assert( (0<=dim1) && (dim1<this->dim()) );
	assert( (0<=dim2) && (dim2<this->dim()) );
	return wDeriv2[dim1][dim2];
   }
   inline double deriv1(int dim) const {
	assert( (0<=dim) && (dim<this->dim()) );
	return wDeriv1[dim];
   }
   inline double deriv0() const {
	return wDeriv0;
   }
   inline double f() const {
	return wF;
   }
};

*/


} // namespace fdm




#endif /* MF_FDM_PDE_HPP  */
