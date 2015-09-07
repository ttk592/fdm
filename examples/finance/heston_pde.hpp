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


#ifndef MF_FDM_HESTON_PDE_HPP
#define MF_FDM_HESTON_PDE_HPP


#include <fdm/pde.hpp>
#include "market.hpp"

namespace fin
{

// log spot Heston PDE
// -------------------
// dim 0 = spot dimension (logarithmic scale)
// dim 1 = vol dimension (variance)
class PDEHeston : public fdm::PDE<2>
{
private:
    MarketHeston    m_market;
public:
    static const int dim_spot=0;
    static const int dim_vol=1;

    void import_params(const MarketHeston& heston)
    {
        m_market=heston;
    }
    double deriv2(int dim1,int dim2) const
    {
        assert( (0<=dim1) && (dim1<this->dim()) );
        assert( (0<=dim2) && (dim2<this->dim()) );
        if(dim1==0) {
            if(dim2==0) {
                // a_xx=0.5*y;
                return 0.5*this->m_x[1];
            } else {
                // a_xy=rho*xi*y;
                return m_market.rho*m_market.xi*m_x[1];
            }
        } else if( (dim1==1) && (dim2==1) ) {
            // a_yy(i,j)=0.5*sqr(xi)*y;
            return 0.5*m_market.xi*m_market.xi*m_x[1];
        }
        return 0.0;
    }
    inline double deriv1(int dim) const
    {
        assert( (0<=dim) && (dim<this->dim()) );
        if(dim==0) {
            // a_x=(rd-rf-0.5*y)
            return (m_market.rd-m_market.rf-0.5*m_x[1]);
        } else {
            // a_y=kappa*(theta-y)
            return m_market.kappa*(m_market.theta-m_x[1]);
        }

    }
    inline double deriv0() const
    {
        return -m_market.rd;
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



} // namespace fin

#endif // MF_FDM_HESTON_PDE_HPP

