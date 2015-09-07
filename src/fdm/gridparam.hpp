/*
 * ---------------------------------------------------------------------
 * Copyright (C) 2009, 2010 MathFinance AG (info@mathfinance.com)
 * Copyright (C) 2009, 2010, 2015 Tino Kluge (ttk448 at gmail.com)
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


#ifndef MF_FDM_GRIDPARAM_H
#define MF_FDM_GRIDPARAM_H

#include <cassert>
#include <fdm/common.hpp>
#include <fdm/svector.hpp>

namespace fdm
{

// this class contains all the parameters of the non-uniform gird
// with grid spacing based on the generating function
// fdm::FunctionMapping1, i.e.one concentration point,
// see tk diploma thesis page 23
template <std::size_t VecDim>
class GridParam
{
private:
    // 0 .. VecDim-1 are spacial dimensions, VecDim is for time dimension
    SVector<int,VecDim+1>      m_num_points;
    SVector<double,VecDim+1>   m_lower_bound;
    SVector<double,VecDim+1>   m_upper_bound;
    SVector<double,VecDim+1>   m_concentration_point;
    SVector<double,VecDim+1>   m_concentration_ratio;

    // we give the user the chance to not specify upper/lower bounds
    // and concentration info but rather delegate other functions to
    // fill that with appropriate default values
    SVector<bool,VecDim+1>     m_isset_lower;
    SVector<bool,VecDim+1>     m_isset_upper;
    SVector<bool,VecDim+1>     m_isset_concentration;
    SVector<bool,VecDim+1>     m_isset_num_points;

public:
    // removes values except for the number of points data
    // (it actually doesn't remove any data but sets the boolean values
    // to indicate we have not overwritten any defaults)
    void unset()
    {
        for(size_t i=0; i<=VecDim; i++) {
            m_isset_lower[i]=false;
            m_isset_upper[i]=false;
            m_isset_concentration[i]=false;
            m_isset_num_points[i]=false;
        }
    }

    GridParam()
    {
        this->unset();
    }

    void set_lower(int dim, double lower, bool overwrite=true)
    {
        assert(0<=dim && dim<=(int)VecDim);
        if(overwrite==true || m_isset_lower[dim]==false) {
            m_lower_bound[dim]=lower;
            m_isset_lower[dim]=true;
        }
    }
    void set_upper(int dim, double upper, bool overwrite=true)
    {
        assert(0<=dim && dim<=(int)VecDim);
        if(overwrite==true || m_isset_upper[dim]==false) {
            m_upper_bound[dim]=upper;
            m_isset_upper[dim]=true;
        }
    }
    void set_bounds(int dim, double lower, double upper, bool overwrite=true)
    {
        assert(lower<upper);
        set_lower(dim,lower,overwrite);
        set_upper(dim,upper,overwrite);
    }
    void set_concentration(int dim, double point, double ratio,
                           bool overwrite=true)
    {
        assert(0<=dim && dim<=(int)VecDim);
        assert(ratio>0.0);
        if(overwrite==true || m_isset_concentration[dim]==false) {
            m_concentration_point[dim]=point;
            m_concentration_ratio[dim]=ratio;
            m_isset_concentration[dim]=true;
        }
    }
    void set_num_points(int dim, int num, bool overwrite=true)
    {
        assert(0<=dim && dim<=(int)VecDim);
        assert(num>0);
        if(overwrite==true || m_isset_num_points[dim]==false) {
            m_num_points[dim]=num;
            m_isset_num_points[dim]=true;
        }
    }

    void get_bounds(int dim, double& lower, double& upper) const
    {
        assert(0<=dim && dim<=(int)VecDim);
        lower=m_lower_bound[dim];
        upper=m_upper_bound[dim];
    }
    void get_concentration(int dim, double& point, double& ratio) const
    {
        assert(0<=dim && dim<=(int)VecDim);
        point=m_concentration_point[dim];
        ratio=m_concentration_ratio[dim];
    }
    int get_num_points(int dim) const
    {
        assert(0<=dim && dim<=(int)VecDim);
        return m_num_points[dim];
    }
    // returns true if all parameters are set
    bool isset(int dim) const
    {
        assert(0<=dim && dim<=(int)VecDim);
        return m_isset_num_points[dim] && m_isset_lower[dim] &&
               m_isset_upper[dim] && m_isset_concentration[dim];
    }

};




} // namespace fdm


#endif // MF_FDM_GRIDPARAM_H
