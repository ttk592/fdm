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


#ifndef MF_FDM_GRIDFACTORY_HPP
#define MF_FDM_GRIDFACTORY_HPP

#include <cassert>
#include <fdm/common.hpp>
#include <fdm/svector.hpp>
#include <fdm/functors.hpp>
#include <fdm/grid.hpp>
#include <fdm/gridparam.hpp>

namespace fdm
{


// create a space grid based on given parameters
template<std::size_t VecDim>
void set_grid(GridC2<VecDim>& grid, const GridParam<VecDim>& param)
{
    fdm::FunctionMapping1   mapping;    // generating function for grid points

    for(size_t dim=0; dim<VecDim; dim++) {
        assert(param.isset(dim)==true);
        double lower, upper, point, ratio;
        param.get_bounds(dim, lower, upper);
        param.get_concentration(dim, point, ratio);
        mapping.set_param(lower, upper, point, ratio);
        grid.set_by_mapping(dim,mapping,param.get_num_points(dim));
    }
    grid.calculate();
}

// create a time grid based on given parameters
// TODO: possibly change to std::vector<double>& t
template<std::size_t VecDim>
void set_time_grid(Grid<1>& grid, const GridParam<VecDim>& param)
{
    assert(param.isset(VecDim)==true);
    fdm::FunctionMapping1   mapping;    // generating function for grid points

    double lower, upper, point, ratio;
    param.get_bounds(VecDim, lower, upper);
    param.get_concentration(VecDim, point, ratio);
    mapping.set_param(lower, upper, point, ratio);
    grid.set_by_mapping(0,mapping,param.get_num_points(VecDim));
}



} // namespace fdm


#endif /* MF_FDM_GRIDFACTORY_HPP */
