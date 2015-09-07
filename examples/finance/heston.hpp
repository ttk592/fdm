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


#ifndef MF_FDM_HESTON_HPP
#define MF_FDM_HESTON_HPP

#include <fdm/gridparam.hpp>
#include "market.hpp"
#include "option.hpp"
#include "result.hpp"

namespace fin
{

// barrier option pricing
GridGreeks heston_fdm(const fin::MarketHeston& market,
                      const fin::OptionBarrier& option,
                      fdm::GridParam<2>& param);

// Asian option pricing
GridGreeks heston_fdm( const fin::MarketHeston& market,
                       const fin::OptionAsian& option,
                       fdm::GridParam<3>& param);



void read_param(fin::MarketHeston& market, fin::OptionBarrier& option,
                fdm::GridParam<2>& param,
                int argc, char** argv);


} // namespace fin

#endif // MF_FDM_HESTON_HPP

