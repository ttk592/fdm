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

#ifndef MF_FDM_RESULT_HPP
#define MF_FDM_RESULT_HPP

namespace fin
{

// structure which can save option price and all greeks directly available
// from the grid
struct GridGreeks {
public:
    double	price;		// V(x)
    double	theta;		// d/dt
    double	delta;		// d/dS
    double	gamma;		// d^2/dS^2
    double	vega;		// d/dsigma
    double	volga;		// d^2/dsigma^2
    double	vanna;		// d^2/(dS dsigma)

    GridGreeks(): price(0.0), theta(0.0), delta(0.0), gamma(0.0),
        vega(0.0), volga(0.0), vanna(0.0) {};
};

} // namespace fin

#endif // MF_FDM_RESULT_HPP
