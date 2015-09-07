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

#include <cstdlib>
#include <cstdio>
#include <boost/timer.hpp>

#include <fdm/gridparam.hpp>
#include "market.hpp"
#include "option.hpp"
#include "result.hpp"
#include "heston.hpp"
#include "heston_analytic.hpp"

using fdm::sqr;


// ======================================================================
//         Name:  main
//  Description:  main function, executed at runtime
// ======================================================================
int main(int argc, char** argv)
{

    if(argc<1) {
        printf("usage: %s ...\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    fdm::enableexcept();    // for testing, raise fpu exceptions

    // market parameters
    fin::MarketHeston   market;
    market.spot=1.0;                // current spot at time t
    market.vol=0.12;                // initial vol
    market.theta=sqr(market.vol);   // long term variance
    market.kappa=3.;                // speed of mean-reversion
    market.xi=0.25;                 // vol of the variance process
    market.rho=0.2;                 // correllation
    market.rd=0.03;                 // domestic rate
    market.rf=0.01;                 // foreign rate / continuous yield

    {
        // Asian option
        // ------------
        const int spotdim=0;    // fin::PDEHeston::dim_spot;
        const int voldim=1;     // fin::PDEHeston::dim_vol;
        const int avgdim=2;
        const int timedim=3;

        fin::OptionAsian    option;
        option.putcall=fin::type::call;
        option.fixed_float=fin::type::fixed_strike;
        option.T=1.0;       // time to maturity
        option.t=0.2;       // current time
        // note market.spot is the spot at time t, not 0
        option.strike=1.0;  // strike or strike factor
        // averaging schedule
        option.dates = {0.75, 0.80, 0.85, 0.90, 0.95, 1.00};
        // for past dates, we need to know spot fixes
        option.spot =  {1.00, 0.00, 0.00, 0.00, 0.00, 0.00};

        // numerical parameters
        fdm::GridParam<3>   gridparam;
        gridparam.set_num_points(spotdim,40);
        gridparam.set_num_points(voldim,15);
        gridparam.set_num_points(avgdim,15);
        gridparam.set_num_points(timedim,20);
        // can set more grid parameters here, but if not they will be
        // filled with defaults

        // FDM calculation
        boost::timer	timer;
        timer.restart();
        fin::GridGreeks     result;
        result=fin::heston_fdm(market, option, gridparam);
        double elapsed=timer.elapsed();
        printf("asian: %ix%ix%ix%i, fdm: %f (%.3fs)\n",
               gridparam.get_num_points(0), gridparam.get_num_points(1),
               gridparam.get_num_points(2), gridparam.get_num_points(3),
               result.price,elapsed);

    }



    {
        // Barrier option
        // --------------
        const int spotdim=0; //fin::PDEHeston::dim_spot;
        const int voldim=1; // fin::PDEHeston::dim_vol;
        const int timedim=2;

        fin::OptionBarrier  option;
        option.putcall=fin::type::call;
        option.T=1.0;               // time to maturity
        option.strike=1.0;	    // strike
        option.barrier_lower=0.0;   // 0.0 if no barrier
        option.barrier_upper=0.0;   // 0.0 if no barrier
        option.rebate_lower=0.0;    // rebate if hitting lower barrier
        option.rebate_upper=0.0;    // rebate if hitting upper barrier

        // numerical parameters
        fdm::GridParam<2>   gridparam;
        gridparam.set_num_points(spotdim,100);
        gridparam.set_num_points(voldim,50);
        gridparam.set_num_points(timedim,30);

        // FDM calculation
        boost::timer	timer;
        timer.restart();
        fin::GridGreeks     result;
        result=fin::heston_fdm(market, option, gridparam);
        double elapsed=timer.elapsed();

        // compare against semi-analytic results
        double exact=0.0;
        if(option.barrier_lower<=0.0 && option.barrier_upper<=0.0) {
            exact=fin::heston_analytic(market.rd,market.rf,market.spot,
                                       option.strike,option.T,market.kappa,
                                       market.theta, market.xi,
                                       sqr(market.vol),market.rho,
                                       option.putcall);
        } else if(market.rd==market.rf && market.rho==0.0 &&
                  option.rebate_lower==0.0 && option.rebate_upper==0.0) {
            exact=fin::heston_lipton(market.rd,market.rf,market.spot,
                                     option.strike,option.T,market.kappa,
                                     market.theta,
                                     market.xi,sqr(market.vol),market.rho,
                                     option.putcall,option.barrier_lower,
                                     option.barrier_upper);
        }


        printf("barrier: %ix%ix%i, fdm: %f, exact: %f, rel err: %.3f%% (%.3fs)\n",
               gridparam.get_num_points(0), gridparam.get_num_points(1),
               gridparam.get_num_points(2),
               result.price,exact, (result.price-exact)/exact*100.0 ,elapsed);
    }




    return EXIT_SUCCESS;
}


