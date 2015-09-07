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


#include "heston.hpp"

#include <cstdio>
#include <cstdlib>

#include <fdm/svector.hpp>
#include <fdm/interpolation.hpp>
#include <fdm/grid.hpp>
#include <fdm/gridparam.hpp>
#include <fdm/gridfactory.hpp>
#include <fdm/solve_pde.hpp>

#include "heston_pde.hpp"
#include "market.hpp"
#include "option.hpp"
#include "result.hpp"
#include "stats.hpp"

namespace fin
{


// barrier option pricing
// ----------------------
// calculates the price of an option and returns price sensitivities
// which are easily obtainable by looking at the grid
// (all except rho_d, rho_f)
GridGreeks heston_fdm( const fin::MarketHeston& market,
                       const fin::OptionBarrier& option,
                       fdm::GridParam<2>& param)
{

#ifndef NDEBUG
    // check Feller condition
    if(market.kappa*market.theta<0.5*market.xi*market.xi) {
        printf("variance process can reach 0! (kappa=%f, theta=%f, xi=%f)\n",
               market.kappa,market.theta,market.xi);
    }
#endif

    const int spotdim=fin::PDEHeston::dim_spot;
    const int voldim=fin::PDEHeston::dim_vol;
    const int timedim=2;

    // central point where we'll read the value from the grid
    fdm::SVector<double,2> x;
    x[spotdim]=log(market.spot);
    x[voldim]=market.vol*market.vol;

    // set default grid parameter in case gridparam didn't have them
    // set already
    {
        // find upper and lower bounds based on properties of the corresponding
        // stochastic process
        double quantile=1e-12;
        double vol_avg=std::max(market.vol,
                                sqrt(fin::q_cir(0.9,option.T, market.vol*market.vol,
                                                market.kappa, market.theta, market.xi)) );

        double var_max=fin::q_cir(quantile,option.T,
                                  market.vol*market.vol, market.kappa,
                                  market.theta, market.xi,true);
        double xrange=-fin::q_norm(quantile)*sqrt(option.T)*vol_avg;

        // default spot parameters (do not overwrite, if values already set)
        param.set_bounds(spotdim, x[spotdim]-xrange, x[spotdim]+xrange,
                         false);
        param.set_concentration(spotdim, x[spotdim], 0.1, false);

        // if option has barriers we have to enforce grid boundaries
        if(option.barrier_lower>0.0) {
            param.set_lower(spotdim, log(option.barrier_lower), true);
        }
        if(option.barrier_upper>0.0) {
            param.set_upper(spotdim, log(option.barrier_upper), true);
        }

        // default vol parameters (do not overwrite)
        param.set_bounds(voldim, 0.0, var_max, false);
        param.set_concentration(voldim, x[voldim], 0.2, false);

        // delfault time bounds (overwrite bounds)
        param.set_bounds(timedim, 0.0, option.T, true);
        param.set_concentration(timedim, 0.0, 0.2, false);
    }


    // set the space and time grid based on parameters of gridparam
    fdm::Grid<1>	grid_time;      // actual time grid
    fdm::GridC2<2>	grid_space;     // actual space grid

    fdm::set_grid(grid_space, param);
    fdm::set_time_grid(grid_time, param);
    // we may not want dt to change every step, so we make
    // the grid piecewise uniform, the second parameter says
    // how many different values of dt we want to have max
    //grid_time.changeGridToPiecewiseUniform(0,10);


    // setting up the pde coefficients
    fin::PDEHeston  pde;
    pde.import_params(market);

    // defining the boundary conditions
    fdm::BoundaryConst<2>	boundary;
    boundary.TypeLower[spotdim]=fdm::BoundaryType::Free;
    boundary.TypeUpper[spotdim]=fdm::BoundaryType::Free;
    if(option.barrier_lower>0.0) {
        boundary.TypeLower[spotdim]=fdm::BoundaryType::Value;
        boundary.ValueLower[spotdim]=option.rebate_lower;
    }
    if(option.barrier_upper>0.0) {
        boundary.TypeUpper[spotdim]=fdm::BoundaryType::Value;
        boundary.ValueUpper[spotdim]=option.rebate_upper;
    }
    boundary.TypeLower[voldim]=fdm::BoundaryType::Free;
    boundary.TypeUpper[voldim]=fdm::BoundaryType::Free;


    // setting up the grid function object with initial conditions
    fdm::GridFunction<2>    u(grid_space);

    u.iterator.begin();
    while(u.iterator.isend()==false) {
        double s=exp(u.coord(spotdim));
        u.value()=option.payoff(s);
        u.iterator++;
    }


    // solving the pde
    std::vector<double>	price(grid_time.size(0)); // price at all times t_i
    price[0]=u(x);
    for(size_t i=1; i<grid_time.size(0); i++) {
        double t1=grid_time.coord(0,i-1);
        double t2=grid_time.coord(0,i);

        // solve pde for one time-step
        double theta = (i<=1) ? 1.0 : 0.5;  // first step implicit, then CN
        bool recalc=true;           // TODO: set to false if dt remains const
        fdm::FdmSolveStep(u,pde,boundary,t1,t2,theta,recalc,
                          fdm::type::pred_corr);
        price[i]=u(x);
    }

    // obtaining price sensitivites
    GridGreeks	value;
    int		n=grid_time.size(0)-1;
    value.price = price.at(n);
    value.theta = (price.at(n) - price.at(n-1)) /
                  (grid_time.coord(0,n) - grid_time.coord(0,n-1));
    // we use normalised coordinates in s, ie x[0]=log(S),
    // so we have to take that into account for calculating derivatives
    value.delta = u.f_1(spotdim,x)/market.spot;
    value.gamma = (u.f_2(spotdim,spotdim,x)-u.f_1(spotdim,x)) /
                  (market.spot*market.spot);
    value.vega  = u.f_1(voldim,x);
    value.volga = u.f_2(voldim,voldim,x);
    value.vanna = u.f_2(spotdim,voldim,x)/market.spot;

    return value;
}



// Asian option pricing
// --------------------
// calculates the price of an option and returns price sensitivities
// which are easily obtainable by looking at the grid

// applies continuity of option price, when a spot for an averaging
// date is observed
// average is updated as follows: a --> (a*n+s)/(n+1)
// with log-space
// U(s,v,a,t-) = U(s,v,(a*n+s)/(n+1),t+)    (t non-reversed time)
void asian_continuity(fdm::GridFunction<3>& u, int n_pastobs)
{
    assert(n_pastobs>=0);
    const int   spotdim=fin::PDEHeston::dim_spot;
    const int   voldim=fin::PDEHeston::dim_vol;
    const int   avgdim=2;
    fdm::SVector<int, 3> idx;
    fdm::Cspline spline;
    for(size_t i=0; i<u.size(spotdim); i++) {
        idx[spotdim]=i;
        double s = exp(u.coord(spotdim,i));
        for(size_t j=0; j<u.size(voldim); j++) {
            idx[voldim]=j;
            // define 1d-spline along avg-axis, for fixed spot,vol
            size_t n=u.size(avgdim);
            spline.resize(n);
            for(size_t k=0; k<n; k++) {
                idx[avgdim]=k;
                // set x-y grid values for spline
                double a = exp(u.coord(avgdim,k));
                spline.set_point(k,a,u.value(idx));
            }
            spline.calculate();
            // interpolate: u(s,v,a,t+)=spline((a*n+s)/(n+1))
            for(size_t k=0; k<u.size(avgdim); k++) {
                idx[avgdim]=k;
                double a = exp(u.coord(avgdim,k));
                //printf("s=%f, v=%f, a=%f, u(x)=%f --> ",
                //        s,u.coord(voldim,idx[voldim]),a,u.value(idx) );
                u.value(idx) = spline((a*n_pastobs+s)/(double)(n_pastobs+1));
                //printf("%f (a=%f)\n", u.value(idx),
                //       (a*n_pastobs+s)/(double)(n_pastobs+1) );
            }
        }
    }
}


GridGreeks heston_fdm( const fin::MarketHeston& market,
                       const fin::OptionAsian& option,
                       fdm::GridParam<3>& param)
{

#ifndef NDEBUG
    // check Feller condition
    if(market.kappa*market.theta<0.5*market.xi*market.xi) {
        printf("variance process can reach 0! (kappa=%f, theta=%f, xi=%f)\n",
               market.kappa,market.theta,market.xi);
    }
#endif

    const int   spotdim=fin::PDEHeston::dim_spot;
    const int   voldim=fin::PDEHeston::dim_vol;
    const int   avgdim=2;
    const int   timedim=3;
    int         n_pastobs, n_obs;
    double      avg;

    // need to distinguish three cases
    // (1) t < obs_dates[0], i.e. no average observation has yet occured
    // (2) obs_dates[0] <= t < obs_dates[n-1], i.e. avg obs has occured
    //     and more are to come
    // (3) obs_dates[n-1] <= t, i.e. no more avg obs to follow
    //     we ignore this case here, because it reduces to a vanilla
    //     option, and we are lazy

    assert(option.t<option.dates.back());   // do not allow case (3)

    option.current_avg(avg, n_pastobs, n_obs);
    if(n_pastobs==0)
        avg=market.spot;

    // central point where we'll read the value from the grid
    fdm::SVector<double,3> x;
    x[spotdim]=log(market.spot);
    x[voldim]=sqr(market.vol);
    x[avgdim]=log(avg);

    // set default grid parameter in case gridparam didn't have them
    // set already
    {
        // find upper and lower bounds based on properties of the corresponding
        // stochastic process
        double quantile=1e-12;
        double vol_avg=std::max(market.vol,
                                sqrt(fin::q_cir(0.9,option.T, sqr(market.vol),
                                                market.kappa, market.theta, market.xi)) );

        double var_max=fin::q_cir(quantile,option.T,
                                  market.vol*market.vol, market.kappa,
                                  market.theta, market.xi,true);
        double xrange=-fin::q_norm(quantile)*sqrt(option.T)*vol_avg;

        // for the auxilliary variable "average" the bounds depend on
        // the future observation schedule
        /*
        double sum_low=0.0;
        double sum_up =0.0;
        for(size_t i=0; i<option.dates.size(); i++){
            if(option.t<option.dates[i]){
                double mu=market.rd-market.rf-0.5*sqr(vol_avg);
                double T=option.dates[i]-option.t;
                double alpha=-fin::q_norm(quantile);
                sum_low += market.spot*exp(mu*T - alpha*vol_avg*sqrt(T));
                sum_up  += market.spot*exp(mu*T + alpha*vol_avg*sqrt(T));
            }
        }
        double logavg_min=log( (avg*n_pastobs+sum_low)/n_obs );
        double logavg_max=log( (avg*n_pastobs+sum_up) /n_obs );
        */
        double logavg_min=x[spotdim]-xrange;        // same as spot
        double logavg_max=x[spotdim]+xrange;        // same as spot

        // default spot parameters (do not overwrite, if values already set)
        param.set_bounds(spotdim, x[spotdim]-xrange, x[spotdim]+xrange, false);
        param.set_concentration(spotdim, x[spotdim], 0.1, false);

        // default vol parameters (do not overwrite)
        param.set_bounds(voldim, 0.0, var_max, false);
        param.set_concentration(voldim, x[voldim], 0.2, false);

        // default average auxilliary variable parameters
        param.set_bounds(avgdim, logavg_min, logavg_max, false);
        param.set_concentration(avgdim, x[avgdim], 0.5, false);

        // delfault time bounds (overwrite bounds)
        param.set_bounds(timedim, 0.0, option.T-option.t, true);
        param.set_concentration(timedim, 0.0, 0.9, false);
    }


    // set the space and time grid based on parameters of gridparam
    fdm::Grid<1>	grid_time;      // actual time grid
    fdm::GridC2<3>	grid_space;     // actual space grid

    fdm::set_grid(grid_space, param);
    fdm::set_time_grid(grid_time, param);

    // add all observation dates to the time grid (reverse time)
    std::vector<double> avg_dates;        // avg dates, reverse time
    for(int i=(int)option.dates.size()-1; i>=0; i--) {
        double ti=option.T-option.dates[i];
        if(0<ti && ti<option.T-option.t) {
            grid_time.insert(0,ti);
        }
        avg_dates.push_back(ti);
    }

    // setting up the pde coefficients
    fin::PDEHeston  pde;
    pde.import_params(market);

    // defining the boundary conditions
    fdm::BoundaryConst<2>	boundary;
    boundary.TypeLower[spotdim]=fdm::BoundaryType::Free;
    boundary.TypeUpper[spotdim]=fdm::BoundaryType::Free;
    boundary.TypeLower[voldim]=fdm::BoundaryType::Free;
    boundary.TypeUpper[voldim]=fdm::BoundaryType::Free;

    // setting up the grid function object with initial conditions,
    // i.e. terminal payoff
    fdm::GridFunction<3>    u(grid_space);
    u.iterator.begin();
    while(u.iterator.isend()==false) {
        double s=exp(u.coord(spotdim));
        double a=exp(u.coord(avgdim));
        u.value()=option.payoff(s,a);
        u.iterator++;
    }

    // applying continuity condition if maturity is an avg obs date
    double n_pastobs_grid=n_obs-1;
    if(option.dates.back()  == option.T) {
        asian_continuity(u,n_pastobs_grid);
        n_pastobs_grid--;
    }

    // solving the pde
    std::vector<double>	price(grid_time.size(0)); // price at all times t_i
    price[0]=u(x);
    size_t i_cont=grid_time.size(0);  // in case we need to continue on 2d grid
    for(size_t i=1; i<grid_time.size(0); i++) {
        double t1=grid_time.coord(0,i-1);
        double t2=grid_time.coord(0,i);

        // solve pde for one time-step
        double theta = (i<=1) ? 1.0 : 0.5;  // first step implicit, then CN
        bool recalc=true;
        fdm::FdmSolveStep(u,pde,boundary,t1,t2,theta,recalc,
                          fdm::type::pred_corr);

        // break if we are at the first avg obs date (last in reverse time)
        if(t2==avg_dates.back()) {
            i_cont=i+1;
            break;
        }

        // if t2 is an observation date for averaging
        if( std::binary_search(avg_dates.begin(), avg_dates.end(), t2) ) {
            asian_continuity(u,n_pastobs_grid);
            n_pastobs_grid--;
        }
        price[i]=u(x);
    }

    GridGreeks	value;

    // in case of (1), i.e. no avg obs have occured yet, we continue
    // with a 2-d grid without auxilliary variables
    if(option.t<option.dates[0]) {
        // create a 2d-grid without auxilliary variables
        fdm::GridC2<2>	grid_space2d(grid_space,3); // spot-vol grid
        fdm::GridFunction<2> v(grid_space2d);       // 2d function v(spot,vol)

        // initialise based on V(spot,vol)=U(spot,vol,spot);
        v.iterator.begin();
        while(v.iterator.isend()==false) {
            double spot=v.coord(spotdim);      // actually logspot
            double vol=v.coord(voldim);        // actually var
            // note u() currently does not allow for extrapolation
            // --> require avg-grid is not smaller than spot-grid
            // (alterantively rewrite using 1d-splines which allow extrapol)
            v.value() = u( fdm::SVector<double,3>({spot, vol, spot}) );
            v.iterator++;
        }

        // solve 2-d pde
        fdm::SVector<double, 2> x2d = {x[0], x[1]};
        for(size_t i=i_cont; i<grid_time.size(0); i++) {
            double t1=grid_time.coord(0,i-1);
            double t2=grid_time.coord(0,i);

            // solve pde for one time-step
            double theta = 0.5;
            bool recalc=true;
            fdm::FdmSolveStep(v,pde,boundary,t1,t2,theta,recalc,
                              fdm::type::pred_corr);
            price[i]=v(x2d);
        }
        // obtaining price sensitivites
        int		n=grid_time.size(0)-1;
        value.price = price.at(n);
        value.theta = (price.at(n) - price.at(n-1)) /
                      (grid_time.coord(0,n) - grid_time.coord(0,n-1));
        // we use normalised coordinates in s, ie x[0]=log(S),
        // so we have to take that into account for calculating derivatives
        value.delta = v.f_1(spotdim,x2d)/market.spot;
        value.gamma = (v.f_2(spotdim,spotdim,x2d)-v.f_1(spotdim,x2d)) /
                      (market.spot*market.spot);
        value.vega  = v.f_1(voldim,x2d);
        value.volga = v.f_2(voldim,voldim,x2d);
        value.vanna = v.f_2(spotdim,voldim,x2d)/market.spot;
    } else {
        // obtaining price sensitivites
        int		n=grid_time.size(0)-1;
        value.price = price.at(n);
        value.theta = (price.at(n) - price.at(n-1)) /
                      (grid_time.coord(0,n) - grid_time.coord(0,n-1));
        // we use normalised coordinates in s, ie x[0]=log(S),
        // so we have to take that into account for calculating derivatives
        value.delta = u.f_1(spotdim,x)/market.spot;
        value.gamma = (u.f_2(spotdim,spotdim,x)-u.f_1(spotdim,x)) /
                      (market.spot*market.spot);
        value.vega  = u.f_1(voldim,x);
        value.volga = u.f_2(voldim,voldim,x);
        value.vanna = u.f_2(spotdim,voldim,x)/market.spot;
    }

    return value;
}








void read_param(fin::MarketHeston& market, fin::OptionBarrier& option,
                fdm::GridParam<2>& param,
                int argc, char** argv)
{
    double num_spot, num_vol, num_time;
    // reading command line input or set default values if none supplied
    if(argc<18) {
        // use some default values
        market.spot=1.38;
        market.vol=sqrt(0.0049);
        market.theta=0.0049;
        market.kappa=2.;
        market.xi=0.33147;
        market.rho=0.2;
        market.rd=0.;
        market.rf=0.;

        option.putcall=type::call;
        option.T=0.1;               // time to maturity
        option.strike=1.38;	    // strike
        option.barrier_lower=0.0;   // 0.0 if no barrier
        option.barrier_upper=0.0;   // 0.0 if no barrier
        option.rebate_lower=0.0;    // rebate if hitting lower barrier
        option.rebate_upper=0.0;    // rebate if hitting upper barrier

        num_spot=100;
        num_vol=50;
        num_time=30;

    } else {
        // use command line arguments
        market.spot=atof(argv[1]);
        market.vol=sqrt(atof(argv[2]));
        market.theta=atof(argv[3]);
        market.kappa=atof(argv[4]);
        market.xi=atof(argv[5]);
        market.rho=atof(argv[6]);
        market.rd=atof(argv[7]);
        market.rf=atof(argv[8]);

        option.putcall=static_cast<type::put_call>(atoi(argv[9]));  // put/call
        option.T=atof(argv[10]);		// time to maturity
        option.strike=atof(argv[11]);           // strike
        option.barrier_lower=atof(argv[12]);    // 0.0 if no barrier
        option.barrier_upper=atof(argv[13]);    // 0.0 if no barrier
        option.rebate_lower=atof(argv[14]);     // rebate
        option.rebate_upper=atof(argv[15]);     // rebate

        num_spot=atof(argv[16]);
        num_vol=atof(argv[17]);
        num_time=atof(argv[18]);
    }

    const int spotdim=fin::PDEHeston::dim_spot;
    const int voldim=fin::PDEHeston::dim_vol;
    const int timedim=2;

    param.set_num_points(timedim,num_time);
    param.set_num_points(spotdim,num_spot);
    param.set_num_points(voldim,num_vol);
}


} // namespace fin
