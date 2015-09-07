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



#ifndef MF_FDM_STATS_HPP
#define MF_FDM_STATS_HPP

#include <fdm/common.hpp>                       // for sqr only
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>

using fdm::sqr;

namespace fin
{

// cumulative distriburion function (cdf)
double cdf_norm(double x)
{
    return boost::math::cdf(boost::math::normal(), x);
}
// quantiles, inverse cdf
double q_norm(double x, bool complement=false)
{
    if(complement==false) {
        return boost::math::quantile(boost::math::normal(), x);
    } else {
        return boost::math::quantile(
                   boost::math::complement(boost::math::normal(),x) );
    }
}

// non-central chi squared
double cdf_nchisq(double x, double k, double lambda)
{
    return boost::math::cdf(
               boost::math::non_central_chi_squared(k,lambda), x);
}
double q_nchisq(double x, double k, double lambda, bool complement=false)
{
    if(complement==false) {
        return boost::math::quantile(
                   boost::math::non_central_chi_squared(k,lambda), x);
    } else {
        return boost::math::quantile(boost::math::complement(
                                         boost::math::non_central_chi_squared(k,lambda), x) );
    }
}


// CIR process: dv = kappa*(theta-v)  dt + xi sqrt(v) dW
// we know exp(kappa t)/c(t) v_t is a non-central chisq(k,lambda) with
// k=4kappa*theta/xi^2, lambda=v_0/c(t)
// c(t)=xi^2/(4*kappa) * (exp(kappa*t)-1)
// cdf of a CIR process: P(v_T<x | v_0)
double cdf_cir(double x, double T, double v0,
               double kappa, double theta, double xi )
{

    double k       = 4.0*kappa*theta/(xi*xi);
    double tau     = xi*xi*(exp(kappa*T)-1.0)/(4.0*kappa);
    double lambda  = v0/tau;
    double y       = x*exp(kappa*T)/tau;
    return cdf_nchisq(y,k,lambda);
}
// quantiles of a CIR process
double q_cir(double x, double T, double v0,
             double kappa, double theta, double xi, bool complement=false)
{
    assert( (x>0.0) && (x<1.0) );
    double k       = 4.0*kappa*theta/(xi*xi);
    double tau     = xi*xi*(exp(kappa*T)-1.0)/(4.0*kappa);
    double lambda  = v0/tau;
    return q_nchisq(x,k,lambda,complement)*tau*exp(-kappa*T);
}

// quantiles of an exponential brownian motion
double q_bs(double x, double T, double s0,
            double mu, double sigma, bool complement=false)
{
    return s0*exp((mu-0.5*sigma*sigma)*T + sqrt(T)*sigma*q_norm(x,complement));
}

double var_barrier(double eps, double v_0, double T, double xi)
{
    double gamma=0.5*xi*xi;
    return std::max(0.01, sqr( sqrt(v_0) +
                               sqrt(-gamma*T*log(0.1*eps*sqrt(4.0*M_PI*gamma*T))) ) );
}
double logspot_range(double eps, double v, double T)
{
    double gamma=0.5*v;
    return sqrt(-4.0*gamma*T*log(eps*sqrt(4.0*M_PI*gamma*T)));
}




// given pairs of index --> x (e.g. fdm grid parameter --> fdm result)
// assume the resulting value x is normally distributed as follows, and
// x_i are independent
//
//    (x-exact)/excact = N(mu, sigma^2) / index^order
//
void estim_error(std::vector<double> index, std::vector<double> x, double order,
                 bool exact_is_known, bool mu_is_known,
                 double& exact, double& mu, double& sigma)
{
    assert(index.size()==x.size() && index.size()>1);
    int n=index.size();
    if(exact_is_known) {
        // if "exact" is know, everything is simple as it can be transformed
        // to iid N(mu, sigma)
        double sum=0.0;
        std::vector<double>  z(n);
        for(int i=0; i<n; i++) {
            z[i]=(x[i]-exact)/exact * pow(index[i],order);
            sum+=z[i];
        }
        if(mu_is_known==false) {
            mu=sum/n;
        }
        sum=0.0;
        for(int i=0; i<n; i++) {
            sum+=(z[i]-mu)*(z[i]-mu);
        }
        if(mu_is_known) {
            sigma=sqrt(sum/n);
        } else {
            sigma=sqrt(sum/(n-1));
        }

    } else {
        // exact is not known
        // apply maximum likelihood method on
        // x = N( exact + mu/index^order, (sigma/index^order)^2 )
        double sum_i2k=0.0;
        double sum_ik=0.0;
        double sum_xi2k=0.0;
        double sum_xik=0.0;
        for(int i=0; i<n; i++) {
            sum_i2k  += pow(index[i],2.0*order);
            sum_ik   += pow(index[i],order);
            sum_xi2k += x[i]*pow(index[i],2.0*order);
            sum_xik  += x[i]*pow(index[i],order);
        }
        if(mu_is_known) {
            exact = (sum_xi2k-mu*sum_ik) /sum_i2k;
        } else {
            exact = (sum_xi2k-sum_xik*sum_ik/n) / (sum_i2k-sum_ik*sum_ik/n);
            mu = (sum_xik-sum_ik*exact) / n;
        }

        double sum=0.0;
        for(int i=0; i<n; i++) {
            sum+=sqr(x[i]-exact-mu/pow(index[i],order)) * pow(index[i],2.0*order);
        }
        sigma=sqrt(sum/(n-1));     // mlh only gives sqrt(sum/n)
        // TODO, check (n-1) is the right thing to do

        // since we assumed a slightly different model for easy mlh
        // we need to divide by the estimated exact value to get the original
        // mu and sigma
        mu/=exact;
        sigma/=exact;
    }
}

} // namespace fin



#endif /* MF_FDM_STATS_HPP */
