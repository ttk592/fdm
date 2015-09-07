#ifndef HESTON_ANALYTIC_HPP
#define HESTON_ANALYTIC_HPP

namespace fin
{

double heston_analytic(double rd, double rf, double spot, double strike,
                       double T, double kappa, double theta, double xi,
                       double v0, double rho, int putcall);

double heston_lipton(double rd, double rf, double spot, double strike,
                     double T, double kappa, double theta, double xi,
                     double v0, double rho, int putcall,
                     double barrier_lower=-1.0, double barrier_upper=-1.0);

}

#endif /* HESTON_ANALYTIC_HPP */
