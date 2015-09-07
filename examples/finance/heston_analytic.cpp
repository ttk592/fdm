// implements the semi-analytic formula for pricing under the Heston model
//
// (1) heston_analytic(), based on QuantLib, for vanilla call/put options,
//     code taken from QuantLib-1.0.1, see license below
//
// (2) heston_lipton(), own code, for barrier options, given rho=0, rd=rf


/*
 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include "heston_analytic.hpp"

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <functional>
#include <limits>


namespace fin
{

// anonymous namespace so we don't export helper functions
namespace
{

// ql/math/distributions/gammadistribution.cpp
double GammaFunction_logValue(double x)
{
    assert(x>0.0);

    const double c1_ = 76.18009172947146;
    const double c2_ = -86.50532032941677;
    const double c3_ = 24.01409824083091;
    const double c4_ = -1.231739572450155;
    const double c5_ = 0.1208650973866179e-2;
    const double c6_ = -0.5395239384953e-5;

    double temp = x + 5.5;
    temp -= (x + 0.5)*std::log(temp);
    double ser=1.000000000190015;
    ser += c1_/(x + 1.0);
    ser += c2_/(x + 2.0);
    ser += c3_/(x + 3.0);
    ser += c4_/(x + 4.0);
    ser += c5_/(x + 5.0);
    ser += c6_/(x + 6.0);
    return -temp+std::log(2.5066282746310005*ser/x);
}



// ql/math/matrixutilities/tqreigendecomposition.hpp
class TqrEigenDecomposition
{
private:
    int iter_;
    std::vector<double> d_;
    std::vector< std::vector<double> > ev_;
    bool offDiagIsZero(int k, const std::vector<double>& e);
public:
    enum EigenVectorCalculation { WithEigenVector,
                                  WithoutEigenVector,
                                  OnlyFirstRowEigenVector
                                };

    enum ShiftStrategy { NoShift,
                         Overrelaxation,
                         CloseEigenValue
                       };

    TqrEigenDecomposition(const std::vector<double>& diag,
                          const std::vector<double>& sub,
                          EigenVectorCalculation calc = WithEigenVector,
                          ShiftStrategy strategy = CloseEigenValue);

    const std::vector<double>& eigenvalues()  const
    {
        return d_;
    }
    const std::vector< std::vector<double> >& eigenvectors() const
    {
        return ev_;
    }

    int iterations() const
    {
        return iter_;
    }

};


// ql/math/matrixutilities/tqreigendecomposition.cpp
TqrEigenDecomposition::TqrEigenDecomposition(
    const std::vector<double>& diag, const std::vector<double>& sub,
    EigenVectorCalculation calc, ShiftStrategy strategy)
    : iter_(0), d_(diag)
{

    int n = diag.size();
    if(calc == WithEigenVector)
        ev_.resize(n);
    else if(calc == WithoutEigenVector)
        ev_.resize(0);
    else
        ev_.resize(1);
    for(size_t i=0; i<ev_.size(); i++) {
        ev_[i].resize(n);
    }


    assert(n == (int) sub.size()+1);

    std::vector<double> e(n, 0.0);
    std::copy(sub.begin(),sub.end(),e.begin()+1);
    for (size_t i=0; i < ev_.size(); ++i) {
        ev_[i][i] = 1.0;
    }

    for (int k=n-1; k >=1; --k) {
        while (!offDiagIsZero(k, e)) {
            int l = k;
            while (--l > 0 && !offDiagIsZero(l,e));
            iter_++;

            double q = d_[l];
            if (strategy != NoShift) {
                // calculated eigenvalue of 2x2 sub matrix of
                // [ d_[k-1] e_[k] ]
                // [  e_[k]  d_[k] ]
                // which is closer to d_[k+1].
                // FLOATING_POINT_EXCEPTION
                const double t1 = std::sqrt(
                                      0.25*(d_[k]*d_[k] + d_[k-1]*d_[k-1])
                                      - 0.5*d_[k-1]*d_[k] + e[k]*e[k]);
                const double t2 = 0.5*(d_[k]+d_[k-1]);

                const double lambda =
                    (std::fabs(t2+t1 - d_[k]) < std::fabs(t2-t1 - d_[k]))?
                    t2+t1 : t2-t1;

                if (strategy == CloseEigenValue) {
                    q-=lambda;
                } else {
                    q-=((k==n-1)? 1.25 : 1.0)*lambda;
                }
            }

            // the QR transformation
            double sine = 1.0;
            double cosine = 1.0;
            double u = 0.0;

            bool recoverUnderflow = false;
            for (int i=l+1; i <= k && !recoverUnderflow; ++i) {
                const double h = cosine*e[i];
                const double p = sine*e[i];

                e[i-1] = std::sqrt(p*p+q*q);
                if (e[i-1] != 0.0) {
                    sine = p/e[i-1];
                    cosine = q/e[i-1];

                    const double g = d_[i-1]-u;
                    const double t = (d_[i]-g)*sine+2*cosine*h;

                    u = sine*t;
                    d_[i-1] = g + u;
                    q = cosine*t - h;

                    for (size_t j=0; j < ev_.size(); ++j) {
                        const double tmp = ev_[j][i-1];
                        ev_[j][i-1] = sine*ev_[j][i] + cosine*tmp;
                        ev_[j][i] = cosine*ev_[j][i] - sine*tmp;
                    }
                } else {
                    // recover from underflow
                    d_[i-1] -= u;
                    e[l] = 0.0;
                    recoverUnderflow = true;
                }
            }

            if (!recoverUnderflow) {
                d_[k] -= u;
                e[k] = q;
                e[l] = 0.0;
            }
        }
    }

    // sort (eigenvalues, eigenvectors),
    // code taken from symmetricSchureDecomposition.cpp
    std::vector<std::pair<double, std::vector<double> > > temp(n);
    std::vector<double> eigenVector(ev_.size());
    for (int i=0; i<n; i++) {
        for(size_t j=0; j<ev_.size(); j++) {
            eigenVector[j]=ev_[j][i];
        }
        temp[i] = std::make_pair(d_[i], eigenVector);
    }
    std::sort(temp.begin(), temp.end(),
              std::greater<std::pair<double, std::vector<double> > >());
    // first element is positive
    for (int i=0; i<n; i++) {
        d_[i] = temp[i].first;
        double sign = 1.0;
        if (ev_.size() > 0 && temp[i].second[0]<0.0)
            sign = -1.0;
        for (size_t j=0; j<ev_.size(); ++j) {
            ev_[j][i] = sign * temp[i].second[j];
        }
    }
}

// ql/math/matrixutilities/tqreigendecomposition.cpp
// see NR for abort assumption as it is
// not part of the original Wilkinson algorithm
bool TqrEigenDecomposition::offDiagIsZero(int k,
        const std::vector<double>& e)
{
    return std::fabs(d_[k-1])+std::fabs(d_[k])
           == std::fabs(d_[k-1])+std::fabs(d_[k])+std::fabs(e[k]);
}



// ql/math/integrals/gaussianquadratures.hpp
//
//! orthogonal polynomial for Gaussian quadratures
/*! References:
    Gauss quadratures and orthogonal polynomials

    G.H. Gloub and J.H. Welsch: Calculation of Gauss quadrature rule.
    Math. Comput. 23 (1986), 221-230

    "Numerical Recipes in C", 2nd edition,
    Press, Teukolsky, Vetterling, Flannery,

    The polynomials are defined by the three-term recurrence relation
    \f[
    P_{k+1}(x)=(x-\alpha_k) P_k(x) - \beta_k P_{k-1}(x)
    \f]
    and
    \f[
    \mu_0 = \int{w(x)dx}
    \f]
*/
class GaussianOrthogonalPolynomial
{
public:
    virtual ~GaussianOrthogonalPolynomial() {}
    virtual double mu_0()        const = 0;
    virtual double alpha(int i) const = 0;
    virtual double beta(int i)  const = 0;
    virtual double w(double x)     const = 0;

    double value(int i, double x) const;
    double weightedValue(int i, double x) const;
};


// ql/math/integrals/gaussianorthogonalpolynomial.hpp
//! Gauss-Laguerre polynomial
class GaussLaguerrePolynomial : public GaussianOrthogonalPolynomial
{
public:
    GaussLaguerrePolynomial(double s = 0.0);

    double mu_0() const;
    double alpha(int i) const;
    double beta(int i) const;
    double w(double x) const;

private:
    const double s_;
};

// ql/math/integrals/gaussianorthogonalpolynomial.cpp
double GaussianOrthogonalPolynomial::value(int n, double x) const
{
    if (n > 1) {
        return  (x-alpha(n-1)) * value(n-1, x)
                - beta(n-1) * value(n-2, x);
    } else if (n == 1) {
        return x-alpha(0);
    }

    return 1;
}

// ql/math/integrals/gaussianorthogonalpolynomial.cpp
double GaussianOrthogonalPolynomial::weightedValue(int n, double x) const
{
    return std::sqrt(w(x))*value(n, x);
}

// ql/math/integrals/gaussianorthogonalpolynomial.cpp
GaussLaguerrePolynomial::GaussLaguerrePolynomial(double s)
    : s_(s)
{
    assert(s > -1.0);
}

// ql/math/integrals/gaussianorthogonalpolynomial.cpp
double GaussLaguerrePolynomial::mu_0() const
{
    return std::exp(GammaFunction_logValue(s_+1.0));
}

// ql/math/integrals/gaussianorthogonalpolynomial.cpp
double GaussLaguerrePolynomial::alpha(int i) const
{
    return 2*i+1+s_;
}

// ql/math/integrals/gaussianorthogonalpolynomial.cpp
double GaussLaguerrePolynomial::beta(int i) const
{
    return i*(i+s_);
}

// ql/math/integrals/gaussianorthogonalpolynomial.cpp
double GaussLaguerrePolynomial::w(double x) const
{
    return std::pow(x, s_)*std::exp(-x);
}


// ql/math/integrals/gaussianquadratures.hpp
//
//! Integral of a 1-dimensional function using the Gauss quadratures method
/*! References:
    Gauss quadratures and orthogonal polynomials

    G.H. Gloub and J.H. Welsch: Calculation of Gauss quadrature rule.
    Math. Comput. 23 (1986), 221-230

    "Numerical Recipes in C", 2nd edition,
    Press, Teukolsky, Vetterling, Flannery,

    \test the correctness of the result is tested by checking it
          against known good values.
*/
class GaussianQuadrature
{
public:
    GaussianQuadrature(int n,
                       const GaussianOrthogonalPolynomial& p);

    template <class F>
    double operator()(const F& f) const
    {
        double sum = 0.0;
        for (int i = order()-1; i >= 0; --i) {
            sum += w_[i] * f(x_[i]);
        }
        return sum;
    }

    size_t order() const
    {
        return x_.size();
    }

private:
    std::vector<double> x_, w_;
};

// ql/math/integrals/gaussianquadratures.hpp
//
//! generalized Gauss-Laguerre integration
/*! This class performs a 1-dimensional Gauss-Laguerre integration.
    \f[
    \int_{0}^{\inf} f(x) \mathrm{d}x
    \f]
    The weighting function is
    \f[
        w(x;s)=x^s \exp{-x}
    \f]
    and \f[ s > -1 \f]
*/
class GaussLaguerreIntegration : public GaussianQuadrature
{
public:
    GaussLaguerreIntegration(int n, double s = 0.0)
        : GaussianQuadrature(n, GaussLaguerrePolynomial(s)) {}
};


// ql/math/integrals/gaussianquadratures.cpp
GaussianQuadrature::GaussianQuadrature(
    int n,
    const GaussianOrthogonalPolynomial& orthPoly)
    : x_(n), w_(n)
{

    // set-up matrix to compute the roots and the weights
    std::vector<double> e(n-1);

    for (int i=1; i < n; ++i) {
        x_[i] = orthPoly.alpha(i);
        e[i-1] = std::sqrt(orthPoly.beta(i));
    }
    x_[0] = orthPoly.alpha(0);

    TqrEigenDecomposition tqr(
        x_, e,
        TqrEigenDecomposition::OnlyFirstRowEigenVector,
        TqrEigenDecomposition::Overrelaxation);

    x_ = tqr.eigenvalues();
    const std::vector< std::vector<double> >& ev = tqr.eigenvectors();

    double mu_0 = orthPoly.mu_0();
    for (int i=0; i<n; ++i) {
        w_[i] = mu_0*ev[0][i]*ev[0][i] / orthPoly.w(x_[i]);
    }
}




// ql/pricingengines/vanilla/analytichestonengine.hpp
// ql/pricingengines/vanilla/analytichestonengine.cpp
//
//! analytic Heston-model engine based on Fourier transform

/*! Integration detail:
    Two algebraically equivalent formulations of the complex
    logarithm of the Heston model exist. Gatherals [2005]
    (also Duffie, Pan and Singleton [2000], and Schoutens,
    Simons and Tistaert[2004]) version does not cause
    discoutinuities whereas the original version (e.g. Heston [1993])
    needs some sort of "branch correction" to work properly.
    Gatheral's version does also work with adaptive integration
    routines and should be preferred over the original Heston version.
*/

/*! References:

    Heston, Steven L., 1993. A Closed-Form Solution for Options
    with Stochastic Volatility with Applications to Bond and
    Currency Options.  The review of Financial Studies, Volume 6,
    Issue 2, 327-343.

    A. Sepp, Pricing European-Style Options under Jump Diffusion
    Processes with Stochastic Volatility: Applications of Fourier
    Transform (<http://math.ut.ee/~spartak/papers/stochjumpvols.pdf>)

    R. Lord and C. Kahl, Why the rotation count algorithm works,
    http://papers.ssrn.com/sol3/papers.cfm?abstract_id=921335

    H. Albrecher, P. Mayer, W.Schoutens and J. Tistaert,
    The Little Heston Trap, http://www.schoutens.be/HestonTrap.pdf

    J. Gatheral, The Volatility Surface: A Practitioner's Guide,
    Wiley Finance

    \ingroup vanillaengines

    \test the correctness of the returned value is tested by
          reproducing results available in web/literature
          and comparison with Black pricing.
*/
#define QL_EPSILON             ((std::numeric_limits<double>::epsilon)())
enum ComplexLogFormula { Gatheral, BranchCorrection };
class heston_Fj : public std::unary_function<double, double>
{
private:
    const int j_;
    const double kappa_, theta_, sigma_, v0_;
    const ComplexLogFormula cpxLog_;

    // helper variables
    const double term_;
    const double x_, sx_, dd_;
    const double sigma2_, rsigma_;
    const double t0_;

    // log branch counter
    mutable int  b_;     // log branch counter
    mutable double g_km1_; // imag part of last log value

public:
    heston_Fj(double kappa, double theta, double sigma, double v0,
              double s0, double rho,
              ComplexLogFormula cpxLog,
              double term, double strike, double ratio, int j);

    double operator()(double phi)      const;
};

heston_Fj::heston_Fj(double kappa, double theta, double sigma, double v0,
                     double s0, double rho,
                     ComplexLogFormula cpxLog,
                     double term, double strike, double ratio, int j) :
    j_(j), kappa_(kappa), theta_(theta), sigma_(sigma),
    v0_(v0), cpxLog_(cpxLog), term_(term), x_(std::log(s0)),
    sx_(std::log(strike)), dd_(x_-std::log(ratio)),
    sigma2_(sigma_*sigma_), rsigma_(rho*sigma_),
    t0_(kappa - ((j== 1)? rho*sigma : 0)),
    b_(0), g_km1_(0) { };


double heston_Fj::operator()(double phi) const
{
    const double rpsig(rsigma_*phi);

    const std::complex<double> t1 = t0_+std::complex<double>(0, -rpsig);
    const std::complex<double> d =
        std::sqrt(t1*t1 - sigma2_*phi
                  *std::complex<double>(-phi, (j_== 1)? 1 : -1));
    const std::complex<double> ex = std::exp(-d*term_);
    const std::complex<double> addOnTerm =  0.0;

    if (cpxLog_ == Gatheral) {
        if (phi != 0.0) {
            if (sigma_ > 1e-5) {
                const std::complex<double> p = (t1-d)/(t1+d);
                const std::complex<double> g
                    = std::log((1.0 - p*ex)/(1.0 - p));

                return
                    std::exp(v0_*(t1-d)*(1.0-ex)/(sigma2_*(1.0-ex*p))
                             + (kappa_*theta_)/sigma2_*((t1-d)*term_-2.0*g)
                             + std::complex<double>(0.0, phi*(dd_-sx_))
                             + addOnTerm
                            ).imag()/phi;
            } else {
                const std::complex<double> td = phi/(2.0*t1)
                                                *std::complex<double>(-phi, (j_== 1)? 1 : -1);
                const std::complex<double> p = td*sigma2_/(t1+d);
                const std::complex<double> g = p*(1.0-ex);

                return
                    std::exp(v0_*td*(1.0-ex)/(1.0-p*ex)
                             + (kappa_*theta_)*(td*term_-2.0*g/sigma2_)
                             + std::complex<double>(0.0, phi*(dd_-sx_))
                             + addOnTerm
                            ).imag()/phi;
            }
        } else {
            // use l'Hospital's rule to get lim_{phi->0}
            if (j_ == 1) {
                const double kmr = rsigma_-kappa_;
                if (std::fabs(kmr) > 1e-7) {
                    return dd_-sx_
                           + (std::exp(kmr*term_)*kappa_*theta_
                              -kappa_*theta_*(kmr*term_+1.0) ) / (2*kmr*kmr)
                           - v0_*(1.0-std::exp(kmr*term_)) / (2.0*kmr);
                } else
                    // \kappa = \rho * \sigma
                    return dd_-sx_ + 0.25*kappa_*theta_*term_*term_
                           + 0.5*v0_*term_;
            } else {
                //TODO: check why is p,g unused
                //const std::complex<double> p = (t1-d)/(t1+d);
                //const std::complex<double> g
                //    = std::log((1.0 - p*ex)/(1.0 - p));

                return dd_-sx_
                       - (std::exp(-kappa_*term_)*kappa_*theta_
                          +kappa_*theta_*(kappa_*term_-1.0))/(2*kappa_*kappa_)
                       - v0_*(1.0-std::exp(-kappa_*term_))/(2*kappa_);
            }
        }
    } else if (cpxLog_ == BranchCorrection) {
        const std::complex<double> p  = (t1+d)/(t1 - d);

        // next term: g = std::log((1.0 - p*std::exp(d*term_))/(1.0 - p))
        std::complex<double> g;

        // the exp of the following expression is needed.
        const std::complex<double> e = std::log(p)+d*term_;

        // does it fit to the machine precision?
        if (std::exp(-e.real()) > QL_EPSILON) {
            g = std::log((1.0 - p/ex)/(1.0 - p));
        } else {
            // use a "big phi" approximation
            g = d*term_ + std::log(p/(p - 1.0));

            if (g.imag() > M_PI || g.imag() <= -M_PI) {
                // get back to principal branch of the complex logarithm
                double im = std::fmod(g.imag(), 2*M_PI);
                if (im > M_PI)
                    im -= 2*M_PI;
                else if (im <= -M_PI)
                    im += 2*M_PI;

                g = std::complex<double>(g.real(), im);
            }
        }

        // be careful here as we have to use a log branch correction
        // to deal with the discontinuities of the complex logarithm.
        // the principal branch is not always the correct one.
        // (s. A. Sepp, chapter 4)
        // remark: there is still the change that we miss a branch
        // if the order of the integration is not high enough.
        const double tmp = g.imag() - g_km1_;
        if (tmp <= -M_PI)
            ++b_;
        else if (tmp > M_PI)
            --b_;

        g_km1_ = g.imag();
        g += std::complex<double>(0, 2*b_*M_PI);

        return std::exp(v0_*(t1+d)*(ex-1.0)/(sigma2_*(ex-p))
                        + (kappa_*theta_)/sigma2_*((t1+d)*term_-2.0*g)
                        + std::complex<double>(0,phi*(dd_-sx_))
                        + addOnTerm
                       ).imag()/phi;
    } else {
        assert(false);
    }

    assert(false);
    return 0.0;
}

} // namespace



// own wrapper
double heston_analytic(double rd, double rf, double spot, double strike,
                       double T, double kappa, double theta, double xi,
                       double v0, double rho, int putcall)
{

    const double dividendDiscount = exp(-rf*T);
    const double riskFreeDiscount = exp(-rd*T);
    const double ratio = riskFreeDiscount/dividendDiscount;
    const double sigma=xi;
    const double spotPrice=spot;
    const double strikePrice=strike;
    const ComplexLogFormula cpxLog(Gatheral);
    int evaluations;



    //const double c_inf = std::min(10.0, std::max(0.0001,
    //        std::sqrt(1.0-rho*rho)/sigma)) *(v0 + kappa*theta*T);

    evaluations = 0;
    static GaussLaguerreIntegration   integration(190);
    heston_Fj                  f1(kappa,theta,sigma,v0,spotPrice,rho, cpxLog, T,
                                  strikePrice, ratio, 1);
    heston_Fj                  f2(kappa,theta,sigma,v0,spotPrice,rho, cpxLog, T,
                                  strikePrice, ratio, 2);

    const double p1 = integration(f1)/M_PI;
    evaluations+= integration.order();

    const double p2 = integration(f2)/M_PI;
    evaluations+= integration.order();

    double value=0.0;
    switch (putcall) {
    case 1:
        value = spot*dividendDiscount*(p1+0.5)
                - strikePrice*riskFreeDiscount*(p2+0.5);
        break;
    case -1:
        value = spotPrice*dividendDiscount*(p1-0.5)
                - strikePrice*riskFreeDiscount*(p2-0.5);
        break;
    default:
        assert(false);
    }
    return value;

}



// analytic solver for barrier options based on Lipton
// only works for rd=rf, rho=0
// this is not part of Quantlib, just an extremely simple implementation
double heston_lipton(double rd, double rf, double spot, double strike,
                     double T, double kappa, double theta, double xi,
                     double v0, double rho, int putcall,
                     double barrier_lower, double barrier_upper)
{

    assert(rd==rf);
    assert(rho==0.0);
    assert(putcall==1);

    if(rd!=rf || rho!=0.0 || putcall!=1) {
        return -1.0;
    }



    const int  max_iter=5000;

    double k_n;		// n-th eigenvalue
    double zeta_k;
    double A,B;		// function values, depending on tau, k_n
    double phi_n;

    double forward_spot	= spot*exp((rd-rf)*T);
    double forward_lower	= barrier_lower*exp((rd-rf)*T);
    double forward_upper	= barrier_upper*exp((rd-rf)*T);
    double main_factor	= exp(-rd*T)*sqrt(forward_spot*strike);
    double mu		= -0.5*kappa;

    double threshold = 1E-15;

    double increase=0.0;
    double value=0.0;
    int  sign=-1;

    if(forward_lower<=0.0) forward_lower=strike/1000.0;
    if(forward_upper<=0.0) forward_upper=strike*1000.0;

    for(int n=1; n<max_iter; n++) {
        sign  *= -1;
        k_n    = M_PI * n / log(forward_upper / forward_lower);
        zeta_k = 0.5 * sqrt(xi*xi * (k_n*k_n + 0.25) + kappa*kappa);
        A      =	-kappa * theta * (mu + zeta_k) * T
                    -kappa * theta *
                    log((-mu + zeta_k + (mu + zeta_k) * exp(-2.0 * zeta_k * T)) /
                        (2.0 * zeta_k));
        B      =	(xi*xi * (k_n*k_n + 0.25) * (1.0 - exp(-2.0 * zeta_k * T))) /
                    (4.0 * (-mu + zeta_k + (zeta_k + mu) * exp(-2.0 * zeta_k * T)));

        phi_n  =	(2.0 * (sign * k_n * (sqrt(forward_upper / strike) -
                                          sqrt(strike / forward_upper)) +
                            sin(k_n * log(forward_lower / strike)))) /
                    ((k_n*k_n + 0.25) * log(forward_upper / forward_lower));

        increase=	main_factor * exp(2.0 * (A - B * v0) / (xi*xi)) * phi_n *
                    sin(k_n * log(forward_spot/forward_lower));


        value += increase;
        if(fabs(increase)<threshold) return value;
    }
    return value;
}

} // namespace fin
