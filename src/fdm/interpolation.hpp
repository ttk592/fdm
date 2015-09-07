/*
 * ---------------------------------------------------------------------
 * Copyright (C) 2015 Tino Kluge (ttk448 at gmail.com)
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

#ifndef MF_FDM_INTERPOLATION_HPP
#define MF_FDM_INTERPOLATION_HPP

#include <cstdio>
#include <cassert>
#include <vector>
#include <boost/array.hpp>

#include <fdm/common.hpp>
#include <fdm/algorithm.hpp>
#include <fdm/matrix.hpp>

#include <boost/numeric/ublas/banded.hpp>

namespace fdm
{


// TODO:
//  - unify both spline interpolation schemes
//  - make them more efficient, ie no unnecessary copies of std::vector
//  - note, though, that Hermite splines may be used more frequently
//    for single point interpolation, where standard cubic splines
//    are used more for grid-->grid interpolation
//  - return derivatives as well as interpolated function values?


// cubic Hermite spline interpolation
// fixes function values, and first derivatives at given two points
// this interpolation is only C1 but only requires two nodes (as well
// as immediately adjacent nodes for derivatives) for polynomial
// coefficients to be determined
// this is in contrast to the standard cubic splines which are C2 but
// require all nodes to determin polynomial coefficients in each interval
class Hspline
{
private:
    std::vector<double> m_x, m_y;
    std::vector<double> m_deriv1;
    bool    m_calculated;
public:
    Hspline()
    {
        m_calculated=false;
    }
    void resize(size_t n)
    {
        m_calculated=false;
        m_x.resize(n);
        m_y.resize(n);
    }
    size_t size() const
    {
        return m_x.size();
    }
    double x(size_t i) const
    {
        assert(i<m_x.size());
        return m_x[i];
    }
    double& x(size_t i)
    {
        assert(i<m_x.size());
#ifndef NDEBUG
        m_calculated=false;
#endif
        return m_x[i];
    }
    double y(size_t i) const
    {
        assert(i<m_y.size());
        return m_y[i];
    }
    double& y(size_t i)
    {
        assert(i<m_y.size());
#ifndef NDEBUG
        m_calculated=false;
#endif
        return m_y[i];
    }

    // calculate all required fields for later evaluation
    void calculate()
    {
        // determining first derivatives at all grid points
        // use the same factors as in GridC2::calculate()
        assert_strictly_increasing(m_x);
        assert(m_x.size()==m_y.size());
        assert(m_x.size()>=3);
        const size_t n=m_x.size();
        double dx1, dx2;
        m_deriv1.resize(n);
        // inner grid points
        for(size_t i=1; i<n-1; i++) {
            dx1=m_x[i]-m_x[i-1];
            dx2=m_x[i+1]-m_x[i];

            m_deriv1[i] = -dx2/((dx1+dx2)*dx1) * m_y[i-1]
                          + (dx2-dx1)/(dx1*dx2)  * m_y[i]
                          + dx1/((dx1+dx2)*dx2)  * m_y[i+1];
        }
        // outer grid points
        dx1=m_x[1]-m_x[0];
        dx2=m_x[2]-m_x[1];
        m_deriv1[0]     = -(2.0*dx1+dx2)/(dx1*(dx1+dx2)) * m_y[0]
                          + (dx1+dx2)/(dx1*dx2)           * m_y[1]
                          - dx1/(dx2*(dx1+dx2))           * m_y[2];

        dx1=m_x[n-2]-m_x[n-3];
        dx2=m_x[n-1]-m_x[n-2];
        m_deriv1[n-1]   = dx2/(dx1*(dx1+dx2))           * m_y[n-3]
                          - (dx1+dx2)/(dx1*dx2)           * m_y[n-2]
                          + (2.0*dx2+dx1)/(dx2*(dx1+dx2)) * m_y[n-1];

        m_calculated=true;
    }
    void set_points(const std::vector<double>& x,
                    const std::vector<double>& y)
    {
        m_x=x;
        m_y=y;
        calculate();
    }

    // function evaluation at x
    // can also provide initial_guess, which is index of the grid point
    // closest to x
    double operator() (double x, int initial_guess=-1) const
    {
        assert(m_calculated==true);
        const size_t n=m_x.size();
        // find the closest point m_x[k] < x, k=0 even if x<m_x[0]
        const int k=std::max(0, fastfind(m_x,x,initial_guess));
        assert(k<(int)n);

        double interpol;
        if(x<=m_x[0]) {
            // linear extrapolation to the left
            interpol = m_y[0] + m_deriv1[0]*(x-m_x[k]);
        } else if(x>=m_x[n-1]) {
            // linear extrapolation to the right
            interpol=m_y[n-1] + m_deriv1[n-1]*(x-m_x[k]);
        } else {
            // interpolation
            interpol=hermite(x,k);
        }
        return interpol;
    }

    // hermite interpolation at x, if grid node is known to be k
    // do not call directly unless index k is known
    inline double hermite(double x, int k) const
    {
        assert(k>=0);
        assert(k<(int)m_x.size()-1);
        assert(m_calculated==true);

        // TODO: for performance consider moving some of this to calculate()
        double dx=m_x[k+1]-m_x[k];
        double z=(x-m_x[k])/dx;
        double z2=z*z;              // z^2
        double z3=z2*z;             // z^3
        // hermite polynomials
        double h00 =  2.0*z3 - 3.0*z2 + 1.0;    // h(0)=1
        double h01 = -2.0*z3 + 3.0*z2;          // h(1)=1
        double h10 =      z3 - 2.0*z2 + z;      // h'(0)=1
        double h11 =      z3 -     z2;          // h'(1)=1

        return m_y[k]*h00 + m_y[k+1]*h01
               + m_deriv1[k]*dx*h10 + m_deriv1[k+1]*dx*h11;
    }

};


// cubic spline interpolation
// generates a C2 function everywhere, but it requires all input points to
// perform the interpolation at any given point; it is therfore not
// efficient enough for higher dimensional interpolation, but very good
// for 1D grid --> grid interpoaltion
// linear extrapolation if default (zero curvature) boundary conditions
// are specified
class Cspline
{
public:
    enum bd_type {
        first_deriv = 1,
        second_deriv = 2
    };

private:
    std::vector<double> m_x,m_y;            // x,y coordinates of points
    // interpolation parameters
    // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
    boost::numeric::ublas::vector<double> m_a,m_b,m_c;  // coefficients
    double  m_b0, m_c0;                     // for left extrapol
    bd_type m_left, m_right;
    double  m_left_value, m_right_value;
    bool    m_force_linear_extrapolation;

public:
    // set default boundary condition to be zero curvature at both ends
    Cspline(): m_left(second_deriv), m_right(second_deriv),
        m_left_value(0.0), m_right_value(0.0),
        m_force_linear_extrapolation(false)
    {
        ;
    }

    // optional, but if called it has to come be before set_points()
    void set_boundary(bd_type left, double left_value,
                      bd_type right, double right_value,
                      bool force_linear_extrapolation=false)
    {
        assert(m_x.size()==0);  // set_points() must not have happened yet
        m_left=left;
        m_right=right;
        m_left_value=left_value;
        m_right_value=right_value;
        m_force_linear_extrapolation=force_linear_extrapolation;

    }
    void set_points(const std::vector<double>& x,
                    const std::vector<double>& y, bool cubic_spline=true)
    {
        assert(x.size()==y.size());
        assert(x.size()>2);
        m_x=x;
        m_y=y;
        calculate(cubic_spline);
    }
    void resize(size_t n)
    {
        m_x.resize(n);
        m_y.resize(n);
    }
    void set_point(size_t i, double x, double y)
    {
        assert(i<m_x.size() && i<m_y.size());
        m_x[i]=x;
        m_y[i]=y;
    }
    void calculate(bool cubic_spline=true)
    {
        assert(m_x.size()==m_y.size());
        int   n=m_x.size();
        // TODO: maybe sort x and y, rather than returning an error
        for(int i=0; i<n-1; i++) {
            assert(m_x[i]<m_x[i+1]);
        }

        if(cubic_spline==true) { // cubic spline interpolation
            // setting up the matrix and right hand side of the equation system
            // for the parameters b[]
            boost::numeric::ublas::banded_matrix<double> A(n,n,1,1);
            boost::numeric::ublas::vector<double> rhs(n);
            for(int i=1; i<n-1; i++) {
                A(i,i-1)=1.0/3.0*(m_x[i]-m_x[i-1]);
                A(i,i)=2.0/3.0*(m_x[i+1]-m_x[i-1]);
                A(i,i+1)=1.0/3.0*(m_x[i+1]-m_x[i]);
                rhs[i]=(m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i])
                       - (m_y[i]-m_y[i-1])/(m_x[i]-m_x[i-1]);
            }
            // boundary conditions
            if(m_left == Cspline::second_deriv) {
                // 2*b[0] = f''
                A(0,0)=2.0;
                A(0,1)=0.0;
                rhs[0]=m_left_value;
            } else if(m_left == Cspline::first_deriv) {
                // c[0] = f', needs to be re-expressed in terms of b:
                // (2b[0]+b[1])(m_x[1]-m_x[0])
                // = 3 ((m_y[1]-m_y[0])/(m_x[1]-m_x[0]) - f')
                A(0,0)=2.0*(m_x[1]-m_x[0]);
                A(0,1)=1.0*(m_x[1]-m_x[0]);
                rhs[0]=3.0*((m_y[1]-m_y[0])/(m_x[1]-m_x[0])-m_left_value);
            } else {
                assert(false);
            }
            if(m_right == Cspline::second_deriv) {
                // 2*b[n-1] = f''
                A(n-1,n-1)=2.0;
                A(n-1,n-2)=0.0;
                rhs[n-1]=m_right_value;
            } else if(m_right == Cspline::first_deriv) {
                // c[n-1] = f', needs to be re-expressed in terms of b:
                // (b[n-2]+2b[n-1])(m_x[n-1]-m_x[n-2])
                // = 3 (f' - (m_y[n-1]-m_y[n-2])/(m_x[n-1]-m_x[n-2]))
                A(n-1,n-1)=2.0*(m_x[n-1]-m_x[n-2]);
                A(n-1,n-2)=1.0*(m_x[n-1]-m_x[n-2]);
                rhs[n-1]=3.0*(m_right_value-(m_y[n-1]-m_y[n-2])
                              /(m_x[n-1]-m_x[n-2]));
            } else {
                assert(false);
            }

            // solve the equation system to obtain the parameters b[]
            m_b.resize(n);
            LRDecompose(A);
            LRSolve(A,m_b,rhs);

            // calculate parameters a[] and c[] based on b[]
            m_a.resize(n);
            m_c.resize(n);
            for(int i=0; i<n-1; i++) {
                m_a[i]=1.0/3.0*(m_b[i+1]-m_b[i])/(m_x[i+1]-m_x[i]);
                m_c[i]=(m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i])
                       - 1.0/3.0*(2.0*m_b[i]+m_b[i+1])*(m_x[i+1]-m_x[i]);
            }
        } else { // linear interpolation
            m_a.resize(n);
            m_b.resize(n);
            m_c.resize(n);
            for(int i=0; i<n-1; i++) {
                m_a[i]=0.0;
                m_b[i]=0.0;
                m_c[i]=(m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i]);
            }
        }

        // for left extrapolation coefficients
        m_b0 = (m_force_linear_extrapolation==false) ? m_b[0] : 0.0;
        m_c0 = m_c[0];

        // for the right extrapolation coefficients
        // f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
        double h=m_x[n-1]-m_x[n-2];
        // m_b[n-1] is determined by the boundary condition
        m_a[n-1]=0.0;
        m_c[n-1]=3.0*m_a[n-2]*h*h+2.0*m_b[n-2]*h+m_c[n-2];// = f'_{n-2}(x_{n-1})
        if(m_force_linear_extrapolation==true)
            m_b[n-1]=0.0;


    }


    double operator() (double x) const
    {
        size_t n=m_x.size();
        // find the closest point m_x[idx] <= x, idx=0 even if x<m_x[0]
        const int idx=std::max(0, std_find(m_x,x));

        double h=x-m_x[idx];
        double interpol;
        if(x<m_x[0]) {
            // extrapolation to the left
            interpol=(m_b0*h + m_c0)*h + m_y[0];
        } else if(x>m_x[n-1]) {
            // extrapolation to the right
            interpol=(m_b[n-1]*h + m_c[n-1])*h + m_y[n-1];
        } else {
            // interpolation
            interpol=((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
        }
        return interpol;
    }
};



} // namespace fdm

#endif /* MF_FDM_INTERPOLATION_HPP */
