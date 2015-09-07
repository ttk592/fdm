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


#include <cstdio>

#include <boost/numeric/ublas/banded.hpp>

//#include <fdm/fdmparameter.hpp>
//#include <fdm/functors.hpp>
//#include <fdm/grid.hpp>
//#include <fdm/matrix.hpp>
//#include <fdm/solve_pde.hpp>
//#include <fdm/svector.hpp>



// LR-Decomposition of a band matrix, we work on the matrix
// returns -1 on failure
int LRDecompose(boost::numeric::ublas::banded_matrix<double>& A)
{
    int		i_max,j_max;
    const double	eps = 1e-300;
    double		x;

    assert(A.size1()==A.size2());

    // preconditioning: normalize column i so that a_ii=1
    // we omit this here because we would need to save the
    // normalising factors, so need another vector
    // caution!

    // Gauss LR-Decomposition
    for(int k=0; k<(int)A.size1(); k++) {
        i_max=std::min(A.size1()-1,k+A.lower());	// n_lower is not a mistake!
        for(int i=k+1; i<=i_max; i++) {
            if(fabs(A(k,k))<eps) {
                printf("LRDecompose(): zero (%e) pivot element at line %i\n",A(k,k),k);
                return -1;
            }
            x=-A(i,k)/A(k,k);
            A(i,k)=-x;				// assembly part of L
            j_max=std::min(A.size1()-1,k+A.upper());
            for(int j=k+1; j<=j_max; j++) {
                A(i,j)=A(i,j)+x*A(k,j);		// assembly part of R
            }
        }
    }
    return 0;
}

// solves Ly=b
void LSolve(const boost::numeric::ublas::banded_matrix<double>& A,
            boost::numeric::ublas::vector<double>& x,
            const boost::numeric::ublas::vector<double>& b)
{
    int	j_start;
    double	sum;

    for(int i=0; i<(int)A.size1(); i++) {
        sum=0.0;
        j_start=std::max(0,i-(int)A.lower());
        for(int j=j_start; j<i; j++) sum += A(i,j)*x(j);
        x(i)=b(i) - sum;
    }
}

// solves Rx=y
void RSolve(const boost::numeric::ublas::banded_matrix<double>& A,
            boost::numeric::ublas::vector<double>& x,
            const boost::numeric::ublas::vector<double>& b)
{
    assert( (A.size1()==A.size2()) );
    assert( (A.size1()==x.size()) && (A.size1()==b.size()) );
    int	j_stop;
    double	sum;
    for(int i=A.size1()-1; i>=0; i--) {
        sum=0.0;
        j_stop=std::min(A.size1()-1,i+A.upper());
        for(int j=i+1; j<=j_stop; j++) {
            sum += A(i,j)*x(j);	// warning, here we are getting underflows!!!
        }
        x(i)=( b(i) - sum ) / A(i,i);
    }
}

// solves Ax=b if A is already LR decomposed
void LRSolve(const boost::numeric::ublas::banded_matrix<double>& A,
             boost::numeric::ublas::vector<double>& x,
             const boost::numeric::ublas::vector<double>& b)
{
    assert( (A.size1()==A.size2()) );
    assert( (A.size1()==x.size()) && (A.size1()==b.size()) );

    boost::numeric::ublas::vector<double>	y(A.size1());
    LSolve(A,y,b);
    RSolve(A,x,y);
}
