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


#ifndef MF_FDM_MATRIX_H
#define MF_FDM_MATRIX_H

#include <boost/numeric/ublas/banded.hpp>

// LR-Decomposition of a band matrix, we work on the matrix
// returns -1 on failure
int LRDecompose(boost::numeric::ublas::banded_matrix<double>& A);

// solves Ly=b
void LSolve(const boost::numeric::ublas::banded_matrix<double>& A,
            boost::numeric::ublas::vector<double>& x,
            const boost::numeric::ublas::vector<double>& b);
// solves Rx=y
void RSolve(const boost::numeric::ublas::banded_matrix<double>& A,
            boost::numeric::ublas::vector<double>& x,
            const boost::numeric::ublas::vector<double>& b);

// solves Ax=b if A is already LR decomposed
void LRSolve(const boost::numeric::ublas::banded_matrix<double>& A,
             boost::numeric::ublas::vector<double>& x,
             const boost::numeric::ublas::vector<double>& b);

// implementation part (could be moved to cpp)
#include <fdm/matrix.inl>

#endif // MF_FDM_MATRIX_H
