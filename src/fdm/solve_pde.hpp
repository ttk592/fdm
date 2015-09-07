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


#ifndef MF_FDM_SOLVE_PDE_H
#define MF_FDM_SOLVE_PDE_H

#include <cstdio>
#include <cassert>
#include <vector>
#include <string>

#include <boost/timer.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp>

#include <fdm/svector.hpp>
#include <fdm/grid.hpp>
#include <fdm/matrix.hpp>
#include <fdm/pde.hpp>


namespace fdm
{


// namespace for enum's
namespace type
{

enum method {
    cn          = 0,        // Crank-Nicolson
    adi         = 1,        // ADI
    pred_corr   = 2,        // Predictor Corrector Method like Craig&Sneyd
    cn_1d       = 3,        // Crank-Nicolson in 1D
};

} // namespace type




// structure for storing timing information
class TimingInfo
{
private:
    std::array<double, 3>       m_timing;
    std::array<boost::timer, 3> m_t;
public:

    enum type {
        assembly    = 0,
        solve       = 1,
        total       = 2
    };

    std::size_t size() const
    {
        assert(m_timing.size()==m_t.size());
        return m_timing.size();
    }

    TimingInfo()
    {
        for(size_t i=0; i<size(); i++) {
            m_timing[i]=0.0;
        }
    }

    void start(type w)
    {
        std::size_t k = static_cast<std::size_t>(w);
        assert(k<size());
        m_t[k].restart();
    }
    void stop(type w)
    {
        std::size_t k = static_cast<std::size_t>(w);
        assert(k<size());
        m_timing[k]+=m_t[k].elapsed();
    }
    void reset(type w)
    {
        std::size_t k = static_cast<std::size_t>(w);
        assert(k<size());
        m_timing[k]=0.0;
    }
    double elapsed(type w) const
    {
        std::size_t k = static_cast<std::size_t>(w);
        assert(k<size());
        return m_timing[k];
    }

    // adding all times together
    TimingInfo operator+(const TimingInfo& a)
    {
        TimingInfo b;
        for(size_t i=0; i<size(); i++) {
            b.m_timing[i]=m_timing[i]+a.m_timing[i];
        }
        return b;
    }

};


// ---------------------------------------------------------------------
// C r a n k - N i c o l s o n    M e t h o d   1 D
// ---------------------------------------------------------------------
// we only implement a Crank-Nicolson in 1d as otherwise we would
// need a direct sparse matrix solver like SuperLU
// (requires redistribution of copyright and other notices) or
// Umfpack (LGPL)
//
// Crank-Nicolson one time step
// returns timing information
template <std::size_t GridDim>
TimingInfo SolveCN1dStep(PDE<GridDim>& UNUSED(pde),
                         GridFunctionBase<GridDim>& UNUSED(u),
                         const Boundary<GridDim>& UNUSED(boundary),
                         double UNUSED(theta), double UNUSED(dt),
                         double UNUSED(t), bool UNUSED(recalc))
{
    //TODO static assert would be better here
    throw std::runtime_error("SolveCN1dStep only supported for 1d grids.");
    return TimingInfo();
}

// 1d-specialisation of Crank-Nicolson
// (note, template specialisation are not automatically inline)
template <>
inline TimingInfo SolveCN1dStep<1>(PDE<1>& pde,
                                   GridFunctionBase<1>& u,
                                   const Boundary<1>& boundary,
                                   double theta, double dt, double t,
                                   bool recalc)
{

    TimingInfo      timing;
    timing.start(TimingInfo::total);

    const int	    vec_dim=u.num_grid_points();

    int		vec_ind;
    int		bdim;
    int		istart;

    // we define the matrices to be static so it can be reused
    // at a later function call (making this function not thread safe)
    static boost::numeric::ublas::banded_matrix<double>	A;

    // checks if band matrices have been initialised
    if((int)A.size1()!=vec_dim) {
        // sets up a bandmatrix of size NxN with a diagonal plus
        // 2 lower and 2 upper diagonals (remainder all zero)
        A.resize(vec_dim,vec_dim,2,2,false);
    }

    boost::numeric::ublas::vector<double>	b(vec_dim);
    boost::numeric::ublas::vector<double>	sol(vec_dim);

    SVector<int, 1>		idx;		// multi-index
    SVector<double, 1>	x;		// coordinates

    // set time t, at which we read pde coefficients
    pde.set_time(t+theta*dt);

    // assembly matrix A and vector b
    // ------------------------------
    timing.start(TimingInfo::assembly);
    if(recalc) A.clear();		// set all entries to zero
    u.iterator.begin();
    while(u.iterator.isend()==false) {
        idx=u.idx();
        x=u.coord();
        vec_ind=idx[0];
        bdim=u.iterator.check_boundary();
        istart=u.deriv_index_start(0,idx[0]);

        // if inner point or if boundary and free condition
        if( (bdim<0) || (boundary.type(u.iterator)==BoundaryType::Free) ) {
            // discretise pde
            pde.set_space(x);
            // assembly vector: b = (I+ (1-theta)dt L) u_{i-1}
            b(vec_ind) = u.value() + (1.0-theta) * dt * (
                             pde.deriv2(0,0)*u.f_2(0,0,idx)
                             + pde.deriv1(0)*u.f_1(0,idx)
                             + pde.deriv0() * u.value(idx)
                         );

            // assembly matrix:
            // (I - theta dt L) u_{i} = b
            if(recalc) {
                for(int i=istart; i<istart+3; i++) {
                    A(vec_ind,vec_ind+i) = - theta * dt * (
                                               pde.deriv2(0,0)*u.fac_deriv2(0,idx,i)
                                               + pde.deriv1(0)*u.fac_deriv1(0,idx,i)
                                           );
                }
                A(vec_ind,vec_ind) += 1.0 - theta*dt*pde.deriv0();
            }

        } else {
            // boundary conditions
            b(vec_ind)=boundary.value(u.iterator,x,t+theta*dt);
            // assemble matrix according to the type of condition

            if(recalc) switch(boundary.type(u.iterator)) {
                case BoundaryType::Value:	// fixed values given
                    A(vec_ind,vec_ind)=1.0;
                    break;
                case BoundaryType::FirstDeriv:	// first derivative given
                    for(int k=istart; k<istart+3; k++) {
                        A(vec_ind,vec_ind+k) = u.fac_deriv1(0,idx,k);
                    }
                    break;
                case BoundaryType::SecondDeriv:	// second derivative given
                    for(int k=istart; k<istart+3; k++) {
                        A(vec_ind,vec_ind+k) = u.fac_deriv2(bdim,idx,k);
                    }
                    break;
                default:
                    abort();
                }
        } // end of assembly all for point idx

        u.iterator++;
    } // end of assembly all
    timing.stop(TimingInfo::assembly);


    // solve equation A x = b
    // ----------------------
    timing.start(TimingInfo::solve);
    if(recalc) LRDecompose(A);
    LRSolve(A,sol,b);
    timing.stop(TimingInfo::solve);

    // writing the solution back into u
    // --------------------------------
    timing.start(TimingInfo::assembly);
    u.iterator.begin();
    while(u.iterator.isend()==false) {
        idx=u.idx();
        vec_ind=idx[0];
        u.value()=sol(vec_ind);
        u.iterator++;
    }
    timing.stop(TimingInfo::assembly);
    timing.stop(TimingInfo::total);

    return timing;
}







// ---------------------------------------------------------------------
// P r e d i c t o r    C o r r e c t o r    M e t h o d
// ---------------------------------------------------------------------
// Predictor: calculates u_new = u_old + dt L u_old + lambda X (u_estim-u_old)
// ie fully explicit scheme
// X depends on the scheme used, L is the linear space operator


// Predictor at an inner point (or boundary where pde discretised)
// this is called as one of the most inner loops of the assembly routines
// and needs to be as efficient as possible
template <std::size_t GridDim>
double PredictorCalculateInnerPoint(const GridFunctionBase<GridDim>& u_old,
                                    const GridFunctionBase<GridDim>& u_estim,
                                    const PDE<GridDim>& pde,
                                    const SVector<int, GridDim>& idx,
                                    double theta, double dt, bool use_estim)
{
    double value=0.0;

    // calculate L u_old for second derivatives only
    for(int dim1=0; dim1<pde.dim(); dim1++) {
        for(int dim2=dim1; dim2<pde.dim(); dim2++) {
            value += pde.deriv2(dim1,dim2)*u_old.f_2(dim1,dim2,idx);
        }
    }
    // calculate L u_old for first derivatives only
    for(int dim1=0; dim1<pde.dim(); dim1++) {
        value += pde.deriv1(dim1)*u_old.f_1(dim1,idx);
    }
    // calculate L u_old for 0th derivative
    value += pde.deriv0() * u_old.value(idx);

    // calculate L u_old for the source term f
    value += pde.f();


    // we add a further correction: L full space matrix, M only mixed deriv
    // Craig & Sneyd:		theta dt M (u_estim-u_old)
    // Modified Craig & Sneyd: 	above + (0.5-theta) dt L (u_estim-u_old)
    // Hundsdorfer & Verwer		0.5 dt L (u_estim-u_old)
    if(use_estim) {
        // assembly second derivatives
        for(int dim1=0; dim1<pde.dim(); dim1++) {
            for(int dim2=dim1+1; dim2<pde.dim(); dim2++) {
                value += theta*pde.deriv2(dim1,dim2)*
                         ( u_estim.f_2(dim1,dim2,idx) - u_old.f_2(dim1,dim2,idx) );
            }
        }
    }

    return u_old.value(idx)+dt*value;
}


// Predictor at a boundary point
// in an explicit scheme the values at the boundary are fully determined
// by the value of the grid's inner point and not the previous time slice
// so we only need u_new as input
// the index point is taken from u_new.iterator
// (do not use for free boundary condition (type 3))
template <std::size_t GridDim>
double PredictorCalculateBoundaryPoint(const GridFunctionBase<GridDim>& u_new,
                                       const Boundary<GridDim>& boundary, double t)
{

    SVector<int, GridDim>	idx=u_new.idx();
    SVector<int, GridDim>	idx2(idx);
    SVector<double,GridDim>	x=u_new.coord();
    int			dim=u_new.iterator.boundary();
    int			type=boundary.type(u_new.iterator);
    double		bvalue=boundary.value(u_new.iterator,x,t);
    int			istart=u_new.deriv_index_start(dim,idx[dim]);
    double	value=0.0;

    assert( (type==0) || (type==1) || (type==2) );
    switch(type) {
    case BoundaryType::Value:		// fixed values given
        value=bvalue;
        break;
    case BoundaryType::FirstDeriv:	// first derivative given
        for(int k=istart; k<istart+2; k++) {
            if(k!=0) {
                idx2[dim]=idx[dim]+k;
                value+=u_new.fac_deriv1(dim,idx,k)*u_new.value(idx2);
            }
        }
        value=bvalue-value/u_new.fac_deriv1(dim,idx,0);
        break;
    case BoundaryType::SecondDeriv:	// second derivative given
        for(int k=istart; k<istart+2; k++) {
            if(k!=0) {
                idx2[dim]=idx[dim]+k;
                value+=u_new.fac_deriv2(dim,idx,k)*u_new.value(idx2);
            }
        }
        value=bvalue-value/u_new.fac_deriv2(dim,idx,0);
        break;
    }
    return value;
}


// Predictor: calculates u_new = u_old + dt L u_old + lambda X (u_estim-u_old)
// ie fully explicit scheme (lambda=0 if use_estim=false, otherwise 1)
// X depends on the scheme used, L is the linear space operator
template <std::size_t GridDim>
TimingInfo PredictorStep(PDE<GridDim>& pde,
                         GridFunction<GridDim>& u_new,
                         const GridFunctionBase<GridDim>& u_old,
                         const GridFunction<GridDim>& u_estim,
                         const Boundary<GridDim>& boundary,
                         double theta, double dt, double t,
                         bool use_estim, TimingInfo& timing)
{

    timing.start(TimingInfo::total);
    timing.start(TimingInfo::assembly);

    SVector<double,GridDim>	x;
    SVector<int,GridDim>	idx;

    // set time t, at which we read pde coefficients
    pde.set_time(t+theta*dt);

    // inner grid points
    u_new.iterator.begin_inner();
    while(u_new.iterator.isend()==false) {
        idx=u_new.idx();
        x=u_new.coord();

        // set coordinates at which PDE coefficients are obtained
        pde.set_space(x);

        // calculate the new value at the index point
        u_new.value() = PredictorCalculateInnerPoint(u_old,u_estim,pde,
                        idx,theta,dt,use_estim);
        u_new.iterator.incr_inner();
    }

    // boundary points
    u_new.iterator.begin();
    while(u_new.iterator.isend()==false) {
        idx=u_new.idx();

        if(boundary.type(u_new.iterator)==BoundaryType::Free) {
            // we treat it as an inner point
            x=u_new.coord();
            pde.set_space(x);
            u_new.value() = PredictorCalculateInnerPoint(u_old,u_estim,pde,
                            idx,theta,dt,use_estim);

        } else {
            // normal boundary condition
            // TODO check whether it should also be t+theta*dt
            u_new.value() = PredictorCalculateBoundaryPoint(u_new,boundary,t);
        }
        u_new.iterator.incr_bound();
    }
    timing.stop(TimingInfo::assembly);
    timing.stop(TimingInfo::total);
    return timing;
}




// calculates for each i=0..GridDim-1:
// (I-theta dt D_i) u = u - theta dt D_i u_old   (overriding u in each step)
// returns timing information
template <std::size_t GridDim>
void CorrectorStep(PDE<GridDim>& pde,
                   GridFunctionBase<GridDim>& u,
                   const GridFunctionBase<GridDim>& u_old,
                   const Boundary<GridDim>& boundary,
                   double theta, double dt, double t,
                   bool recalc, TimingInfo& timing)
{

    timing.start(TimingInfo::total);

    const int	vec_dim=u.num_grid_points();

    int		vec_ind;
    int		bdim;
    int		istart;

    // we define the matrices to be static so it can be reused
    // at a later function call
    static boost::numeric::ublas::banded_matrix<double>	M[GridDim];

    // checks if band matrices have been initialised
    if((int)M[0].size1()!=vec_dim) {
        for(int dim=0; dim<pde.dim(); dim++) {
            // sets up a bandmatrix of size NxN with a diagonal plus
            // 2 lower and 2 upper diagonals (remainder all zero)
            M[dim].resize(vec_dim,vec_dim,2,2,false);
        }
    }

    boost::numeric::ublas::vector<double>	b(vec_dim);
    boost::numeric::ublas::vector<double>	sol(vec_dim);

    SVector<int, GridDim>		idx;		// multi-index
    SVector<double, GridDim>	x;		// coordinates

    // set time t, at which we read pde coefficients
    pde.set_time(t+theta*dt);

    // going through all dimensions
    // ----------------------------
    for(int adi_dim=0; adi_dim<pde.dim(); adi_dim++) {
        boost::numeric::ublas::banded_matrix<double>	&A=M[adi_dim];
        if(recalc) A.clear();		// set all entries to zero

        // assembly matrix A and vector b
        // ------------------------------
        timing.start(TimingInfo::assembly);
        u.iterator.begin();
        while(u.iterator.isend()==false) {
            idx=u.idx();
            x=u.coord();
            vec_ind=u.iterator.index_embed(idx,adi_dim);
            bdim=u.iterator.check_boundary();
            istart=u.deriv_index_start(adi_dim,idx[adi_dim]);

            if( (bdim<0) || (boundary.type(u.iterator)==BoundaryType::Free) ) {
                // discretise pde
                pde.set_space(x);
                // assembly vector: b = u_{i-1} - theta dt D_i u_old
                b(vec_ind) = u.value() - theta * dt * (
                                 pde.deriv2(adi_dim,adi_dim)*u_old.f_2(adi_dim,adi_dim,idx)
                                 + pde.deriv1(adi_dim)*u_old.f_1(adi_dim,idx)
                                 + pde.deriv0() * u_old.value(idx)/pde.dim()
                             );

                // assembly matrix:
                // (I - theta dt D_i) u_{i} = b
                if(recalc) {
                    for(int i=istart; i<istart+3; i++) {
                        A(vec_ind,vec_ind+i) = - theta * dt * (
                                                   pde.deriv2(adi_dim,adi_dim)*u.fac_deriv2(adi_dim,idx,i)
                                                   + pde.deriv1(adi_dim)*u.fac_deriv1(adi_dim,idx,i)
                                               );
                    }
                    A(vec_ind,vec_ind) += 1.0 - theta*dt*pde.deriv0()/pde.dim();
                }

            } else {
                // boundary conditions
                // TODO should it be t+theta*dt
                b(vec_ind)=boundary.value(u.iterator,x,t);
                // assemble matrix according to the type of condition
                switch(boundary.type(u.iterator)) {
                case BoundaryType::Value:	// fixed values given
                    if(recalc) A(vec_ind,vec_ind)=1.0;
                    break;
                case BoundaryType::FirstDeriv:	// first derivative given
                    // only if derivative is in the ADI direction
                    if(bdim==adi_dim) {
                        if(recalc)
                            for(int k=istart; k<istart+3; k++) {
                                A(vec_ind,vec_ind+k) = u.fac_deriv1(bdim,idx,k);
                            }
                    } else {
                        // assign value of previous iteration (improve!)
                        if(recalc) A(vec_ind,vec_ind)=1.0;
                        b(vec_ind)=u.value(idx);
                    }
                    break;
                case BoundaryType::SecondDeriv:	// second derivative given
                    // only if derivative is in the ADI direction
                    if(bdim==adi_dim) {
                        if(recalc)
                            for(int k=istart; k<istart+3; k++) {
                                A(vec_ind,vec_ind+k) = u.fac_deriv2(bdim,idx,k);
                            }
                    } else {
                        // assign value of previous iteration (improve!)
                        if(recalc) A(vec_ind,vec_ind)=1.0;
                        b(vec_ind)=u.value(idx);
                    }
                    break;
                default:
                    abort();
                }
            } // end of assembly all for point idx

            u.iterator++;
        } // end of assembly all
        timing.stop(TimingInfo::assembly);


        // solve equation A x = b
        // ----------------------
        timing.start(TimingInfo::solve);
        if(recalc) LRDecompose(A);
        LRSolve(A,sol,b);
        timing.stop(TimingInfo::solve);

        // writing the solution back into u
        // --------------------------------
        timing.start(TimingInfo::assembly);
        u.iterator.begin();
        while(u.iterator.isend()==false) {
            idx=u.idx();
            vec_ind=u.iterator.index_embed(idx,adi_dim);
            u.value()=sol(vec_ind);
            u.iterator++;
        }
        timing.stop(TimingInfo::assembly);

    } // end of going through all adi dimensions

    timing.stop(TimingInfo::total);
}


// ------------------------------------------------------------
// Predictor Corrector splitting method
// ------------------------------------------------------------
// u    ... spacial function at time t, is taken as input value
//          and is overwritten with calculated new function at t+dt
// t    ... old time (time of input u)
// t+dt ... new time (time of output u)
template <std::size_t GridDim>
TimingInfo SolvePredCorStep(PDE<GridDim>& pde,
                            GridFunctionBase<GridDim>& u,
                            const Boundary<GridDim>& boundary,
                            double theta, double dt, double t,
                            bool recalc)
{

    TimingInfo	timing1,timing2,timing3,timing4;
    GridFunction<GridDim>	u_new(u), u_estim(u);

    // (I - dt theta D) u_new = (I - dt theta D) u_old + dt L u_old + dt f
    PredictorStep(pde,u_new,u,u_estim,boundary,theta,dt,t,false,timing1);
    CorrectorStep(pde,u_new,u,boundary,theta,dt,t,recalc,timing2);
    u_estim.assign_values(u_new);

    // (I - dt theta D) u_new = (I - dt theta D) u_old + dt L u_old + dt f
    //				+ dt theta  M (u_estim-u_old)
    PredictorStep(pde,u_new,u,u_estim,boundary,theta,dt,t,true,timing3);
    CorrectorStep(pde,u_new,u,boundary,theta,dt,t,false,timing4);
    u.assign_values(u_new);

    return timing1+timing2+timing3+timing4;
}



// a crude output function, TODO: retire
// it displays a few timing information on screen and can output
// the whole function u in a file
template <std::size_t GridDim>
void fdm_print(const std::string& method, GridFunctionBase<GridDim>& u,
               const SVector<double,GridDim>& x, double t, int tn,
               const TimingInfo& timing,
               bool out_screen, bool out_file)
{

    char	filename[30];
    SVector<int, GridDim>	idx=u.find_nearest(x);

    // output on screen
    if(out_screen) {
        printf("%s: t=%5.3f: u(x)=%f, u(", method.c_str(), t, u(x));
        for(int i=0; i<u.dim(); i++) {
            printf("%i",idx[i]);
            if(i<u.dim()-1) printf(",");
        }
        printf(")=%f, stime=%.2f  atime=%.2f\n", u.value(idx),
               timing.elapsed(TimingInfo::solve),
               timing.elapsed(TimingInfo::assembly));
    }

    // output into file
    if(out_file) {
        sprintf(filename,"out/fdm-%.5i.dat",tn);
        u.export_points(filename);
    }
}



// TODO: retire  (cannot deal with auxilliary variables)
// high level routine for one time step solver
// takes care of whether to recalc matrix, selects the solver,
// makes sure the first step is an implicit one,
// executes the one step solver, returns solution in specified point x
//
// given u^(k-1) we calculate u^(k), super index indicating time index, ie
// given u at t_(k-1) we calculate u at t_(k)
template <std::size_t GridDim, std::size_t PdeDim, std::size_t PdeAuxDim>
double FDMSolveStep(PDE<PdeDim, PdeAuxDim>& pde,
                    const Boundary<PdeDim>& boundary,
                    GridFunctionBase<GridDim>& u,
                    const Grid<1>& grid_time,
                    const SVector<double,GridDim>& x,
                    int k,
                    type::method method=type::pred_corr,
                    bool out_screen=false, bool out_file=false)
{
    static_assert(PdeDim<=GridDim, "PdeDim cannot be greater than GridDim");
    static_assert(PdeAuxDim==0 || PdeAuxDim==GridDim-PdeDim,
                  "PdeAuxDim incompatible");

    std::string method_str;
    double	theta;		// 0 explicit, 0.5 Crank Nicolson, 1 implicit
    TimingInfo  timing;         // timing of one time step
    bool	recalc;			// recalculation of the matrix
    double	t, dt, dt_diff;
    const double	eps=1E-15;

    assert( (0<k) && (k<grid_time.size(0)) );

    // solving the pde in one time step
    // ----------------------------
    dt = grid_time.coord(0,k)-grid_time.coord(0,k-1);
    t  = grid_time.coord(0,k-1);

    // we start with one implicit step and then use CN
    if(k==1)
        theta=1.0;	// implicit
    else
        theta=0.5;	// Crank Nicolson

    // check if we are forced to recalculate the matrix
    if( (k<=2) || (pde.is_time_indep()==false) ) {
        recalc=true;
    } else {
        // check if dt has changed
        dt_diff=grid_time.coord(0,k-1)-grid_time.coord(0,k-2)-dt;
        if(fabs(dt_diff)<5*eps*grid_time.coord(0,k))	recalc=false;
        else						recalc=true;
    }

    // solving the PDE in one step
    switch(method) {
    case type::pred_corr:	// predictor corrector method (Craig&Snaid)
        timing=SolvePredCorStep(pde,u,boundary,theta,dt,t,recalc);
        method_str="PC";
        break;
    case type::cn_1d:     // 1d Crank-Nicolson
        timing=SolveCN1dStep(pde,u,boundary,theta,dt,t,recalc);
        method_str="CN";
        break;
    default:
        printf("method %i not implemented\n",method);
        abort();
    }
    //out_file=true;
    fdm_print(method_str,u,x,t,k,timing,out_screen,out_file);
    return u(x);
}



template <std::size_t GridDim, std::size_t PdeAuxDim>
TimingInfo FdmSolveStepHelper(GridFunctionBase<GridDim>& u,
                              PDE<GridDim, PdeAuxDim>& pde,
                              const Boundary<GridDim>& boundary,
                              double t, double dt,
                              double theta, bool recalc, type::method method)
{
    TimingInfo  timing;
    switch(method) {
    case type::pred_corr:	// predictor corrector method (Craig&Snaid)
        timing=SolvePredCorStep(pde,u,boundary,theta,dt,t,recalc);
        break;
    case type::cn_1d:     // 1d Crank-Nicolson
        timing=SolveCN1dStep(pde,u,boundary,theta,dt,t,recalc);
        break;
    default:
        printf("method %i not implemented\n",method);
        abort();
    }
    return timing;
}

// new
// can deal with auxilliary variables
template <std::size_t GridDim, std::size_t PdeDim, std::size_t PdeAuxDim>
void FdmSolveStep(GridFunction<GridDim>& u,
                  PDE<PdeDim, PdeAuxDim>& pde,
                  const Boundary<PdeDim>& boundary,
                  double t1, double t2,
                  double theta=0.5,
                  bool recalc=true,
                  type::method method=type::pred_corr)
{
    static_assert(PdeDim<=GridDim, "PdeDim cannot be greater than GridDim");
    static_assert(PdeAuxDim==0 || PdeAuxDim==GridDim-PdeDim,
                  "PdeAuxDim incompatible");

    TimingInfo  timing;     // timing of one time step


    // with auxilliary variables, we first create "views"
    // which only sees the lower-dimensional grid with fixed
    // auxilliary variable, i.e.
    //  for all auxilliary variables a
    //      V(x) = U(x,a), for all x
    //      solve pde on V

    SVector<int, GridDim-PdeDim> shape;
    for(size_t i=0; i<GridDim-PdeDim; i++) {
        shape[i] = u.shape()[i+PdeDim];
    }
    GridIterator<GridDim-PdeDim>    it(shape);
    GridFunctionView<PdeDim> v(u);

    it.begin();
    while(!it.isend()) {
        // set PDE (only executed if PdeAuxDim>0, i.e. pde coefficients
        // depend on auxilliary variables)
        SVector<double, PdeAuxDim> y;
        for(size_t i=0; i<PdeAuxDim; i++) {
            y[i]=u.coord(GridDim-PdeAuxDim+i,it.idx()[i]);
            pde.set_aux(y);
            recalc=true;    // pde changes, must recalc
        }
        // solves PDE on the specific slice pointed to by it.idx()
        v.reassign_view(u,it.idx());
        FdmSolveStepHelper(v,pde,boundary,t1,t2-t1,theta,recalc,method);
        recalc=false;       // if pde doesn't change, no need to recalc
        it++;
    }
}
// template specialisation for no auxilliary variables
template <std::size_t GridDim>
void FdmSolveStep(GridFunction<GridDim>& u,
                  PDE<GridDim>& pde,
                  const Boundary<GridDim>& boundary,
                  double t1, double t2,
                  double theta=0.5,
                  bool recalc=true,
                  type::method method=type::pred_corr)
{
    FdmSolveStepHelper(u,pde,boundary,t1,t2-t1,theta,recalc,method);

}



// solves the PDE over the entire time
// executes the one time step solver many times
// returns vector of results (u(x,t0), u(x,t1), ..., u(x,tn))
// TODO: retire
template <std::size_t GridDim>
std::vector<double> FDMSolve(
    PDE<GridDim>& pde, const Boundary<GridDim>& boundary,
    GridFunctionBase<GridDim>& u,
    const Grid<1>& grid_time, const SVector<double,GridDim>& x,
    type::method method=type::pred_corr,
    bool out_screen=false, bool out_file=false)
{

    TimingInfo      timing;     // timing of one time step
    std::vector<double>	result(grid_time.size(0));
    //out_file=true
    fdm_print("INIT",u,x,0.0,0,timing,out_screen,out_file);

    // solving the pde step by step
    // ----------------------------
    result[0]=u(x);
    for(size_t k=1; k<grid_time.size(0); k++) {
        FDMSolveStep(pde,boundary,u,grid_time,x,k,method,out_screen,out_file);
        result[k]=u(x);
    }
    return result;
}


} // namespace fdm


#endif // MF_FDM_SOLVE_PDE_H
