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



#ifndef MF_FDM_GRID_H
#define MF_FDM_GRID_H

#include <cstdio>
#include <cassert>
#include <cstring>
#include <vector>
#include <boost/array.hpp>

#include <fdm/common.hpp>
#include <fdm/svector.hpp>
#include <fdm/algorithm.hpp>
#include <fdm/interpolation.hpp>

namespace fdm
{

// class for the representation of non-uniform, but regular grids
// (tensor grids)
// - it contains N 1d grids
// - the resulting N-dim grid is the direct product of all the 1d grids
// - example, let {x_1, ..., x_n}, and {y_1,...,y_m} be 2 1d grids
//   then {(x_i,y_j)}_{i,j=1}^{n,m} are the resulting 2d grid nodes
//
template <std::size_t GridDim>
class Grid
{
protected:
    // coordinates of the grid, stored as a series of 1d grid points, i.e.
    // m_coord[i] contains all grid points in dimension i
    boost::array< std::vector<double>, GridDim >    m_coord;

public:
    // constructors
    Grid() {};
    virtual ~Grid()
    {
        ;
    }
    // assignment
    virtual void operator = (const Grid& grid)
    {
        m_coord=grid.m_coord;
    }
    // dimension of the gird, ie number of 1-d grids
    inline int dim() const
    {
        return GridDim;
    }
    // number of grid points of the dim-th grid
    inline size_t size(int dim) const
    {
        assert( (0<=dim) && (dim<this->dim()) );
        return m_coord[dim].size();
    }
    SVector<int, GridDim> shape() const
    {
        SVector<int, GridDim>	sh;
        for(int dim=0; dim<this->dim(); dim++) {
            sh[dim]=this->size(dim);
        }
        return sh;
    }
    // total number of grid points
    size_t num_grid_points() const
    {
        size_t	prod=1;
        for(int dim=0; dim<this->dim(); dim++) {
            prod *= this->size(dim);
        }
        return prod;
    }
    size_t size() const
    {
        return num_grid_points();
    }

    // sets and returns all coordinates in one dimension
    inline std::vector<double> coord(int dim) const
    {
        assert( (0<=dim) && (dim<this->dim()) );
        return m_coord[dim];
    }
    inline std::vector<double> & coord(int dim)
    {
        assert( (0<=dim) && (dim<this->dim()) );
        return m_coord[dim];
    }
    // sets and returns single point coordinate in one dimension
    // could be removed as equivalent to above and coord(dim)[i]
    inline double coord(int dim, int i) const
    {
        assert( (0<=dim) && (dim<this->dim()) );
        assert( (0<=i) && (i<(int)m_coord[dim].size()) );
        return m_coord[dim][i];
    }
    inline double & coord(int dim, int i)
    {
        assert( (0<=dim) && (dim<this->dim()) );
        assert( (0<=i) && (i<(int)m_coord[dim].size()) );
        return m_coord[dim][i];
    }
    // returns (x_i,y_j, ...) coordinates of grid at index (i,j, ...)
    SVector<double, GridDim> coord(const SVector<int, GridDim>& idx) const
    {
        SVector<double, GridDim>	x;
        for(int dim=0; dim<(int)GridDim; dim++) {
            x[dim]=coord(dim,idx[dim]);
        }
        return x;
    }
    // lower and upper bounds
    inline double lower(int dim) const
    {
        return coord(dim, 0);
    }
    inline double upper(int dim) const
    {
        return coord(dim, this->size(dim)-1);
    }

    // returns the index of the grid point with coordinate less or equal x
    // returns -1, if no such grid point exists, i.e. x is smaller than all
    // grid points
    int find_leq(int dim, double x, int initial_guess=-1) const
    {
        assert( (0<=dim) && (dim<this->dim()) );
        return fastfind(m_coord[dim],x,initial_guess);
    }
    SVector<int, GridDim> find_leq(const SVector<double, GridDim>& x) const
    {
        SVector<int, GridDim>	idx;
        for(int dim=0; dim<this->dim(); dim++) {
            idx[dim]=find_leq(dim,x[dim]);
        }
        return idx;
    }


    // returns the index of the grid point with coordinate closest to x
    int find_nearest(int dim, double x, int initial_guess=-1) const
    {
        assert( (0<=dim) && (dim<this->dim()) );
        int i;
        i = fastfind(m_coord[dim],x,initial_guess);
        // check whether v[i] or v[i+1] is closest, unless i==-1, i==last index
        if(i<0) {
            i=0;
        } else if(i>=(int)size(dim)-1) {
            i=size(dim)-1;
        } else {
            if(fabs(m_coord[dim][i]-x)>fabs(m_coord[dim][i+1]-x)) {
                i++;
            }
        }
        return i;
    }
    SVector<int, GridDim> find_nearest(const SVector<double, GridDim>& x) const
    {
        SVector<int, GridDim>	idx;
        for(int dim=0; dim<this->dim(); dim++) {
            idx[dim]=find_nearest(dim,x[dim]);
        }
        return idx;
    }

    // insert a grid point into an existing grid
    void insert(int dim, double x)
    {
        assert( (0<=dim) && (dim<this->dim()) );
        int i=find_leq(dim,x);
        // only add if grid point doesn't already exist
        if(i<0 || m_coord[dim].at(i) !=x) {
            i++;
            if(i<(int)size(dim)) {
                m_coord[dim].insert(m_coord[dim].begin()+i, x);
            } else {
                m_coord[dim].push_back(x);
            }
        }
    }


    // set a uniform grid, with n points in the interval [a,b]
    void set_uniform(int dim, double a, double b, int n)
    {
        assert( (a<b) && (n>2) );
        double	dx=(b-a)/(n-1);
        m_coord.at(dim).resize(n);
        for(int i=0; i<n-1; i++) {
            m_coord[dim][i] =  a + dx*i;
        }
        m_coord[dim][n-1] = b;    // make sure bounds are exact
    }
    // set a grid based on a "distance ratio" function
    // g(x) is a measure on how distant neighbouring grid points are near x
    // TODO: need to implement a very good ODE solver
    //   void set_by_density(int dim, const F& g, int n);


    // defines grid points in one dimension by a given mapping
    // g:[0,1] --> [a,b] and monotonic increasing
    // the function class needs to define
    //      double operator()(double) const;
    //      double lower() const;
    //      double upper() const;
    template<class Function>
    void set_by_mapping(int dim, const Function& g, int n)
    {
        assert( (0<=dim) && (dim<this->dim()) && (n>2) );
        m_coord[dim].resize(n);
        double	x=0.0;
        double	dx=1.0/(double)(n-1);
        for(int i=0; i<n; i++) {
            x=dx*i;
            m_coord[dim][i]=g(x);
        }
        // due to numerical errors we might not get the exact bounds [a,b]
        // in that case we linearly stretch the grid, as exact bounds are vital
        if( m_coord[dim][0]!=g.lower() || m_coord[dim][n-1]!=g.upper() ) {
            double error_lo=m_coord[dim][0]-g.lower();
            double error_up=m_coord[dim][n-1]-g.upper();
            double stretch=(error_lo-error_up)/(m_coord[dim][n-1]-m_coord[dim][0]);
            double offset=g.lower()-m_coord[dim][0];
            //TODO throw if errors too big
            for(int i=1; i<n-1; i++) {
                m_coord[dim][i] += offset + stretch*(m_coord[dim][i]-m_coord[dim][0]);
            }
            m_coord[dim][0]   = g.lower();
            m_coord[dim][n-1] = g.upper();
        }
        assert(valid_coord(dim)==true);
    }

    // slightly modifies a grid so it is piecewise uniform
    // might need a double check!!, but seems to work
    void change_to_piecewise_uniform(int dim, int n_uniform)
    {
        assert( (0<=dim) && (dim<this->dim()) && (n_uniform>1) );
        int start,count,n_const;
        double	alpha;
        start=0;
        // number of constant slices
        n_const=(int) floor((double) (size(dim)/n_uniform));
        count=(size(dim)-1)%n_uniform;	// initial number of adjacent points
        if(count==0) count=n_uniform;
        for(int k=0; k<n_const; k++) {
            start=k*n_uniform;
            for(int i=1; i<=count ; i++) {
                alpha=(m_coord[dim][start+count]-m_coord[dim][start])/count;
                m_coord[dim][start+i]=alpha+m_coord[dim][start+i-1];
            }
            count=n_uniform;
        }
        assert(valid_coord(dim)==true);
    }

    // check validity of grid points, eg they have to be strictly monotic
    // increasing
    bool valid_coord(int dim) const
    {
        assert( (0<=dim) && (dim<this->dim()) && this->size(dim)>2 );
        double	tmp=m_coord[dim][0];
        for(size_t i=1; i<this->size(dim); i++) {
            if(tmp>=coord(dim,i)) {
                printf("Grid::valid_coord(%i): (%lu, %f)-->(%lu, %f) "
                       "not increasing\n", dim,i-1,tmp,i,coord(dim,i));
                return false;
            }
            tmp=this->coord(dim,i);
        }
        return true;
    }
    // check whether the given multi-index is within the grid bounds
    bool check_valid_idx(const SVector<int,GridDim>& idx) const
    {
        for(int dim=0; dim<this->dim(); dim++) {
            if( (idx[dim]<0) || (idx[dim]>=(int)this->size(dim)) ) {
                return false;
            }
        }
        return true;
    }

    void export_points(const char* filename)
    {
        FILE *fp;
        fp=fopen(filename,"w");
        assert(fp!=NULL);

        for(int dim=0; dim<this->dim(); dim++) {
            fprintf(fp,"# dimension %i\n",dim);
            for(int i=0; i<this->size(dim); i++) {
                fprintf(fp,"%i, %f\n",i,this->coord(dim,i));
            }
            fprintf(fp,"\n");
        }
        fclose(fp);
    }

};


// derived grid class
//
// additionally, this grid contains information on how to approximate
// derivatives of a function over that grid, using 3 grid points
// C2 is short for twice continiously differentiable
// TODO: possibly merge Grid with GridC2; set m_calculated=false as soon
//       as grid coordinates are changed
template <std::size_t GridDim>
class GridC2 : public Grid<GridDim>
{
protected:
    // whether the grid is in a ready/calculated state?
    bool  m_calculated;

    // factors for approximating derivatives in at point i within dimension dim
    // with a=m_deriv[dim][i], we do
    // a[0]*f[i-2] + a[1]*f[i-1] + a[2]*f[i] + a[3]*f[i+1] + a[4]*f[i+2]
    // but we only use 3 points, usually a[0]=a[4]=0, this is a small waste
    // of memory but we can more easily and efficiently implement derivatives
    // at boundaries
    boost::array< std::vector< boost::array<double,5> >, GridDim > m_deriv1;
    boost::array< std::vector< boost::array<double,5> >, GridDim > m_deriv2;

    // this defines which 3 out of the 5 points are used to approx derivatives
    // e.g. if m_deriv_index_start[dim][i]=-1, we approximate derivatives like
    // a[1]*f[i-1] + a[2]*f[i] + a[3]*f[i+1]
    // i.e. we start with index i+m_deriv_index_start,
    // till i+m_deriv_index_start+2
    // hence allowed values are -2,-1,0, only
    boost::array< std::vector<int>, GridDim >	m_deriv_index_start;

public:
    // constructors
    GridC2()
    {
        m_calculated=false;
    }

    // friend own class
    template <size_t N>
    friend class GridC2;
    // copy constructor from higher to lower dimension
    // require dim as additional (redundand) input to avoid
    // unintentinal casting between different dimensions
    // copies the first dimensions:
    //      grid: 0, ..., GridDim-1, ..., N-1
    //      this: 0, ..., GridDim-1
    // TODO: make it general to allow any list of dimensions, not
    // necessarily the first
    // TODO: better even allow views, rather than copying
    template <size_t N>
    GridC2(const GridC2<N>& grid, size_t UNUSED(dim))
    {
        //static_assert(N==dim, "GridC2(): dimensions must match");
        assert(N==dim);
        static_assert(GridDim<=N, "GridC2(): dimensions incompatible");
        // copy first this depends on the order
        m_calculated=grid.m_calculated;
        for(size_t i=0; i<GridDim; i++) {
            this->m_coord[i]=grid.m_coord[i];
            m_deriv1[i]=grid.m_deriv1[i];
            m_deriv2[i]=grid.m_deriv2[i];
            m_deriv_index_start[i]=grid.m_deriv_index_start[i];
        }
    }
    // assignment
    void operator = (const GridC2& grid)
    {
        this->m_coord=grid.m_coord;
        m_calculated=grid.m_calculated;
        m_deriv1=grid.m_deriv1;
        m_deriv2=grid.m_deriv2;
        m_deriv_index_start=grid.m_deriv_index_start;
    }
    // factor for approximating finite differences
    inline int deriv_index_start(int dim, int i) const
    {
        assert(m_calculated==true);
        assert( (0<=dim) && (dim<this->dim()) );
        assert( (i>=0) && (i<(int)this->size(dim)) );
        return m_deriv_index_start[dim][i];
    }
    inline double fac_deriv1(int dim, int i, int k) const
    {
        // this is only 2nd order accurate if we have a 1-d grid (cd-01)
        assert(m_calculated==true);
        assert( (0<=dim) && (dim<this->dim()) );
        assert( (i>=0) && (i<(int)this->size(dim)) && (-2<=k) && (k<=2) );
        return m_deriv1[dim][i][k+2];
    }
    inline double fac_deriv2(int dim, int i, int k) const
    {
        assert(m_calculated==true);
        assert( (0<=dim) && (dim<this->dim()) );
        assert( (i>=0) && (i<(int)this->size(dim)) && (-2<=k) && (k<=2) );
        return m_deriv2[dim][i][k+2];
    }
    double fac_deriv2(int dim1, int dim2, int i, int j, int k, int l) const
    {
        if(dim1==dim2) {
            if(l==0)	{
                return fac_deriv2(dim1,i,k);
            } else		{
                return 0.0;
            }
        } else {
            return fac_deriv1(dim1,i,k)*fac_deriv1(dim2,j,l);
        }
    }
    // for convenience the same function but with a multi index idx as input
    double fac_deriv1(int dim, const SVector<int, GridDim>& idx, int k) const
    {
        assert( (0<=dim) && (dim<this->dim()) );
        return fac_deriv1(dim,idx[dim],k);
    }
    double fac_deriv2(int dim, const SVector<int, GridDim>& idx, int k) const
    {
        assert( (0<=dim) && (dim<this->dim()) );
        return fac_deriv2(dim,idx[dim],k);
    }
    double fac_deriv2(int dim1, int dim2, const SVector<int, GridDim>& idx,
                      int k, int l) const
    {
        assert( (0<=dim1) && (dim1<this->dim()) );
        assert( (0<=dim2) && (dim2<this->dim()) );
        return fac_deriv2(dim1,dim2,idx[dim1],idx[dim2],k,l);
    }


    // calculate all the derivative approximation factors
    // this method needs to be called after all grid points are set
    void calculate()
    {
        // TODO: only calculate if m_calculated is false
        //       at the moment m_calculated does not get invalidated
        //       if new coordinates are set (coord(), set_uniform(), ...)

        double	dx1,dx2;
        double	x1,x2,x3;
        size_t  i;
        for(int dim=0; dim<this->dim(); dim++) {
            m_deriv1[dim].resize(this->size(dim));
            m_deriv2[dim].resize(this->size(dim));
            m_deriv_index_start[dim].resize(this->size(dim));

            // factors for inner points
            for(i=1; i<this->size(dim)-1; i++) {
                x1=this->m_coord[dim][i-1];
                x2=this->m_coord[dim][i];
                x3=this->m_coord[dim][i+1];

                dx1=x2-x1;
                dx2=x3-x2;

                // determine factors for derivative approximation
                // this is the heart of the fdm scheme in order to obtain
                // a high degree of accuracy
                // --------------------------------------------------
                m_deriv_index_start[dim][i]=-1;

                m_deriv1[dim][i][0]= 0.0;
                m_deriv1[dim][i][1]=-dx2/((dx1+dx2)*dx1);
                m_deriv1[dim][i][2]= (dx2-dx1)/(dx1*dx2);
                m_deriv1[dim][i][3]= dx1/((dx1+dx2)*dx2);
                m_deriv1[dim][i][4]= 0.0;

                m_deriv2[dim][i][0]= 0.0;
                m_deriv2[dim][i][1]= 2.0/((dx1+dx2)*dx1);
                m_deriv2[dim][i][2]=-2.0/(dx1*dx2);
                m_deriv2[dim][i][3]= 2.0/((dx1+dx2)*dx2);
                m_deriv2[dim][i][4]= 0.0;
                // --------------------------------------------------
            }
            // factors for boundary points
            i=0;
            dx1=this->m_coord[dim][i+1]-this->m_coord[dim][i];
            dx2=this->m_coord[dim][i+2]-this->m_coord[dim][i+1];
            // --------------------------------------------------
            m_deriv_index_start[dim][i]=0;

            m_deriv1[dim][i][0]=0.0;
            m_deriv1[dim][i][1]=0.0;
            m_deriv1[dim][i][2]=-(2.0*dx1+dx2)/(dx1*(dx1+dx2));
            m_deriv1[dim][i][3]=(dx1+dx2)/(dx1*dx2);
            m_deriv1[dim][i][4]=-dx1/(dx2*(dx1+dx2));

            m_deriv2[dim][i][0]= 0.0;
            m_deriv2[dim][i][1]= 0.0;
            m_deriv2[dim][i][2]= 2.0/(dx1*(dx1+dx2));
            m_deriv2[dim][i][3]=-2.0/(dx1*dx2);
            m_deriv2[dim][i][4]=2.0/(dx2*(dx1+dx2));
            // --------------------------------------------------

            i=this->size(dim)-1;
            dx1=this->m_coord[dim][i-1]-this->m_coord[dim][i-2];
            dx2=this->m_coord[dim][i]  -this->m_coord[dim][i-1];
            // --------------------------------------------------
            m_deriv_index_start[dim][i]=-2;

            m_deriv1[dim][i][0]=dx2/(dx1*(dx1+dx2));	// check again
            m_deriv1[dim][i][1]=-(dx1+dx2)/(dx1*dx2);
            m_deriv1[dim][i][2]=(2.0*dx2+dx1)/(dx2*(dx1+dx2));
            m_deriv1[dim][i][3]=0.0;
            m_deriv1[dim][i][4]=0.0;

            m_deriv2[dim][i][0]= 2.0/(dx1*(dx1+dx2));
            m_deriv2[dim][i][1]=-2.0/(dx1*dx2);
            m_deriv2[dim][i][2]= 2.0/(dx2*(dx1+dx2));
            m_deriv2[dim][i][3]=0.0;
            m_deriv2[dim][i][4]=0.0;
            // --------------------------------------------------
        }
        m_calculated=true;
    }

};



// GridIterator: iterates over multi-indices of a particular grid
// and contains mechanism for multi-index --> 1d index embedding
// can also check if at boundary and if so at which one
template <std::size_t GridDim>
class GridIterator
{
protected:
    SVector<int, GridDim>   m_shape; // number of grid points per dim
    SVector<int, GridDim>   m_idx;   // current multi index
    int	                    m_ind;   // corresponding embedded 1d index
    bool                    m_end;   // have we reached the end
    // boundary identification: if inner point m_bound_dim=-1
    // otherwise it will contain the smallest dim where
    // m_idx[dim]=0 or m_idx[dim]=m_shape[dim]-1
    int m_bound_dim;
    int m_bound_up_down;    // 0 lower, 1 upper boundary
public:
    // constructors
    GridIterator() { };
    void set_shape(const SVector<int, GridDim>& shape)
    {
        for(int i=0; i<(int)GridDim; i++) {
            assert(shape[i]>0);
        }
        m_shape=shape;
    }
    GridIterator(const SVector<int, GridDim>& shape)
    {
        set_shape(shape);
    };

    // dimension of the grid
    int dim() const
    {
        return GridDim;
    }
    // check whether we are at a boundary point, returns dimension
    // at which we are at a boundary, -1 if inner point
    // note: here we are using a simplified boundary model insofar as
    //   we only return one dimension but we could be at a boundary at
    //   many more dimensions (eg vertex=0d point, edge=1d line, face=2d)
    //   and so if we wanted to return the full picture we would need to
    //   return an array of bool telling whether we are at boundary at
    //   each dimension
    int check_boundary()
    {
        m_bound_dim=-1;
        for(int dim=0; dim<this->dim(); dim++) {
            if(m_idx[dim]==0) {
                m_bound_dim=dim;
                m_bound_up_down=0;
                return m_bound_dim;
            } else if(m_idx[dim]==m_shape[dim]-1) {
                m_bound_dim=dim;
                m_bound_up_down=1;
                return m_bound_dim;
            }
        }
        return m_bound_dim;
    }
    // returns the boundary dimension, assuming it has been calculated
    // beforehand
    int boundary () const
    {
        assert(m_bound_dim>=-1);
        return m_bound_dim;
    }
    // 0 = lower boundary, 1 = upper boundary
    int boundary_updown () const
    {
        assert(m_bound_dim>=0);
        return m_bound_up_down;
    }

    // sometimes we might need to know at which other dimensions we are
    // at a boundary, this function tells us this
    // -1 inner point at dim, 0 lower bound, 1 upper bound
    int is_boundary(int dim) const
    {
        assert( (0<=dim) && (dim<(int)GridDim) );
        if(m_idx[dim]==0) {
            return 0;
        } else if(m_idx[dim]==m_shape[dim]-1) {
            return 1;
        } else {
            return -1;
        }
    }

    // set index to zero
    void begin()
    {
        for(int i=0; i<dim(); i++) {
            m_idx[i]=0;
        }
        m_ind=0;
        m_end=false;
        m_bound_dim=0;
        m_bound_up_down=0;
    }
    void begin_inner()
    {
        for(int i=0; i<dim(); i++) {
            m_idx[i]=1;
        }
        m_ind=this->ind(m_idx);
        m_end=false;
        m_bound_dim=-1;
        m_bound_up_down=0;
    }
    bool isend() const
    {
        return m_end;
    }
    // return the current index
    SVector<int, GridDim> idx() const
    {
        return m_idx;
    }
    int ind() const
    {
        return m_ind;
    }
    // check whether the given multi-index is within the shape bounds
    bool check_valid_idx(const SVector<int,GridDim>& idx) const
    {
        for(int dim=0; dim<this->dim(); dim++) {
            if( (idx[dim]<0) || (idx[dim]>=m_shape[dim]) ) {
                return false;
            }
        }
        return true;
    }
    // set current index
    void set_idx(const SVector<int, GridDim>& idx)
    {
        assert(check_valid_idx(idx)==true);
        m_idx=idx;
        check_boundary();
        // TODO: also check if m_end is reached?
    }

    // multi-index --> 1d-index embedding
    // this is how data is stored internally in a vector
    int ind(const SVector<int, GridDim>& idx) const
    {
        return index_embed_little(m_shape,idx);
    }
    // multi-index --> 1d-index embedding
    // the leading dimension ldim can be specified
    // (adjacent points in dimension ldim are adjacent in the
    //  embedded 1d vector)
    // this is useful for ADI methods, but this is
    // not how data is stored internally in the vector
    int index_embed(const SVector<int, GridDim>& idx, int ldim) const
    {
        return index_embed_little(m_shape,idx,ldim);
    }

    // goes through all possible points within the range given by m_shape
    // TODO: change the order, to make it consistent with going in steps
    // of one with the embeded vector
    void operator ++(int)
    {
        bool	success=false;
        int	dim=this->dim()-1;

        // makes the indicator, whether we're at a boundary, invalid
        m_bound_dim=-2;

        // tries to add a "one" to the last dimension until we reach
        // the first dimension (dim=0)
        while((success==false)&&(dim>=0)) {
            if(m_idx[dim]<m_shape[dim]-1) {
                m_idx[dim]++;
                success=true;
            } else {
                m_idx[dim]=0;
            }
            dim--;
        }
        // increase 1d-index by 1
        m_ind++;
        m_ind=this->ind(m_idx);       // TODO
        assert(m_ind == this->ind(m_idx));  // check consistency
        if(success==false) {
            m_end=true;
            m_ind=0;
        }
    }
    // goes through all possible inner points
    void incr_inner()
    {
        bool	success=false;
        int	dim=this->dim()-1;

        // makes the indicator whether we're at a boundary invalid
        m_bound_dim=-2;

        // tries to add a "one" to the last dimension until we reach
        // the first dimension (dim=0)
        while((success==false)&&(dim>=0)) {
            if(m_idx[dim]<m_shape[dim]-2) {
                m_idx[dim]++;
                success=true;
            } else {
                m_idx[dim]=1;
            }
            dim--;
        }
        // increase 1d-index by 1
        if(dim==this->dim()-2) {
            m_ind++;
            m_ind=this->ind(m_idx); //TODO
            assert(m_ind == this->ind(m_idx));  // check consistency
        } else {
            m_ind=this->ind(m_idx);
        }
        if(success==false) m_end=true;
    }

    // goes through all possible boundary points
    void incr_bound()
    {
        int	dim=this->dim()-1;
        bool	success=false;

        // consider example with m_bound_dim=3 (ie position 4)
        // [x,x,x, 0, y,y]
        // then all x's have to be within {1 ... n-2} and
        // all y's within {0 ... n-1}

        // tries to add a "one" to the last dimension until we reach
        // the boundary dimension (m_bound_dim)
        while((success==false)&&(dim>m_bound_dim)) {
            if(m_idx[dim]<m_shape[dim]-1) {
                m_idx[dim]++;
                success=true;
            } else {
                m_idx[dim]=0;
            }
            dim--;
        }
        // tries to add a "one" to the dimensions before the boundary
        // dimension (wDirection) until we reach dimension 0
        // only generate inner points
        if(success==false) {
            dim=m_bound_dim-1;
            while((success==false)&&(dim>=0)) {
                if(m_idx[dim]<m_shape[dim]-2) {
                    m_idx[dim]++;
                    success=true;
                } else {
                    m_idx[dim]=1;
                }
                dim--;
            }
        }

        // changes the boundary
        if(success==false) {
            if(m_idx[m_bound_dim]==0) {
                // changes to upper boundary
                m_idx[m_bound_dim]=m_shape[m_bound_dim]-1;
                m_bound_up_down=1;
            } else {
                // new boundary
                m_idx[m_bound_dim]=1;
                m_bound_dim++;
                if(m_bound_dim>=this->dim()) {
                    // we reached the end --> set m_index to zero
                    m_bound_dim=-1;
                    m_end=true;
                } else {
                    assert(m_idx[m_bound_dim]==0);
                    m_bound_up_down=0;
                }
            }
        }
        m_ind=this->ind(m_idx);
        assert(m_bound_dim==this->check_boundary());
    }


};



// base class for the representation of scalar fields, i.e.
// functions f:R^n->R
template <std::size_t GridDim>
class ScalarField
{
protected:
public:
    // returns the value at the point
    virtual double operator() (const SVector<double,GridDim>& x) const=0;
    virtual int dim() const
    {
        return GridDim;
    }
    virtual ~ScalarField()
    {
        ;
    }
};



// scalar field, f:R^n-->R, given only values on grid points
// using second order interpolation in between;
// abstract base class
// needed to provide same interface to main class as well views/slices
template <std::size_t GridDim>
class GridFunctionBase : public ScalarField<GridDim>
{
protected:
    // we shouldn't own any data, but may have to own m_grid for simplicity
    GridC2<GridDim>	m_grid;     // grid coordinates
    double*             m_value;    // pointer to values at grid points
public:
    // iterator which can go through all grid points
    // also provides n-D --> 1-d vector embedding
    // TODO: redesign to move it private or external
    GridIterator<GridDim>	iterator;

    // constructor
    GridFunctionBase()
    {
        m_value=NULL;
    }
    // copy constructor
    GridFunctionBase(const GridFunctionBase<GridDim>& f) :
        m_grid(f.m_grid), iterator(f.iterator)
    {
        m_value=NULL;
    }
    // friend my own class otherwise won't have access to protected/private
    // members of same class but with different dimension
    template <std::size_t N>
    friend class GridFunctionBase;
    // copy constructor for lower dimensions
    template <std::size_t N>
    GridFunctionBase(const GridFunctionBase<N>& f) :
        m_grid(f.m_grid,N)
    {
        static_assert(N>=GridDim, "dimensions do not match");
        m_value=NULL;
    }


    // copies function values only, no check that grids are compatible
    void assign_values(const GridFunctionBase<GridDim>& f)
    {
        assert(size()==f.size());
        assert(m_value!=f.m_value);
        std::memcpy(m_value, f.m_value, this->size() * sizeof(double));
    }

    // exporting GridC2 functions
    SVector<double, GridDim> coord(const SVector<int, GridDim>& idx) const
    {
        return m_grid.coord(idx);
    }
    double coord(int dim, int i) const
    {
        return m_grid.coord(dim,i);
    }
    SVector<double, GridDim> coord() const
    {
        return m_grid.coord(this->iterator.idx());
    }
    double coord(int dim) const
    {
        assert( (0<=dim) && (dim<this->dim()) );
        return m_grid.coord(dim,this->iterator.idx()[dim]);
    }
    SVector<int, GridDim> idx() const
    {
        return this->iterator.idx();
    }
    SVector<int, GridDim> shape() const
    {
        return this->m_grid.shape();
    }
    virtual size_t size(int dim) const
    {
        return this->m_grid.size(dim);
    }
    virtual size_t size() const
    {
        return this->m_grid.size();
    }

    size_t num_grid_points() const
    {
        return this->m_grid.num_grid_points();
    }
    SVector<int, GridDim> find_nearest(const SVector<double, GridDim>& x) const
    {
        return m_grid.find_nearest(x);
    }
    int deriv_index_start(int dim, int i) const
    {
        return m_grid.deriv_index_start(dim,i);
    }
    double fac_deriv1(int dim, int i, int k) const
    {
        return m_grid.fac_deriv1(dim,i,k);
    }
    double fac_deriv2(int dim, int i, int k) const
    {
        return m_grid.fac_deriv2(dim,i,k);
    }
    double fac_deriv2(int dim1, int dim2, int i, int j, int k, int l) const
    {
        return m_grid.fac_deriv2(dim1,dim2,i,j,k,l);
    }
    double fac_deriv1(int dim, const SVector<int, GridDim>& idx, int k) const
    {
        return m_grid.fac_deriv1(dim,idx,k);
    }
    double fac_deriv2(int dim, const SVector<int, GridDim>& idx, int k) const
    {
        return m_grid.fac_deriv2(dim,idx,k);
    }
    double fac_deriv2(int dim1, int dim2, const SVector<int, GridDim>& idx,
                      int k, int l) const
    {
        return m_grid.fac_deriv2(dim1,dim2,idx,k,l);
    }

    // read/write function values at grid points
    double value(const SVector<int, GridDim>& idx) const
    {
        assert(m_grid.check_valid_idx(idx)==true);
        int	i=iterator.ind(idx);
        assert( (0<=i) && (i<(int)size()) );
        return	m_value[i];
    }
    double & value(const SVector<int, GridDim>& idx)
    {
        assert(m_grid.check_valid_idx(idx)==true);
        int	i=iterator.ind(idx);
        assert( (0<=i) && (i<(int)size()) );
        return	m_value[i];
    }
    // function value at the current grid point, pointed to by the iterator
    inline double value() const
    {
        int	i=iterator.ind();
        assert( (0<=i) && (i<size()) );
        return	m_value[i];
    }
    inline double & value()
    {
        int	i=iterator.ind();
        assert( (0<=i) && (i<(int)size()) );
        return	m_value[i];
    }

    // approximations of derivative values at grid points
    // --------------------------------------------------
    // first derivative, (d f)/(d_dim)
    double f_1(int dim, const SVector<int, GridDim>& idx) const
    {
        assert( (0<=dim) && (dim<this->dim()) );
        assert(m_grid.check_valid_idx(idx)==true);
        int	istart = m_grid.deriv_index_start(dim,idx[dim]);
        double	value=0.0;
        SVector<int, GridDim>	idx2(idx);
        for(int k=istart; k<istart+3; k++) {
            idx2[dim]=idx[dim]+k;
            value += m_grid.fac_deriv1(dim,idx,k)*this->value(idx2);
        }
        return value;
    }
    // second derivative, (d^2 f) / (d_dim1 d_dim2)
    double f_2(int dim1, int dim2, const SVector<int, GridDim>& idx) const
    {
        assert( (0<=dim1) && (dim1<this->dim()) );
        assert( (0<=dim2) && (dim2<this->dim()) );
        assert(m_grid.check_valid_idx(idx)==true);
        int	istart1 = m_grid.deriv_index_start(dim1,idx[dim1]);
        int	istart2 = m_grid.deriv_index_start(dim2,idx[dim2]);
        double	value=0.0;
        SVector<int, GridDim>	idx2(idx);
        if(dim1==dim2) {
            // TODO: g++ -O3 will sometimes give the following:
            // warning: assuming signed overflow does not occur when assuming
            // that (X - c) <= X is always true [-Wstrict-overflow]
            for(int k=istart1; k<istart1+3; k++) {
                idx2[dim1]=idx[dim1]+k;
                value += m_grid.fac_deriv2(dim1,idx,k)*this->value(idx2);
            }
        } else {
            for(int k=istart1; k<istart1+3; k++) {
                idx2[dim1]=idx[dim1]+k;
                for(int l=istart2; l<istart2+3; l++) {
                    idx2[dim2]=idx[dim2]+l;
                    value += m_grid.fac_deriv2(dim1,dim2,idx,k,l)*this->value(idx2);
                }
            }
        }
        return value;
    }

    // interpolation outside of grid points
    // ------------------------------------

    // hermite cubic spline interpolation
    double hermite_spline(const SVector<double, GridDim>& x) const
    {
        SVector<int, GridDim>   idx;
        const double            eps=1e-14;

        // find the nearest grid point
        idx=m_grid.find_leq(x);

        return hermite_spline_helper<GridDim>(x,idx,eps,
                                              Dummy<GridDim>());
    }

    // interpolation via Taylor series
    double taylor(const SVector<double, GridDim>& x) const
    {
        SVector<int, GridDim>		idx;
        SVector<double, GridDim>	dx;
        double				interpol;

        // find the nearest grid point
        idx=m_grid.find_nearest(x);
        dx=x-m_grid.coord(idx);

        // interpolation using constant trial function
        interpol=this->value(idx);
        // interpolation using linear trial function
        for(int dim=0; dim<this->dim(); dim++) {
            interpol += f_1(dim,idx)*dx[dim];
        }
        // interpolation using quadratic trial function
        for(int dim=0; dim<this->dim(); dim++) {
            for(int dim2=0; dim2<this->dim(); dim2++) {
                interpol += 0.5*(f_2(dim,dim2,idx)*dx[dim]*dx[dim2]);
            }
        }
        return interpol;
    }

    // interpolated value
    double operator() (const SVector<double, GridDim>& x) const
    {
        // check if point x is inside the domain (do not allow extrapolation)
        for(size_t i=0; i<GridDim; i++) {
            assert(m_grid.lower(i)<=x[i] && x[i]<=m_grid.upper(i));
        }
        return hermite_spline(x);
        //return taylor(x);
    }
    // interpolated first derivative
    double f_1(int dim, const SVector<double, GridDim>& x) const
    {
        assert( (0<=dim) && (dim<(int)GridDim) );

        SVector<int, GridDim>		idx;
        SVector<double, GridDim>	dx;
        double				interpol;

        // find the nearest grid point
        idx=m_grid.find_nearest(x);
        dx=x-m_grid.coord(idx);

        // interpolation using constant trial function
        interpol=this->f_1(dim,idx);

        // interpolation using linear trial function
        for(int dim2=0; dim2 < this->dim(); dim2++) {
            interpol += f_2(dim,dim2,idx)*dx[dim2];
        }

        return interpol;
    }

    // interpolated second derivative
    double f_2(int dim1, int dim2, const SVector<double, GridDim>& x) const
    {
        assert( (0<=dim1) && (dim1<this->dim()) );
        assert( (0<=dim2) && (dim2<this->dim()) );

        SVector<int, GridDim>		idx;
        double				interpol;

        // find the nearest grid point
        idx=m_grid.find_nearest(x);

        // interpolation using constant trial function
        interpol=this->f_2(dim1,dim2,idx);

        return interpol;
    }

    // exports function values of all grid points
    // modifies the iterator
    void export_points(const char* filename)
    {
        SVector<double, GridDim>	x;
        FILE *fp;
        fp=fopen(filename,"w");
        if(fp==NULL) {
            printf("exportpoints(%s): cannot create file\n",filename);
            abort();
        }

        this->iterator.begin();
        while(this->iterator.isend()==false) {
            x=this->coord();
            for(int dim=0; dim<this->dim(); dim++) {
                fprintf(fp,"%e ",x[dim]);
            }
            fprintf(fp,"%e\n",this->value());
            this->iterator++;
        }
        fclose(fp);
    }

protected:
    // we hide some detail at the bottom

    // this is ugly, but due to the C++ standard limitation on
    // template specialisation within classes, we make use of
    // overrding rather than specialisation, this requires a dummy class
    template <size_t N>
    struct Dummy {};

    // hermite cubic spline interpolation:
    // we use recursion for the n-dim case, which may not be optimal
    // performance wise, but easy to implement
    // this function returns the interpolation for the point
    // (x_0, ..., x_k, y_{k+1}, ..., y_{GridDim-1})
    // where y=coord(idx), k=Level-1
    // i.e. the first "k+1" coordinates are taken from the input x,
    // and the remaining from the grid point, given by n-dim index idx
    template <size_t Level>
    double hermite_spline_helper(const SVector<double, GridDim>& x,
                                 const SVector<int, GridDim>& idx,
                                 double eps,
                                 Dummy<Level> UNUSED(dummy)) const
    {
        static_assert(Level<=GridDim, "hermite_spline_helper: incorrect Level");
        assert(m_grid.check_valid_idx(idx)==true);
        size_t dim=Level-1;           // dimension/direction we go along
        assert(idx[dim]>=0 && idx[dim] < (int)this->size(dim));
        // if x[dim] corresponds to a grid coordinate, no interpolation
        // needed in that dimension
        if( fabs(x[dim]-coord(dim,idx[dim])) <= fabs(eps*x[dim]) ) {
            return hermite_spline_helper(x, idx, eps, Dummy<Level-1>());
        }
        // else do full one-dimensional interpolation
        assert(idx[dim]<(int)this->size(dim)); // otherwise extrapolation
        Hspline s;
        SVector<int, GridDim> idx1 = idx;
        // check for boundary
        int istart = (idx[dim]>0) ? -1 : 0;
        int istop  = (idx[dim]+2 < (int)this->size(dim)) ? 2 : 1;
        s.resize(istop-istart+1);
        for(int i=istart; i<=istop; i++) {
            idx1[dim]=idx[dim]+i;
            s.x(i-istart) = coord(dim, idx1[dim]);
            s.y(i-istart) = hermite_spline_helper(x, idx1, eps,
                                                  Dummy<Level-1>());
        }
        s.calculate();
        return s.hermite(x[dim],0-istart);
    }
#if 0
    // note: the C++ standard does not allow template specialization within
    //       a class (very annoying) but can only be in namespace
    //       scope and so we can't easily terminate the template recursion.
    // the following code would cause
    //        error: explicit specialization in non-namespace scope
    // see also:
    //        http://stackoverflow.com/questions/3052579/ (second answer)
    //        http://stackoverflow.com/questions/2537716/
    //
    template <>
    double hermite_spline_helper<0>(const SVector<double, GridDim>& x,
                                    const SVector<int, GridDim>& idx,
                                    double eps) const
    {
        // Level=0, return the value at the grid point
        return f(idx);
    }
#endif
    // workaround: use overloading rather than template specialisation
    double hermite_spline_helper(const SVector<double, GridDim>& UNUSED(x),
                                 const SVector<int, GridDim>& idx,
                                 double UNUSED(eps),
                                 Dummy<0> UNUSED(dummy) ) const
    {
        // Level=0, return the value at the grid point
        return value(idx);
    }


};

// forward declare GridFunctionView
template <std::size_t N>
class GridFunctionView;

// scalar field, given only values on grid points
// using second order interpolation in between;
// assuming f in C^2(R^n)
template <std::size_t GridDim>
class GridFunction : public GridFunctionBase<GridDim>
{
private:
    std::vector<double>	m_value_storage; // values at grid points (embedded)
public:
    // template inheritance makes base members invisible unless referenced
    // via this-->..., or "using" statements as below
    // http://stackoverflow.com/questions/4643074/
    using GridFunctionBase<GridDim>::m_grid;
    using GridFunctionBase<GridDim>::m_value;
    using GridFunctionBase<GridDim>::iterator;

    size_t size() const
    {
        // TODO: check this is used
        return m_value_storage.size();
    }
    // the above size() function will also hide any function with the name
    // from base, to unhide, call "using"
    // http://stackoverflow.com/questions/2161462/
    using GridFunctionBase<GridDim>::size;


    // constructor
    GridFunction() { };
    GridFunction(const GridC2<GridDim>& grid)
    {
        m_grid=grid;
        m_grid.calculate();
        m_value_storage.resize(m_grid.num_grid_points());
        m_value = &m_value_storage.at(0);
        iterator.set_shape(m_grid.shape());
    }
    GridFunction(const GridFunction<GridDim>& f)
    {
        m_grid=f.m_grid;
        m_value_storage=f.m_value_storage;
        if(m_value_storage.size()>0)
            m_value = &m_value_storage[0];
        iterator=f.iterator;
    }
    GridFunction(const GridFunctionBase<GridDim>& f) :
        GridFunctionBase<GridDim>(f)
    {
        m_value_storage.resize(f.size());
        if(m_value_storage.size()>0) {
            m_value = &m_value_storage[0];
            this->assign_values(f);     // only call after pointer is set
        }
    }

    // need to friend the GridFunctionView constructor as access to
    // m_value_storage is required
    //template <std::size_t N>
    //friend void GridFunctionView<GridDim>::set_view(GridFunction<N>& f,
    //                                     SVector<int,N-GridDim> idx2);
    template <size_t N>
    friend class GridFunctionView;

};

// view/slice of a GridFunction (possibly lower dimensional)
// does not own any data but links to another GridFunction
template <std::size_t GridDim>
class GridFunctionView : public GridFunctionBase<GridDim>
{
public:
    // template inheritance makes base members invisible unless referenced
    // via this-->..., or "using" statements as below
    // http://stackoverflow.com/questions/4643074/
    using GridFunctionBase<GridDim>::m_grid;
    using GridFunctionBase<GridDim>::m_value;
    using GridFunctionBase<GridDim>::iterator;

    size_t size() const
    {
        return m_grid.size();
    }

    // constructor
    GridFunctionView() { };

private:
    // private helper function (so it can be used in two constructors)
    //
    // create a view from an existing GridFunction
    // possibly from higher to lower dimension
    //      f:      0, ..., GridDim-1, ..., N-1
    //      this:   0, ..., GridDim-1
    //      idx2:   GridDim, ..., N-1
    //
    // TODO: this only works for a view of the first GridDim dimensions,
    //      and requires that vector embedding is continuous in the
    //      first dimensions, and also the grid copy constructor
    //      must do the same
    //      (consider better, more general and robust view/slice design?)

    template <std::size_t N>
    void set_view(GridFunction<N>& f, SVector<int,N-GridDim> idx2)
    {
        // requires that m_grid is already properly copied!
        assert(m_grid.dim()>0 && m_grid.size(0)>0);
        static_assert(GridDim<=N, "GridFunctionView() dimensions incompatible");
        assert(f.m_grid.size()==f.m_value_storage.size());
        assert(f.size()>0);

        // set iterator extensions (assumes m_grid is already correctly set)
        iterator.set_shape(m_grid.shape());

        // find embedded index k where
        // idx = (idx1, idx2), with idx1=(0,...,0), and idx2 as input,
        // points to
        SVector<int,N> idx;
        for(size_t i=0; i<N; i++) {
            if(i<GridDim)
                idx[i]=0;
            else
                idx[i]=idx2[i-GridDim];
        }
        assert(f.iterator.check_valid_idx(idx));
        int k = f.iterator.ind(idx);

        // point to this memory location
        // note, this requires access to private member of other class,
        // --> need to friend this function in GridFunction
        m_value = &f.m_value_storage[k];

        // make sure we will never exceed memory bounds
        assert( k+size() <= f.size() );
    }

public:
    // view/slice
    template <std::size_t N>
    GridFunctionView(GridFunction<N>& f, SVector<int,N-GridDim> idx2) :
        GridFunctionBase<GridDim>(f)    // also calls copy constr of m_grid
    {
        set_view(f,idx2);
    }
    // view/slice based on where the iterator points to
    template <std::size_t N>
    GridFunctionView(GridFunction<N>& f) :
        GridFunctionBase<GridDim>(f)
    {
        SVector<int, N-GridDim> idx2;
        for(size_t i=0; i<N-GridDim; i++) {
            idx2[i] = f.iterator.idx()[i+GridDim];
        }
        set_view(f,idx2);
    }
    template <std::size_t N>
    void reassign_view(GridFunction<N>& f, SVector<int,N-GridDim> idx2)
    {
        set_view(f,idx2);
    }

};




class BoundaryType
{
public:
    // ideally do not change these values as some derived classes of
    // Boundary migth rely on the order of the numbers, e.g.
    // BoundaryConst gives priority to conditions with a lower numerical
    // value on lower dimensional boundaries
    static const int	Value		= 0;	// values given at the boundary
    static const int	FirstDeriv	= 1;	// first deriv given at boundary
    static const int	SecondDeriv	= 2;	// second deriv
    static const int	Free		= 3;	// free boundary, ie pde discretised
};

// class Boundary (abstract)
// it should be able to return the type of boundary condition and
// the value given a boundary grid point (defined by multi-index idx)
// however, a multi-index will only make sense to the boundary class
// if it is aware of the shape of the grid (ie numbers of grid points
// in every dimension), so I find it is more appropriate to supply
// a GridIterator as input argument
template <std::size_t GridDim>
class Boundary
{
public:
    int dim() const
    {
        return GridDim;
    }
    // type of derivative
    // 0 = given value, 1 = given first derivative, 2 = second deriv
    // 3 = free condition, ie implementation of the pde at these points
    virtual int type(const GridIterator<GridDim>& g_iter) const = 0;
    // boundary value, the input parameters are grid point,
    // its coordinates x, and time t
    virtual double value(const GridIterator<GridDim>& g_iter,
                         const SVector<double, GridDim>& x, double t) const = 0;

    // this slower but for convenience
    int type(const GridFunction<GridDim>& u,
             const SVector<int, GridDim>& idx) const
    {
        GridIterator<GridDim>	g_iter;
        g_iter.set_shape(u.size());
        g_iter.set_idx(idx);
        return this->type(g_iter);
    }
    double value(const GridFunction<GridDim>& u,
                 const SVector<int, GridDim>& idx, double t) const
    {
        GridIterator<GridDim>		g_iter;
        SVector<double, GridDim>	x;
        g_iter.set_shape(u.size());
        g_iter.set_idx(idx);
        x=u.coord(idx);
        return this->value(g_iter,x,t);
    }
};


// free Boundary (returns type 3 at every boundary)
// ie the pde is implemented at each boundary
// warning: this will cause the solution to become unstable
// if the pde is of a type that there is incoming flow at any of
// the boundaries
template <std::size_t GridDim>
class BoundaryFree : public Boundary<GridDim>
{
public:
    inline int type(const GridIterator<GridDim>& UNUSED(g_iter)) const
    {
        return 3;
    }
    inline double value(const GridIterator<GridDim>& UNUSED(g_iter),
                        const SVector<double, GridDim>& UNUSED(x),
                        double UNUSED(t)) const
    {
        return 0.0;
    }
};


// boundary, constant type and value at each (n-1) dimensional boundary
// in an n-dim space we have n lower and n upper (n-1) dim boundaries
// lower dimensional boundaries get assigned to the boundary whatever
// the GridIterator says
template <std::size_t GridDim>
class BoundaryConst : public Boundary<GridDim>
{
public:
    SVector<int, GridDim>	TypeLower;
    SVector<int, GridDim>	TypeUpper;
    SVector<double, GridDim>	ValueLower;
    SVector<double, GridDim>	ValueUpper;
    int type(const GridIterator<GridDim>& g_iter) const
    {
        int	UNUSED(dim)=-1, UNUSED(updown)=-1; // only used for asserts
        int	condition_type=2048;               // set to large invalid type
        // if we are at a lower dimensional boundary (ie we are at a
        // boundary at more than one dim) we want to prioritise
        // boundary conditions of type 0 (const) and so we
        // can't use dim=g_iter.boundary();
        for(int i=0; i<(int)GridDim; i++) {
            if(g_iter.is_boundary(i)>=0) {
                updown=g_iter.is_boundary(i);
                if( ((updown==0) && (TypeLower[i]<condition_type)) ) {
                    dim=i;
                    condition_type=TypeLower[i];
                } else if( ((updown==1) && (TypeUpper[i]<condition_type)) ) {
                    dim=i;
                    condition_type=TypeUpper[i];
                }
            }
        }
        assert( (0<=dim) && (dim<(int)GridDim) );
        assert( (0<=condition_type) && (condition_type<=3) );
        return condition_type;
    }
    double value(const GridIterator<GridDim>& g_iter,
                 const SVector<double, GridDim>& UNUSED(x) ,
                 double UNUSED(t) ) const
    {
        int	UNUSED(dim)=-1, UNUSED(updown)=-1; // only used for asserts
        int	condition_type=2048;               // set to large invalid type
        double	value=0.0;
        // if we are at a lower dimensional boundary (ie we are at a
        // boundary at more than one dim) we want to prioritise
        // boundary conditions of type 0 (const) and so we
        // can't use dim=g_iter.boundary();
        for(int i=0; i<(int)GridDim; i++) {
            if(g_iter.is_boundary(i)>=0) {
                updown=g_iter.is_boundary(i);
                if( ((updown==0) && (TypeLower[i]<condition_type)) ) {
                    dim=i;
                    condition_type=TypeLower[i];
                    value=ValueLower[i];
                } else if( ((updown==1) && (TypeUpper[i]<condition_type)) ) {
                    dim=i;
                    condition_type=TypeUpper[i];
                    value=ValueUpper[i];
                }
            }
        }
        assert( (0<=dim) && (dim<(int)GridDim) );
        assert( (0<=condition_type) && (condition_type<=3) );
        return value;
    }
};




} // namespace fdm

#endif // MF_FDM_GRID_H
