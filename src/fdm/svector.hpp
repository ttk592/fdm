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


#ifndef MF_FDM_SVECTOR_H
#define MF_FDM_SVECTOR_H

#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include <boost/array.hpp>

#include <fdm/common.hpp>



namespace fdm
{

// class SVector: static / simple / short vector
// also used in the code to indicate the vector is of low dimensionality
// e.g. 1,2, or 3 dim and represents a point in space
// as opposed to a long vector (100 or 1000's dim) to solve a lin eq system
template <class Type, std::size_t VecDim>
class SVector
{
private:
    Type m_X[VecDim];
public:
    // constructor (default)
    SVector () { };
#ifdef MF_FDM_USE_C11
    // constructor with initializer list
    // so we can instantiate like SVector<int,2> x = {3,7};
    // requires C++11
    SVector (const std::initializer_list<Type>& x)
    {
        assert(x.size()==VecDim);
        size_t i=0;
        for (const Type* it = x.begin(); it != x.end(); it++) {
            this->operator[](i) = *it;
            i++;
        }
    }
#endif
    // copy constructor
    SVector (const SVector& x)
    {
        for(std::size_t i=0; i<VecDim; i++) m_X[i]=x.m_X[i];
    }
    // read the dimension
    inline std::size_t size() const
    {
        return VecDim;
    }
    // access operator
    inline Type operator [] (std::size_t dim) const
    {
        assert(dim<VecDim);   // no need to check for 0<=dim, as size_t >= 0
        return m_X[dim];
    }
    // access operator (write access)
    inline Type& operator [] (std::size_t dim)
    {
        assert(dim<VecDim);
        return m_X[dim];
    }
    // access operator with bounds check
    Type at(std::size_t dim) const
    {
        if(dim>=VecDim) {
            printf("SVector<%lu>.at(%lu) Dimenstion invalid\n",VecDim,dim);
            abort();
        }
        return operator[] (dim);
    }
    // access operator (write) with bounds check
    Type& at(std::size_t dim)
    {
        if(dim>=VecDim) {
            printf("SVector<%lu>.at(%lu) Dimenstion invalid\n",VecDim,dim);
            abort();
        }
        return operator[] (dim);
    }

    // asignment operator
    SVector& operator = (const SVector& x)
    {
#ifndef NDEBUG
        if(this==&x) {
            printf("SVector = : self assignment not allowed\n");
            abort();
        }
#endif
        for(std::size_t i=0; i<VecDim; i++) m_X[i]=x.m_X[i];
        return *this;
    }

    // export to bool::array
    boost::array<Type, VecDim> boost_array()
    {
        boost::array<Type, VecDim>	x;
        for(std::size_t i=0; i<VecDim; i++) {
            x[i]=m_X[i];
        }
        return x;
    }


    // check of equality
    bool operator == (const SVector& x)
    {
        for(std::size_t i=0; i<VecDim; i++) if(m_X[i]!=x.m_X[i]) return false;
        return true;
    }
    // addition operator
    friend SVector operator + (const SVector& x, const SVector& y)
    {
        SVector z;
        for(std::size_t i=0; i<VecDim; i++) z[i]=x[i]+y[i];
        return z;
    }
    // inverse operator
    friend SVector operator - (const SVector& x)
    {
        SVector z;
        for(std::size_t i=0; i<VecDim; i++) z[i]=-x[i];
        return z;
    }
    // substraction operator
    friend SVector operator - (const SVector& x, const SVector& y)
    {
        SVector z;
        for(std::size_t i=0; i<VecDim; i++) z[i]=x[i]-y[i];
        return z;
    }

    // norm of the vector
    double norm() const
    {
        double	val=0.0;
        for(std::size_t i=0; i<VecDim; i++) {
            val += this->operator[](i) * this->operator[](i);
        }
        return sqrt(val);
    }

};


// if the vector degenerates to a one dimensional object
// it should be as efficient as a 1d Type and be interchangable
template <class Type>
class SVector<Type,1>
{
private:
    Type	m_X;
public:
    // constructor (default)
    SVector () { };
    // copy constructor
    SVector (const SVector& x)
    {
        m_X=x.m_X;
    }
    // constructor so that simple types (double/int) can cast to SVector
    SVector (Type x)
    {
        m_X=x;
    }
    // read the dimension
    inline std::size_t size() const
    {
        return 1;
    }
    // access operator
    inline Type operator [] (std::size_t UNUSED(dim)) const
    {
        assert((dim==0));
        return m_X;
    }
    // access operator (write access)
    inline Type& operator [] (std::size_t UNUSED(dim))
    {
        assert((dim==0));
        return m_X;
    }
    // access operator with bounds check
    Type at(std::size_t dim) const
    {
        if((dim<0)||(dim>=size())) {
            printf("SVector<%i>.at(%i) Dimenstion invalid\n",size(),dim);
            abort();
        }
        return operator[] (dim);
    }
    // access operator (write) with bounds check
    Type& at(std::size_t dim)
    {
        if((dim<0)||(dim>=size())) {
            printf("SVector<%i>.at(%i) Dimenstion invalid\n",size(),dim);
            abort();
        }
        return operator[] (dim);
    }

    // asignment operator
    SVector& operator = (const SVector& x)
    {
#ifndef NDEBUG
        if(this==&x) {
            printf("SVector = : self assignment not allowed\n");
            abort();
        }
#endif
        m_X=x.m_X;
        return *this;
    }
    // cast from Type
    SVector& operator = (const Type& x)
    {
        m_X=x;
        return *this;
    }

    // export to bool::array
    boost::array<Type,1> boost_array()
    {
        boost::array<Type,1>	x;
        x[0]=m_X;
        return x;
    }

    // check of equality
    bool operator == (const SVector& x)
    {
        return m_X==x.m_X;
    }
    // addition operator
    friend SVector operator + (const SVector& x, const SVector& y)
    {
        SVector z;
        z[0]=x[0]+y[0];
        return z;
    }
    // inverse operator
    friend SVector operator - (const SVector& x)
    {
        SVector z;
        z[0]=-x[0];
        return z;
    }
    // substraction operator
    friend SVector operator - (const SVector& x, const SVector& y)
    {
        SVector z;
        z[0]=x[0]-y[0];
        return z;
    }
    // norm of the vector
    double norm() const
    {
        return fabs((double)m_X);
    }


};





// index transformation for n-d array --> 1d vector embedding
// "ldim" or leading dim specifies the dimension at which adjacent points
// should be adjacent in the embedded 1d vector
// big endian embedding, (a,b,c) --> (a*n1 + b) * n2 + c
template <std::size_t VecDim>
int index_embed_big(const SVector<int, VecDim>& shape,
                    const SVector<int, VecDim>& idx, int ldim=VecDim-1)
{
    // checks
    static_assert(VecDim>=1, "index_embed_big(): vector dim must be >=1");
    assert((ldim>=0) && (ldim<(int)VecDim));
    for(int dim=0; dim<(int)VecDim; dim++) {
        assert( (idx[dim]>=0) && (idx[dim]<shape[dim]) );
    }

    // n-d index --> 1d index embedding
    if(VecDim==1) {
        // nothing to do in the 1d case
        return idx[0];
    }
    int sum;
    if(ldim==VecDim-1) {
        // default case
        sum=idx[0];
        for(int dim=1; dim<(int)VecDim; dim++)
            sum = sum * shape[dim] + idx[dim];
    } else if(ldim==0) {
        // leading dim is oposite of default
        sum=idx[1];
        for(int dim=2; dim<(int)VecDim; dim++)
            sum = sum * shape[dim] + idx[dim];
        sum = sum * shape[ldim] + idx[ldim];
    } else {
        // leading dim is somewhere in the middle
        sum=idx[0];
        for(int dim=1; dim<ldim; dim++)
            sum = sum * shape[dim] + idx[dim];
        for(int dim=ldim+1; dim<(int)VecDim; dim++)
            sum = sum * shape[dim] + idx[dim];
        sum = sum * shape[ldim] + idx[ldim];
    }
    return sum;
}

// index transformation for n-d array --> 1d vector embedding
// "ldim" or leading dim specifies the dimension at which adjacent points
// should be adjacent in the embedded 1d vector
// little endian embedding, (a,b,c) --> a + (b*n2 + c) * n3
template <std::size_t VecDim>
int index_embed_little(const SVector<int, VecDim>& shape,
                       const SVector<int, VecDim>& idx, int ldim=0)
{
    // checks
    static_assert(VecDim>=1, "index_embed_little(): vector dim must be >=1");
    assert((ldim>=0) && (ldim<(int)VecDim));
    for(int dim=0; dim<(int)VecDim; dim++) {
        assert( (idx[dim]>=0) && (idx[dim]<shape[dim]) );
    }

    // n-d index --> 1d index embedding
    if(VecDim==1) {
        // nothing to do in the 1d case
        return idx[0];
    }
    int sum;
    if(ldim==VecDim-1) {
        // leading dim is oposite of default
        sum=idx[VecDim-2];
        for(int dim=(int)VecDim-3; dim>=0; dim--)
            sum = sum * shape[dim] + idx[dim];
        sum = sum*shape[VecDim-1] + idx[VecDim-1];
    } else if(ldim==0) {
        // default case
        sum=idx[VecDim-1];
        for(int dim=(int)VecDim-2; dim>=0; dim--)
            sum = sum * shape[dim] + idx[dim];
    } else {
        // leading dim is somewhere in the middle
        sum=idx[VecDim-1];
        for(int dim=(int)VecDim-2; dim>ldim; dim--)
            sum = sum * shape[dim] + idx[dim];
        for(int dim=ldim-1; dim>=0; dim--)
            sum = sum * shape[dim] + idx[dim];
        sum = sum * shape[ldim] + idx[ldim];
    }
    return sum;
}


} // namespace fdm

#endif // MF_FDM_SVECTOR_H
