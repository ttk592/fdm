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


#ifndef MF_FDM_ALGORITHM_HPP
#define MF_FDM_ALGORITHM_HPP

#include <cassert>
#include <vector>
#include <algorithm>

#include <fdm/common.hpp>

namespace fdm
{


// check that vector is sorted strictly increasing
template <typename T>
void assert_strictly_increasing(const std::vector<T> v)
{
    for(size_t i=0; i<v.size()-1; i++) {
        assert(v[i]<v[i+1]);
    }
}


// fast search function through a sorted std::vector
// using bisection, O(log(N))
// returns the first index so that v[index] <= x
// returns -1 if x < v[0]
template <class T>
int fastfind(const std::vector<T>& v, T x, int initial_guess=-1)
{
    assert(v.size()>0);
#ifdef DEBUG
    assert_strictly_increasing(v);
#endif

    const int n = v.size()-1;

    // check if x is at the left or right side of all vector values
    if(x<=v[0]) {
        if(x<v[0])
            return -1;
        if(x==v[0])
            return 0;
    }
    if(x>=v[n]) {
        return n;
    }

    int down, up, middle;

    // check whether an initial guess for the index is given
    if(initial_guess<0) {
        down = 0;
        up   = n;
    } else {
        down = initial_guess;
        up   = initial_guess + 2;
        // initial guess for down may be incorrect, and we need to reduce
        int dist=1;
        while( x < v[down] ) {
            up=down;
            down=std::max(0, down-dist);
            dist*=2;
        }
        // initial guess for up may be incorrect, and we need to increase
        // dist=1;      // no need as at most one while loop will be enterred
        while( v[up] < x ) {
            down=up;
            up=std::min(up+dist,n);
            dist*=2;
        }
    }

    // x is somewhere in the middle, search via bisection
    while(up>down+1) {
        // middle index: 0.5*(down+up) rounded up
        middle = (up+down)/2 + (up+down)%2;
        if(x==v[middle]) {
            return middle;
        } else if(x<v[middle]) {
            up=middle;
        } else down=middle;
    }
    return down;
}

// std implementation of a fast find, given a sorted vector as input
// returns the first index so that v[index] <= x
// returns -1 if x < v[0]
template <typename T>
int std_find(const std::vector<T>& v, T x, int UNUSED(initial_guess)=-1)
{
    typename std::vector<T>::const_iterator it;
    it=std::upper_bound(v.begin(),v.end(),x);
    return int(it-v.begin())-1;
}


} // namespace fdm


#endif /* MF_FDM_ALGORITHM_HPP */
