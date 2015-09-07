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


#ifndef MF_FDM_COMMON_H
#define MF_FDM_COMMON_H

// trap floating point exceptions
#ifdef USE_FENV
#  include <cfenv>
#endif

// to suppress compiler warning issued by g++ -Wunused-parameter:
// this can be useful, e.g. when templates are instantiated for
// special cases (n=1) where certain variables may not be used
#ifdef UNUSED
#  error "Preprocessor macro UNUSED already defined elsewhere."
#elif defined(__GNUC__)
#  define UNUSED(x) x __attribute__((unused))
#else
#  define UNUSED(x) x
#endif

// check C++11 compatibility
#if __cplusplus >= 201103L
#  define MF_FDM_USE_C11
#else
#  warning compiler does not support C++11, certain features will be disabled,\
e.g. use g++ -std=c++11
#endif




namespace fdm
{

namespace
{

// square of a number sqr(x)=x^2
template <typename T>
inline T sqr(T x)
{
    return x*x;
}


// trap floating point exceptions, and generate SIGFPE signal
// (GNU libc only)
inline bool enableexcept()
{
#ifdef USE_FENV
    // programme will terminate if any of the following errors occur
    // so we catch nan, inf, etc
    //feenableexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
    feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
    return true;
#else
    return false;
#endif
}

} // anonymous namespace

} // namespace fdm


#endif /* MF_FDM_COMMON_H */
