/*
 * ---------------------------------------------------------------------
 * Copyright (C) 2002, 2009, 2010, 2015 Tino Kluge (ttk448 at gmail.com)
 * Copyright (C) 2009, 2010 Mathfinance AG (info@mathfinance.com)
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


#ifndef MF_PNG_IO_H
#define MF_PNG_IO_H

#include <boost/multi_array.hpp>


namespace png
{

struct rgb {
public:
    unsigned short int	r,g,b;
};

bool write_png(const boost::multi_array<double,2>& value,
               const char* filename);

bool write_png(const boost::multi_array<rgb,2>& pixel,
               const char* filename);



} // namespace

#endif // MF_PNG_IO_H
