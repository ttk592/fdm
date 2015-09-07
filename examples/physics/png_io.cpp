/*
 * ---------------------------------------------------------------------
 * Copyright (C) 2009, 2010, 2015 Tino Kluge (ttk448 at gmail.com)
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

#include "png_io.h"

#include <cstdio>
#include <cmath>
#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <libpng12/png.h>
#include <fdm/common.hpp>

namespace png
{



// for error handling we have to define a function like this
static void writepng_error_handler(png_structp UNUSED(png_ptr),
                                   png_const_charp msg )
{
    // i really don't understand this error handling stuff and i don't want to
    // so let's just quit as soon as we hit an error
    // (it can actually jump back to the main programme but too creepy for me)
    fprintf(stderr, "writepng libpng error: %s\n", msg);
    fflush(stderr);
    abort();
}


// writes a 2d array containing rgb information into a png file
bool write_png(const boost::multi_array<png::rgb,2>& pixel,
               const char* filename)
{
    FILE *fp=NULL;
    fp = fopen(filename, "wb");
    if(fp==NULL) {
        printf("write_png: error writing file %s\n", filename);
        return false;
    }

    // we need to copy the data into an array
    png_byte	**Array;
    int	height, width;
    height=pixel.shape()[0];
    width=pixel.shape()[1];

    Array = new png_byte*[height];
    for(int i=0; i<height; i++) {
        Array[i] = new png_byte[3*width];
        for(int j=0; j<width; j++) {
            assert( pixel[i][j].r<=255 );
            assert( pixel[i][j].g<=255 );
            assert( pixel[i][j].b<=255 );
            Array[i][3*j]=(png_byte) pixel[i][j].r;
            Array[i][3*j+1]=(png_byte) pixel[i][j].g;
            Array[i][3*j+2]=(png_byte) pixel[i][j].b;
        }
    }


    // initialise png_struct and png_info
    png_structp	png_ptr;
    png_infop	info_ptr;

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                      NULL,writepng_error_handler,NULL);
    if(!png_ptr) {
        printf("write_png: error initialising png_ptr\n");
        return(false);
    }

    info_ptr = png_create_info_struct(png_ptr);
    if(!info_ptr) {
        printf("write_png: error creating info_ptr\n");
        png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
        return(false);
    }

    // set error handling (very strange!)
    if(setjmp(png_jmpbuf(png_ptr))) {
        printf("write_png: general error\n");
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(fp);
        return(false);
    }

    // making sure fp is opened in binary mode
    png_init_io(png_ptr, fp);

    // set image parameters
    png_set_compression_level(png_ptr, Z_BEST_COMPRESSION);
    png_set_IHDR(png_ptr, info_ptr, width, height,
                 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    // writing the header
    png_write_info(png_ptr, info_ptr);

    // writing the actual data, provided by pointers to rows
    png_write_image(png_ptr, Array);

    // finishing up
    png_write_end(png_ptr, info_ptr);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    for(int i=0; i<height; i++) {
        delete[] Array[i];
    }
    delete[] Array;
    fclose(fp);

    return true;

}


// this my own thinking of how colour perception works
double retina_response(double mu, double sigma, double x)
{
    return exp(-(x-mu)*(x-mu)/(2*sigma*sigma));
}

png::rgb wavelength2rgb(double x)
{
    png::rgb	pixel;
    boost::array<double, 3>	rp;	// eye retina response to a spectrum x
    boost::array<double, 3>	rgb_val;

    // these are my estimates for the colour receptors
    rp[0] = retina_response(440.0,23.0,x);	// S receptor
    rp[1] = retina_response(540.0,40.0,x); // M receptor
    rp[2] = retina_response(570.0,44.0,x); // L receptor

    // based on R=615, G=550, B=465
    // we want to have the same retina response of R+G+B spectra
    // as a single spectrum would have
    // coefficients have been determined by inverting a matrix
    rgb_val[0] =  0.428*rp[0] - 2.152*rp[1] + 2.313*rp[2];
    rgb_val[1] = -0.397*rp[0] + 1.415*rp[1] - 0.411*rp[2];
    rgb_val[2] =  1.805*rp[0] - 0.000*rp[1] + 0.000*rp[2];

    // making sure all rgb values are within [0,1]
    // here we should think of a better way, maybe a projection
    // into the space alpha e1 + beta e2 + gamma e3, with coeff in [0,1]
    double max=0.0;
    for(int i=0; i<3; i++) {
        if(rgb_val[i]<0.0) rgb_val[i]=0.0;
        if(rgb_val[i]>1.0) rgb_val[i]=1.0;
        if(rgb_val[i]>max) max=rgb_val[i];
    }
    // normalising to get maximum brightness
    for(int i=0; i<3; i++) {
        //   rgb_val[i]/=max;
    }

    pixel.r=(int) (rgb_val[0]*255.0);
    pixel.g=(int) (rgb_val[1]*255.0);
    pixel.b=(int) (rgb_val[2]*255.0);

    return pixel;
}


bool write_png(const boost::multi_array<double,2>& value,
               const char* filename)
{

    int	n0=value.shape()[0];
    int	n1=value.shape()[1];
    boost::multi_array<png::rgb,2>	pixels(boost::extents[n0][n1]);

    double	x;
    for(int i=0; i<n0; i++) {
        for(int j=0; j<n1; j++) {
            x=value[i][j];
            if(x<0.0) x=0.0;
            if(x>1.0) x=1.0;
            pixels[i][j]=wavelength2rgb(420.0+200.0*x);
        }
    }
    return write_png(pixels,filename);
}



} // namespace
