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


//#include <fdm/functors.hpp>

namespace fdm
{

inline int sign(double x)
{
    if(x<0)  return -1;
    if(x==0) return 0;
    return 1;
}


// i have not rewritten this routine
// should only be used in experimental setting
/*
bool FunctionDistance::adapt_scale_factor_nonlin(double a, double b,
			int n_points) {

 double	f0,f1,f_deriv;
 double	x0,x1;
 double	x_mid, f_mid;
 int	i=0;
 const int	max_fixpoint=20;
 const int	max_newton=10;
 const double	accuracy_fix=1E-2;
 const double	accuracy=1E-10;
 const int	integr_points_fix=20;
 const int	integr_points=500;

 // we need to satisfy int_a^b 1/g(x) dx = b-a
 // first try
 wScaleFactor=1.0; f0=integrate_simpson_inv(a,b,integr_points);
 wScaleFactor=f0/(b-a); f1=integrate_simpson_inv(a,b,integr_points);

 if(fabs(f1-(b-a))<=1E-10*(b-a)){
   return true;
 }

 x0=0.0; x1=1.0;
 wScaleFactor=x0; f0=integrate_simpson_inv(a,b,integr_points_fix)-(b-a);
 if(f0<0){
    printf("adapt_scale_factor: can not adapt since wrong parameters have been chosen\n");
    printf("                    using uniform mesh\n");
//    abort();
 }
 if(fabs(f0)<accuracy*(b-a)) return true;

 wScaleFactor=x1; f1=integrate_simpson_inv(a,b,integr_points_fix)-(b-a);
 // using fixpoint iteration to obtain a first approximation
 while((sign(f0)==sign(f1)) && (i<max_fixpoint) ){
   x1 *= 4;
   wScaleFactor=x1; f1=integrate_simpson_inv(a,b,integr_points_fix)-(b-a);
   i++;
 }
 while( (fabs((x1-x0)/x1)>accuracy_fix) && (i<max_fixpoint) ){
   x_mid=0.5*(x0+x1);
   wScaleFactor=x_mid; f_mid=integrate_simpson_inv(a,b,integr_points_fix)-(b-a);
   if(sign(f0)!=sign(f_mid)){
     x1=x_mid;
     f1=f_mid;
   } else {
     x0=x_mid;
     f0=f_mid;
   }
   i++;
 }
 if(i==max_fixpoint){
   printf("adapt_scale_factor: warning: fixpoint method failed\n");
 }
 // using secand method to obtain the right value for the parameter wScaleFactor
 for(i=0;i<max_newton;i++){
    f_deriv=(f1-f0)/(x1-x0);
    x0=x1; f0=f1;
    x1=x0-f0/f_deriv;		// newton step
    wScaleFactor=x1; f1=integrate_simpson_inv(a,b,integr_points)-(b-a);
    if(fabs(f1)<accuracy*(b-a)) break;
 }

 if(fabs(f1)>accuracy*(b-a)){
   printf("adapt_scale_factor: warning adoption failed, remaining error is still %f\n",f1/(b-a));
   return false;
 }
 return true;
}
*/



// using fixpoint interation to obtain wParam so that f(1)=b
// this is slower but more robust than a newton iteration
// we only require that f(b) is monotonic in wParam
bool FunctionMapping::adapt_parameter()
{
    double	p0,p1,p2;
    double	f0,f1,f2;

    const int	fixmax=3*18;	// should give accuracy of about 1e-18

    p0=wParam;
    f0=this->operator()(1.0)-wB;
    p2=p0;

    // finding two parameters p0, p2 which have opposite signs
    // assuming f(b) is monotonic in wParam
    if(f0>0.0) {
        for(int i=0; i<fixmax; i++) {
            p2*=0.5;
            wParam=p2;
            f2=this->operator()(1.0)-wB;
            if(f2<0.0) {
                break;
            }
        }
    } else {
        for(int i=0; i<fixmax; i++) {
            p2*=2.0;
            wParam=p2;
            f2=this->operator()(1.0)-wB;
            if(f2>0.0) {
                break;
            }
        }
    }
    if(sign(f0)==sign(f2)) {
        return false;
    }

    // applying fixpoint
    for(int i=0; i<fixmax; i++) {
        if(p2==p0)	{
            break;
        } else	{
            p1=0.5*(p0+p2);
        }
        wParam=p1;
        f1=this->operator()(1.0)-wB;
        if(f1==0)		{
            return true;
        }

        // chose the right interval where a root is expected
        if(sign(f0)!=sign(f1))	{
            p2=p1;
            f2=f1;
        } else			{
            p0=p1;
            f0=f1;
        }
    }
    //printf("param=%f, f(1)=%f (b=%f)\n",wParam,operator()(1.0),wB);
    return true;
}


// numerical solution of an ode x' = g(x)
std::vector<double> solve_ode_euler(const FunctionREAL& g,
                                    double a, double b, int n_points, double x_0)
{
    std::vector<double>	x(n_points);
    double			dt=(b-a)/(n_points-1);

    x[0]=x_0;
    for(int i=1; i<n_points; i++) {
        x[i]=x[i-1]+g(x[i-1])*dt;
    }
    return x;
}

// numerical solution of an ode x' = g(x)
std::vector<double> solve_ode_runge(const FunctionREAL& g,
                                    double a, double b, int n_points, double x_0)
{

    std::vector<double>	x(n_points);
    double			dt=(b-a)/(n_points-1);
    double			k1,k2,k3,k4;

    x[0]=x_0;
    for(int i=0; i<n_points-1; i++) {
        k1=g(x[i])*dt;
        k2=g(x[i]+0.5*k1)*dt;
        k3=g(x[i]+0.5*k2)*dt;
        k4=g(x[i]+k3)*dt;
        x[i+1]=x[i]+k1/6.0+k2/3.0+k3/3.0+k4/6.0;
    }
    return x;
}

} //  namespace fdm
