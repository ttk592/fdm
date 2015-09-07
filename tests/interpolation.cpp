#include <fdm/interpolation.hpp>

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE InterpolationTest
#include <boost/test/unit_test.hpp>
#include <boost/timer.hpp>
#include <boost/random.hpp>

// the test is based on a well defined function, which is then
// discretised by grid points, and then interpolated at the grid points
//
// (a) define function f
// (b) discretise at grid points, y_i = f(x_i)
// (c) interpolate grid points (x_i,y_i) by a spline
// (d) compare spline(x) with f(x), as well as derivatives


// function f()
class Function
{
public:
    double operator() (double x)  const
    {
        return sin(x)/(1.0+x*x);
    }
};

// first derivative: of any function object F
template <typename F>
double deriv1(const F& f, double x, double h=1e-8)
{
    double x1=x-h;
    double x2=x+h;
    return (f(x2)-f(x1))/(2.0*h);
}

// second derivative: of any function object F
template <typename F>
double deriv2(const F& f, double x, double h=1e-6)
{
    double x1=x-h;
    double x2=x+h;
    return (f(x1) - 2.0*f(x) + f(x2))/(h*h);
}

// difference between two functions ||f1-f2||
// returns the L2-norm, max-norm, and avg of abs differences
// for function value, 1st and 2nd derivative
template <typename F1, typename F2>
void distance(boost::array<double,3>& l2, boost::array<double,3>& max,
              boost::array<double,3>& avg, const F1& f1, const F2& f2,
              double a, double b, int n)
{
    assert(n>1);
    double dx=(b-a)/(n-1);
    for(int j=0; j<3; j++) {
        l2[j]=0.0;
        max[j]=0.0;
        avg[j]=0.0;
    }
    boost::array<double,3> y;
    for(int i=0; i<n; i++) {
        double x=a+dx*i;
        y[0]=fabs(f1(x)-f2(x));
        y[1]=fabs(deriv1(f1,x)-deriv1(f2,x));
        y[2]=fabs(deriv2(f1,x)-deriv2(f2,x));
        for(int j=0; j<3; j++) {
            l2[j]+=y[j]*y[j];
            avg[j]+=y[j];
            if(max[j]<y[j])
                max[j]=y[j];
        }
    }
    for(int j=0; j<3; j++) {
        avg[j]=avg[j]/(double)n;
        l2[j]=sqrt(l2[j]);
    }
}


// plot grid points for visual inspection
void plot(const std::string& filename,
          const std::vector<double> x, std::vector<double> y)
{
    assert(x.size()==y.size());
    FILE* fp=fopen(filename.c_str(), "w");
    assert(fp!=NULL);
    for(size_t i=0; i<x.size(); i++) {
        fprintf(fp, "%f %f\n",x[i],y[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
}
// plot functions/interpolations for visual inspection
// output into a text file which can be read by gnuplot
template <typename F>
void plot(const std::string& filename, const F& f, double a, double b, int n)
{
    assert(n>1);
    FILE* fp=fopen(filename.c_str(), "a");
    assert(fp!=NULL);
    double dx=(b-a)/(n-1);
    for(int i=0; i<n; i++) {
        double x=a+dx*i;
        double y=f(x);
        double y1=deriv1(f,x);
        double y2=deriv2(f,x);
        fprintf(fp, "%f %f %f %f\n",x,y,y1,y2);
    }
    fprintf(fp,"\n");
    fclose(fp);
}



// random number generator
double runif()
{
    static boost::mt19937   rng;
    static boost::uniform_01<boost::mt19937&> myrand(rng);
    return myrand();
}
double runif(double a, double b)
{
    assert(a<=b);
    return a + (b-a)*runif();
}
double rnorm()
{
    static boost::mt19937   rng;
    static boost::variate_generator<boost::mt19937&,
           boost::normal_distribution<> >
           myrand(rng, boost::normal_distribution<>());
    return myrand();
}
double rnorm(double mu, double sigma)
{
    return mu + sigma*rnorm();
}

// random vector for x-coordinates
void rand_vector_x(std::vector<double>& x, double a, double b, size_t n)
{
    assert(n>1);
    assert(b>a);
    x.resize(n);
    x[0]=0.0;
    for(size_t i=1; i<n; i++) {
        x[i]=x[i-1]+runif();
    }
    // rescale all points so that x[0]=a and x[n-1]=b
    double scale=(b-a)/x[n-1];
    for(size_t i=0; i<n; i++) {
        x[i]=a+scale*x[i];
    }
}
// random vector for y-coordinates
void rand_vector_y(std::vector<double>& x, size_t n)
{
    assert(n>1);
    x.resize(n);
    for(size_t i=0; i<n; i++) {
        x[i]=rnorm();
    }
}


// benchmarking spline routines
template <class S>
boost::array<double,3> bench(int num_points, int num_ops)
{
    // timing information
    boost::array<double, 3> timing;
    // create empty spline
    S spline;
    // prepare input grid points
    double  a=0.0, b=1.0;
    std::vector<double> x, y, x2, y2; // (x[i], y[i]) --> (x2[i], y2[i])
    rand_vector_x(x,a,b,num_points);
    rand_vector_x(x2,a+0.05,b+0.05,num_points);
    rand_vector_y(y,num_points);
    y2.resize(num_points);
    // timer
    boost::timer    t;

    // (1) set points, including solving of eq-systems
    int maxiter=num_ops/num_points;
    t.restart();
    for(int i=0; i<maxiter; i++)
        spline.set_points(x,y);
    timing[0]=t.elapsed()/maxiter;

    // (2) random access function evaluation
    maxiter=2*num_ops/(int)ceil(log(num_points));
    double  sum=0.0;            // to avoid compiler optimising loop away
    double  c=sqrt(2.0)-1.0;    // for quasi-random sequence
    double  z=0.0;
    t.restart();
    for(int i=0; i<maxiter; i++) {
        z+=c;
        z-=(int) z;         // fast quasi random sequence
        sum+=spline(z);
    }
    timing[1]=t.elapsed()/maxiter+0.0*sum;

    // (3) interpolation from one grid to another
    maxiter=num_ops/num_points;
    t.restart();
    for(int i=0; i<maxiter; i++) {
        for(int k=0; k<num_points; k++) {
            y2[k]=spline(x2[k]);
        }
    }
    timing[2]=t.elapsed()/maxiter;

    return timing;
}



// ---------------------------------------------------------------------
// begin test cases
// ---------------------------------------------------------------------


// for plotting a spline over a random grid
BOOST_AUTO_TEST_CASE( InterpolRandom )
{
    const size_t    n = 6;
    double          a = -3.0;
    double          b = 4.5;
    std::vector<double> x,y;

    rand_vector_x(x,a,b,n);
    rand_vector_y(y,n);

    fdm::Hspline hs;
    hs.set_points(x,y);

    fdm::Cspline cs;
    cs.set_points(x,y);

    const int num_plot=500;
    plot("out/interpol_rand.csv",x,y);
    plot("out/interpol_rand.csv",hs,a-1.0, b+1.0,num_plot);
    plot("out/interpol_rand.csv",cs,a-1.0, b+1.0,num_plot);

    //TODO: check md5sum of output

    // check that spline(x[i]) = y[i]
    for(size_t i=0; i<n; i++) {
        BOOST_CHECK_EQUAL( hs(x[i]), y[i] );
        BOOST_CHECK_EQUAL( cs(x[i]), y[i] );
    }
}



// for plotting a spline over a discretised function
BOOST_AUTO_TEST_CASE( InterpolFunction )
{
    const size_t    n = 20;
    double          a = -3.0;
    double          b = 4.5;
    std::vector<double> x,y;

    Function f;
    rand_vector_x(x,a,b,n);
    y.resize(n);
    for(size_t i=0; i<n; i++) {
        y[i]=f(x[i]);
    }

    fdm::Hspline hs;
    hs.set_points(x,y);

    fdm::Cspline cs;
    cs.set_points(x,y);

    const int num_plot=500;
    plot("out/interpol_func.csv",x,y);
    plot("out/interpol_func.csv",hs,a-1.0, b+1.0,num_plot);
    plot("out/interpol_func.csv",cs,a-1.0, b+1.0,num_plot);
    plot("out/interpol_func.csv",f,a-1.0, b+1.0,num_plot);

    //TODO: check md5sum of output
}



// check convergence properties for approximating a function by a spline
// expected 4th-order convergence in function values
BOOST_AUTO_TEST_CASE( InterpolFunctionConvergence )
{
    const int       nmax  = 10000;
    const int       ncompare = 5*nmax;
    double          a = -6.0;
    double          b = 7.5;
    std::vector<double> x,y;
    boost::array<double,3> avg,l2,max;

    Function f;

    int n=5;
    while(n<nmax) {
        // uniform x-grid
        x.resize(n);
        y.resize(n);
        double dx=(b-a)/(n-1);
        for(int i=0; i<n; i++) {
            x[i]=a+dx*i;
            y[i]=f(x[i]);
        }
        // set splines
        fdm::Hspline hs;
        fdm::Cspline cs;
        hs.set_points(x,y);
        cs.set_points(x,y);

        // check distance between spline interpolation and original function
        distance(l2,max,avg,f,hs,a,b,ncompare);
        printf("%i\n",n);
        printf(" hs\tavg:\t%f\t%f\t%f\n", avg[0],avg[1],avg[2]);
        printf("\tmax:\t%f\t%f\t%f\n", max[0],max[1],max[2]);

        distance(l2,max,avg,f,cs,a,b,ncompare);
        printf(" cs\tavg:\t%f\t%f\t%f\n", avg[0],avg[1],avg[2]);
        printf("\tmax:\t%f\t%f\t%f\n", max[0],max[1],max[2]);
        n*=2;
    }
    // TODO: check convergence properties

}


// speed benchmark of interpolation
BOOST_AUTO_TEST_CASE( InterpolBenchmark )
{
    std::vector<int>  num = {100, 1000, 10000, 100000};
    const int       num_ops = 1000000;   // order of operations for benchmark
    boost::array<double,3> timing;

    printf("\n\t\tcreate\t\tf(x)\t\tgrid interpol\n");
    for(size_t i=0; i<num.size(); i++) {
        int n=num[i];
        timing=bench<fdm::Cspline>(n,num_ops);
        printf("Cspline: %i\t%7.3fms\t%7.3fns\t%7.3fms\n", n,
               timing[0]*1000.0, timing[1]*1000000.0, timing[2]*1000.0);
    }
    printf("\n");
    for(size_t i=0; i<num.size(); i++) {
        int n=num[i];
        timing=bench<fdm::Hspline>(n,num_ops);
        printf("Hspline: %i\t%7.3fms\t%7.3fns\t%7.3fms\n", n,
               timing[0]*1000.0, timing[1]*1000000.0, timing[2]*1000.0);
    }


}
