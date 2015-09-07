#include <fdm/grid.hpp>

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>

#include <fdm/svector.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GridInterpolationTest
#include <boost/test/unit_test.hpp>
#include <boost/timer.hpp>
#include <boost/random.hpp>


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
// random coordinates
template <size_t GridDim>
void runif(fdm::SVector<double,GridDim>& x,
           const fdm::SVector<double,GridDim>& a,
           const fdm::SVector<double,GridDim>& b)
{
    for(size_t dim=0; dim<GridDim; dim++) {
        x[dim] = runif(a[dim], b[dim]);
    }
}


// exact function we will derive the grid from
template <size_t GridDim>
class Function
{
public:
    double operator() (const fdm::SVector<double, GridDim>& x)  const
    {
        return sin(x.norm());
    }
};

// first derivative: of any function object F
// df/d_dim
template <typename F, size_t GridDim>
double deriv1(const F& f, int dim, const fdm::SVector<double, GridDim>& x,
              double h=1e-8)
{
    assert(0<=dim && dim<(int)GridDim);
    fdm::SVector<double, GridDim> x1(x), x2(x);
    x1[dim]=x[dim]-h;
    x2[dim]=x[dim]+h;
    return (f(x2)-f(x1))/(2.0*h);
}
// second derivative: of any function object F
// d^2 f / (d_dim1 d_dim2)
template <typename F, size_t GridDim>
double deriv2(const F& f, int dim1, int dim2,
              const fdm::SVector<double, GridDim>& x,
              double h=1e-6)
{
    assert(0<=dim1 && dim1<(int)GridDim);
    assert(0<=dim2 && dim2<(int)GridDim);
    if(dim1==dim2) {
        fdm::SVector<double, GridDim> x1(x), x2(x);
        x1[dim1]-=h;
        x2[dim1]+=h;
        return (f(x1) - 2.0*f(x) + f(x2))/(h*h);
    } else {
        fdm::SVector<double, GridDim> x11(x), x12(x), x21(x), x22(x);
        x11[dim1]-=h;
        x11[dim2]-=h;
        x12[dim1]-=h;
        x12[dim2]+=h;
        x21[dim1]+=h;
        x21[dim2]-=h;
        x22[dim1]+=h;
        x21[dim2]+=h;
        return (f(x22)-f(x12)-f(x21)+f(x11))/(4.0*h*h);
    }
}


// initialises the grid with boundaries [a_i,b_i], and n_i numbers of
// grid points, then sets all grid points to the provided function
// then checks accuracy of interpolation
template <typename F, size_t GridDim>
double test_grid( const fdm::SVector<double,GridDim>& a,
                  const fdm::SVector<double,GridDim>& b,
                  const fdm::SVector<int,GridDim>& n,
                  const F& f, int num_trials)
{
    // set grid coordinates
    fdm::GridC2<GridDim> grid;
    for(size_t i=0; i<GridDim; i++) {
        grid.set_uniform(i, a[i], b[i], n[i]);
    }
    // set grid values
    boost::timer t;
    t.restart();
    fdm::GridFunction<GridDim> u(grid);
    fdm::SVector<double, GridDim> x;
    u.iterator.begin();
    while(u.iterator.isend()==false) {
        x=u.coord();
        u.value()=f(x);
        u.iterator++;
    }
    printf(" %lu-dim: set values (%lu nodes)\t%.3fs\n",
           GridDim, grid.num_grid_points(), t.elapsed());

    // check interpolation at random points
    t.restart();
    double max=0.0;
    double avg=0.0;
    for(int i=0; i<num_trials; i++) {
        runif(x,a,b);   // x = random coordinate in [a,b]
        double diff = fabs(u(x)-f(x));
        avg += diff;
        max = std::max(max, diff);
    }
    avg /= (double) num_trials;
    printf(" %lu-dim: random access (%i trials)\t%.3fs\t(avg=%f, max=%f)\n",
           GridDim, num_trials, t.elapsed(), avg, max);

    // check first derivatives at random points
    t.restart();
    double h=1e-4;
    max=0.0;
    avg=0.0;
    int num_trials2=num_trials/GridDim/2;
    for(int i=0; i<num_trials2; i++) {
        runif(x,a,b);   // x = random coordinate in [a,b]
        double diff=0.0;
        for(size_t dim=0; dim<GridDim; dim++) {
            double h1 = std::min( std::min(h, x[dim]-a[dim]), b[dim]-x[dim] );
            assert(h1>0);
            double diff_i = deriv1(u,dim,x,h1) - deriv1(f,dim,x,h1);
            diff += diff_i*diff_i;
        }
        diff = sqrt(diff);
        avg += diff;
        max = std::max(max, diff);
    }
    avg /= (double) num_trials2;
    printf(" %lu-dim: first deriv (%i trials)\t%.3fs\t(avg=%f, max=%f)\n",
           GridDim, num_trials2, t.elapsed(), avg, max);

    // check second derivatives at random points
    t.restart();
    h=1e-4;
    max=0.0;
    avg=0.0;
    num_trials2=num_trials/GridDim/GridDim/4;
    for(int i=0; i<num_trials2; i++) {
        runif(x,a,b);   // x = random coordinate in [a,b]
        double diff=0.0;
        for(size_t dim1=0; dim1<GridDim; dim1++) {
            for(size_t dim2=0; dim2<GridDim; dim2++) {
                double h1,h2;
                h1 = std::min(std::min(h,x[dim1]-a[dim1]),b[dim1]-x[dim1]);
                h2 = std::min(std::min(h,x[dim2]-a[dim2]),b[dim2]-x[dim2]);
                h1 = std::min(h1,h2);
                assert(h1>0);
                double diff_i = deriv2(u,dim1,dim2,x,h1)
                                - deriv2(f,dim1,dim2,x,h1);
                diff += diff_i*diff_i;
            }
        }
        diff = sqrt(diff);
        avg += diff;
        max = std::max(max, diff);
    }
    avg /= (double) num_trials2;
    printf(" %lu-dim: second deriv (%i trials)\t%.3fs\t(avg=%f, max=%f)\n",
           GridDim, num_trials2, t.elapsed(), avg, max);



    return max;
}



BOOST_AUTO_TEST_CASE( GridInterpol1D )
{
    fdm::SVector<double, 1> a = {-1.3};
    fdm::SVector<double, 1> b = {2.2};
    fdm::SVector<int, 1>    n = {120};
    int num_trials=500000;
    Function<1> f;
    test_grid(a,b,n,f,num_trials);

}

BOOST_AUTO_TEST_CASE( GridInterpol3D )
{
    fdm::SVector<double, 3> a = {-1.3, -2.0, 0.0};
    fdm::SVector<double, 3> b = {2.2, 4.1, 1.2};
    fdm::SVector<int, 3>    n = {120, 90, 30};
    int num_trials=500000;
    Function<3> f;
    test_grid(a,b,n,f,num_trials);

}


// comparing 1d grid interpolation against the interpolator
BOOST_AUTO_TEST_CASE( GridVsInterpol )
{

}
