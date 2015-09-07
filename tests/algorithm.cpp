#include <fdm/algorithm.hpp>

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE AlgorithmTest
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
// quasi random generator
inline double qrunif()
{
    const double c = sqrt(2.0)-1.0;
    static double x=0.0;
    x += c;
    x -= floor(x);
    return x;
}
inline double qrunif(double a, double b)
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

// wrap the find functions into functors, to simplify benchmark
class Fastfind
{
public:
    inline int operator() (const std::vector<double>& v, double x,
                           int initial_guess=-1) const
    {
        return fdm::fastfind(v,x,initial_guess);
    }
};
class Std_find
{
public:
    inline int operator() (const std::vector<double>& v, double x,
                           int initial_guess=-1) const
    {
        return fdm::std_find(v,x,initial_guess);
    }
};

// benchmarking find functions
template <class F>
double bench(int num_points, int num_ops, bool rand_seek)
{
    // prepare input grid points
    double  a=0.0, b=1.0;
    std::vector<double> v;
    rand_vector_x(v,a,b,num_points);

    // find function
    F myfind;

    // timer
    boost::timer    t;

    // find idx corresponding to value
    int maxiter=num_ops/ceil(log(num_points));
    int idx=0;
    double dx=2.0*(b-a)/num_points;
    double x=a-0.1;
    t.restart();
    for(int i=0; i<maxiter; i++) {
        if(rand_seek) {
            x=qrunif(a-0.1, b+0.1);
            idx=myfind(v,x);
        } else {
            if(x>b)
                x=a-0.1;
            x+=dx*qrunif();
            idx=myfind(v,x,idx);
        }
    }
    return t.elapsed()/maxiter + 0.0*idx;
}


// ---------------------------------------------------------------------
// begin test cases
// ---------------------------------------------------------------------


// simple manual test to check accuracy of find() functions
BOOST_AUTO_TEST_CASE( TestFind )
{
    std::vector<double>  v = {-7.2, -1.1, 0.0, 0.1, 4.3, 144.7};

    BOOST_CHECK_EQUAL( fdm::fastfind(v,  -7.3) , -1 );
    BOOST_CHECK_EQUAL( fdm::fastfind(v,  -7.2) ,  0 );
    BOOST_CHECK_EQUAL( fdm::fastfind(v,   0.0) ,  2 );
    BOOST_CHECK_EQUAL( fdm::fastfind(v,   3.2) ,  3 );
    BOOST_CHECK_EQUAL( fdm::fastfind(v,   4.3) ,  4 );
    BOOST_CHECK_EQUAL( fdm::fastfind(v,   7.3) ,  4 );
    BOOST_CHECK_EQUAL( fdm::fastfind(v, 144.7) ,  5 );
    BOOST_CHECK_EQUAL( fdm::fastfind(v, 341.6) ,  5 );

    BOOST_CHECK_EQUAL( fdm::std_find(v,  -7.3) , -1 );
    BOOST_CHECK_EQUAL( fdm::std_find(v,  -7.2) ,  0 );
    BOOST_CHECK_EQUAL( fdm::std_find(v,   0.0) ,  2 );
    BOOST_CHECK_EQUAL( fdm::std_find(v,   3.2) ,  3 );
    BOOST_CHECK_EQUAL( fdm::std_find(v,   4.3) ,  4 );
    BOOST_CHECK_EQUAL( fdm::std_find(v,   7.3) ,  4 );
    BOOST_CHECK_EQUAL( fdm::std_find(v, 144.7) ,  5 );
    BOOST_CHECK_EQUAL( fdm::std_find(v, 341.6) ,  5 );
}

// comparing fastfind() and std_find()
BOOST_AUTO_TEST_CASE( TestFind2 )
{
    // prepare input grid points
    const int num_points = 10000;
    const int num_tests  = 100000;
    double x;
    int idx1, idx2;
    double  a=-1.0, b=1.0;
    std::vector<double> v;
    rand_vector_x(v,a,b,num_points);

    // random seek
    for(int i=0; i<num_tests; i++) {
        x = rnorm();
        idx1 = fdm::fastfind(v,x);
        idx2 = fdm::std_find(v,x);
        BOOST_CHECK_EQUAL ( idx1, idx2 );
    }
    // random seek, but providing a bad initial guess
    for(int i=0; i<num_tests; i++) {
        x = rnorm();
        idx1 = fdm::fastfind(v,x,idx1);
        idx2 = fdm::std_find(v,x);
        BOOST_CHECK_EQUAL ( idx1, idx2 );
    }
    // linear seek
    double dx = (b-a) / (num_points-1);
    x = a;
    for(int i=0; i<num_tests; i++) {
        idx1 = fdm::fastfind(v,x,idx1);
        idx2 = fdm::std_find(v,x,idx2);
        BOOST_CHECK_EQUAL ( idx1, idx2 );
        x+=dx;
    }
}



// speed benchmark of find()
BOOST_AUTO_TEST_CASE( BenchFind )
{
    std::vector<int>  num = {100, 1000, 10000, 100000};
    const int       num_ops = 1000000;   // order of operations for benchmark

    printf("\n\t\t\trandom\t\tlinear\n");
    for(size_t i=0; i<num.size(); i++) {
        int n=num[i];
        double timerand=bench<Fastfind>(n,num_ops,true);
        double timelinear=bench<Fastfind>(n,num_ops,false);
        printf("fastfind: %9i\t%7.3fns\t%7.3fns\n", n,
               timerand*1e6, timelinear*1e6);
    }
    printf("\n");
    for(size_t i=0; i<num.size(); i++) {
        int n=num[i];
        double timerand=bench<Std_find>(n,num_ops,true);
        double timelinear=bench<Std_find>(n,num_ops,false);
        printf("std_find: %9i\t%7.3fns\t%7.3fns\n", n,
               timerand*1e6, timelinear*1e6);
    }
}
