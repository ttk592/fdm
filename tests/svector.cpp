#include <cstdio>
#include <cstdlib>

#include <fdm/svector.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SVectorTest
#include <boost/test/unit_test.hpp>


// check 0d case
// e.g. SVector<double, 0> should work, not at all useful
// but in some edge cases we just want this to exist
BOOST_AUTO_TEST_CASE( Specialisation0D )
{
    fdm::SVector<double,0> x,y;
    BOOST_CHECK ( x==y );
    BOOST_CHECK ( x+y==y-x );
}



// check 1d specialisation functions
// e.g. SVector<double, 1> should be usable similar to a double
BOOST_AUTO_TEST_CASE( Specialisation1D )
{
    fdm::SVector<double,1> x = 1.3;
    BOOST_CHECK_EQUAL ( x[0], 1.3 );

    x = x + 2.5;
    BOOST_CHECK_EQUAL ( x[0], 1.3+2.5 );

    x[0] = 1.5;
    BOOST_CHECK_EQUAL ( x[0], 1.5 );
}


// check initialiser list
BOOST_AUTO_TEST_CASE( InitializerList )
{
    fdm::SVector<double,3> x = {1.0, 1.1, 1.2};
    BOOST_CHECK_EQUAL ( x[0], 1.0 );
    BOOST_CHECK_EQUAL ( x[1], 1.1 );
    BOOST_CHECK_EQUAL ( x[2], 1.2 );


    fdm::SVector<double,2> y;
    y = fdm::SVector<double,2>({2.0, 2.1});
    BOOST_CHECK_EQUAL ( y[0], 2.0 );
    BOOST_CHECK_EQUAL ( y[1], 2.1 );

    y = {3.0, 3.1};
    BOOST_CHECK_EQUAL ( y[0], 3.0 );
    BOOST_CHECK_EQUAL ( y[1], 3.1 );

    // having too many or too few arguments, this must fail
    // fdm::SVector<double,3> z = {1.0, 1.1};

}


