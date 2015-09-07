#include <fdm/grid.hpp>

#include <cstdio>
#include <cstdlib>
#include <vector>

#include <fdm/svector.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GridTest
#include <boost/test/unit_test.hpp>
#include <boost/random.hpp>


double runif()
{
    static boost::mt19937   rng;
    static boost::uniform_01<boost::mt19937&> myrand(rng);
    return myrand();
}


class MappingQuadratic
{
private:
    double  m_a, m_b;
public:
    MappingQuadratic(double a, double b): m_a(a), m_b(b) {};
    double lower() const
    {
        return m_a;
    }
    double upper() const
    {
        return m_b;
    }
    double operator() (double x) const
    {
        return  m_a + (m_b-m_a)*x*x;
    }
};


BOOST_AUTO_TEST_CASE( Grid )
{
    std::vector<double> X = {-2.3, -1.2, 0.15, 7.3, 125.0};
    MappingQuadratic    g(-5.0, 5.0);
    fdm::Grid<4> grid;

    // set the grid in each dimension individually using a different method
    grid.coord(0) = { 1.1*X[0], 1.1*X[1], 1.1*X[2], 1.1*X[3], 1.1*X[4] };
    grid.coord(1) = X;
    grid.set_uniform(2, -10.0, 10.0, 21);
    grid.set_by_mapping(3, g, 40);

    // check size is right
    BOOST_CHECK_EQUAL ( grid.dim(), 4 );
    BOOST_CHECK_EQUAL ( grid.size(0), grid.coord(0).size() );
    BOOST_CHECK_EQUAL ( grid.size(0), X.size() );
    BOOST_CHECK_EQUAL ( grid.size(1), X.size() );
    BOOST_CHECK_EQUAL ( grid.size(2), 21 );
    BOOST_CHECK_EQUAL ( grid.size(3), 40 );

    // check correct coordinates
    for(size_t i=0; i<grid.size(0); i++) {
        BOOST_CHECK_EQUAL ( grid.coord(0,i), grid.coord(0)[i] );
        BOOST_CHECK_EQUAL ( grid.coord(0,i), 1.1*X[i] );
    }
    for(size_t i=0; i<grid.size(1); i++)
        BOOST_CHECK_EQUAL ( grid.coord(1,i), X[i] );
    for(size_t i=0; i<grid.size(2); i++)
        BOOST_CHECK_EQUAL ( grid.coord(2,i), -10.0 + (double) i );
    for(size_t i=0; i<grid.size(3); i++) {
        double dx = 1.0/(grid.size(3)-1);
        double x = dx*i;
        BOOST_CHECK_EQUAL ( grid.coord(3,i), g(x) );
    }

    // find_nearest
    BOOST_CHECK_EQUAL ( grid.find_nearest(0,1.2), 2 );
    BOOST_CHECK_EQUAL ( grid.find_nearest(1,-20.3), 0 );
    BOOST_CHECK_EQUAL ( grid.find_nearest(1,-2.0), 0 );
    BOOST_CHECK_EQUAL ( grid.find_nearest(1,-1.7), 1 );

    // on a uniform grid we know exactly the result of find_nearest
    // to simplify we reset grid in dim 3
    grid.set_uniform(3, 0.0, 1.0, 101);
    {
        size_t n=10000;
        double a=-1.0, b=2.0;
        double dx=(b-a)/(n-1);
        for(size_t i=0; i<n; i++) {
            double x=a+dx*i;
            int idx= std::min(std::max((int) round(100.0*x), 0), 100);
            BOOST_CHECK_EQUAL ( grid.find_nearest(3,x), idx );
        }
    }

    // coordinates need to be strictly monotonic increasing
    // confirm that valid_coord() finds any violations
    BOOST_CHECK_EQUAL ( grid.valid_coord(0), true  );
    BOOST_CHECK_EQUAL ( grid.valid_coord(1), true  );
    BOOST_CHECK_EQUAL ( grid.valid_coord(2), true  );
    BOOST_CHECK_EQUAL ( grid.valid_coord(3), true  );
    grid.coord(0,3) = grid.coord(0,2);
    BOOST_CHECK_EQUAL ( grid.valid_coord(0), false );
    grid.coord(0,3) += 1e-5 ;
    BOOST_CHECK_EQUAL ( grid.valid_coord(0), true  );
    grid.coord(0,1) = -4.2;
    BOOST_CHECK_EQUAL ( grid.valid_coord(0), false );

}


BOOST_AUTO_TEST_CASE( GridC2 )
{
}



BOOST_AUTO_TEST_CASE( GridIterator )
{
    fdm::GridC2<4> grid;
    grid.set_uniform(0, 0.0, 1.0, 13);
    grid.set_uniform(1, 0.0, 1.0, 7);
    grid.set_uniform(2, 0.0, 1.0, 16);
    grid.set_uniform(3, 0.0, 1.0, 5);
    fdm::GridFunction<4> g(grid);
    fdm::SVector<int, 4> idx;
    fdm::SVector<int, 4> shape=g.shape();

    BOOST_CHECK_EQUAL ( shape[0], 13 );
    BOOST_CHECK_EQUAL ( shape[1],  7 );
    BOOST_CHECK_EQUAL ( shape[2], 16 );
    BOOST_CHECK_EQUAL ( shape[3],  5 );


    // iterator test, going through all points
    g.iterator.begin();
    while(g.iterator.isend()==false) {
        // check embedded index ind() and ind(idx) are the same
        idx=g.iterator.idx();
        BOOST_CHECK_EQUAL ( g.iterator.ind(), g.iterator.ind(idx) );
        // set grid values
        g.value()=1.3;
        g.iterator++;
    }
    for(idx[0]=0; idx[0]<shape[0]; idx[0]++) {
        for(idx[1]=0; idx[1]<shape[1]; idx[1]++) {
            for(idx[2]=0; idx[2]<shape[2]; idx[2]++) {
                for(idx[3]=0; idx[3]<shape[3]; idx[3]++) {
                    BOOST_CHECK_EQUAL ( g.value(idx), 1.3 );
                }
            }
        }
    }

    // iterator test, going through all inner points
    g.iterator.begin_inner();
    while(g.iterator.isend()==false) {
        // check embedded index ind() and ind(idx) are the same
        idx=g.iterator.idx();
        BOOST_CHECK_EQUAL ( g.iterator.ind(), g.iterator.ind(idx) );
        BOOST_CHECK_EQUAL ( g.iterator.check_boundary(), -1 );
        // set grid values
        g.value()=0.3;
        g.iterator.incr_inner();
    }
    for(idx[0]=0; idx[0]<shape[0]; idx[0]++) {
        for(idx[1]=0; idx[1]<shape[1]; idx[1]++) {
            for(idx[2]=0; idx[2]<shape[2]; idx[2]++) {
                for(idx[3]=0; idx[3]<shape[3]; idx[3]++) {
                    if( idx[0]==0 || idx[1]==0 || idx[2]==0 || idx[3]==0
                            || idx[0]==shape[0]-1 || idx[1]==shape[1]-1
                            || idx[2]==shape[2]-1 || idx[3]==shape[3]-1 ) {
                        // boundary point
                        BOOST_CHECK_EQUAL ( g.value(idx), 1.3 );
                        g.iterator.set_idx(idx);
                        BOOST_CHECK_GE ( g.iterator.check_boundary(), 0 );
                    } else {
                        // inner point
                        BOOST_CHECK_EQUAL ( g.value(idx), 0.3 );
                        g.iterator.set_idx(idx);
                        BOOST_CHECK_EQUAL ( g.iterator.check_boundary(), -1 );
                    }
                }
            }
        }
    }


    // iterator test, going through all boundary points
    g.iterator.begin();
    while(g.iterator.isend()==false) {
        // check embedded index ind() and ind(idx) are the same
        idx=g.iterator.idx();
        BOOST_CHECK_EQUAL ( g.iterator.ind(), g.iterator.ind(idx) );
        BOOST_CHECK_GE ( g.iterator.check_boundary(), 0 );
        // set grid values
        g.value()=7.2;
        g.iterator.incr_bound();
    }
    for(idx[0]=0; idx[0]<shape[0]; idx[0]++) {
        for(idx[1]=0; idx[1]<shape[1]; idx[1]++) {
            for(idx[2]=0; idx[2]<shape[2]; idx[2]++) {
                for(idx[3]=0; idx[3]<shape[3]; idx[3]++) {
                    if( idx[0]==0 || idx[1]==0 || idx[2]==0 || idx[3]==0
                            || idx[0]==shape[0]-1 || idx[1]==shape[1]-1
                            || idx[2]==shape[2]-1 || idx[3]==shape[3]-1 ) {
                        // boundary point
                        BOOST_CHECK_EQUAL ( g.value(idx), 7.2 );
                        g.iterator.set_idx(idx);
                        BOOST_CHECK_GE ( g.iterator.check_boundary(), 0 );
                    } else {
                        // inner point
                        BOOST_CHECK_EQUAL ( g.value(idx), 0.3 );
                        g.iterator.set_idx(idx);
                        BOOST_CHECK_EQUAL ( g.iterator.check_boundary(), -1 );
                    }
                }
            }
        }
    }
}


BOOST_AUTO_TEST_CASE( GridFunctionView )
{
    fdm::GridC2<3> grid;
    grid.set_uniform(0, 0.0, 1.0, 15);
    grid.set_uniform(1, 0.0, 1.0, 8);
    grid.set_uniform(2, 0.0, 1.0, 23);

    fdm::GridFunction<3>    u(grid);

    u.iterator.begin();
    while(u.iterator.isend()==false) {
        u.value()=runif();
        u.iterator++;
    }

    fdm::GridFunctionView<2> v(u);
    fdm::SVector<int,3>     idx;
    fdm::SVector<int,2>     idx2d;
    fdm::SVector<int,1>     idx1d;

    for(size_t i=0; i<grid.size(2); i++) {
        idx1d[0]=i;
        idx[2]=i;
        v.reassign_view(u,idx1d);
        // check the slice of u and v are identical
        for(size_t j=0; j<grid.size(0); j++) {
            for(size_t k=0; k<grid.size(1); k++) {
                idx2d[0] = j;
                idx[0]   = j;
                idx2d[1] = k;
                idx[1]   = k;
                BOOST_CHECK_EQUAL(v.coord(0,idx2d[0]), u.coord(0,idx[0]));
                BOOST_CHECK_EQUAL(v.coord(1,idx2d[1]), u.coord(1,idx[1]));
                BOOST_CHECK_EQUAL(v.value(idx2d), u.value(idx));
                BOOST_CHECK_EQUAL(v.f_1(0,idx2d), u.f_1(0,idx));
                BOOST_CHECK_EQUAL(v.f_1(1,idx2d), u.f_1(1,idx));
                BOOST_CHECK_EQUAL(v.f_2(0,0,idx2d), u.f_2(0,0,idx));
                BOOST_CHECK_EQUAL(v.f_2(0,1,idx2d), u.f_2(0,1,idx));
                BOOST_CHECK_EQUAL(v.f_2(1,0,idx2d), u.f_2(1,0,idx));
                BOOST_CHECK_EQUAL(v.f_2(1,1,idx2d), u.f_2(1,1,idx));
            }
        }

        // check iterator of v goes through all values
        v.iterator.begin();
        while(v.iterator.isend()==false) {
            double x0=v.coord(0);
            double x1=v.coord(1);
            v.value()=x0*100.0+x1-(double)idx[2];
            v.iterator++;
        }
        for(size_t j=0; j<grid.size(0); j++) {
            for(size_t k=0; k<grid.size(1); k++) {
                idx2d[0] = j;
                idx[0]   = j;
                idx2d[1] = k;
                idx[1]   = k;
                double x0=u.coord(0,idx[0]);
                double x1=v.coord(1,idx[1]);
                BOOST_CHECK_EQUAL(v.value(idx2d), u.value(idx));
                BOOST_CHECK_EQUAL(u.value(idx),x0*100.0+x1-(double)idx[2]);
            }
        }

    }

}

