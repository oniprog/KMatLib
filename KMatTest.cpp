/////////////////////////////////////////////////////////////////////////////
/** @file
    @brief
    @author oniprog
*/
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"

#include "gtest/gtest.h"

#ifndef NDEBUG
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/io.hpp"
#include "boost/timer.hpp"
#endif

#include "KMat.h"

/////////////////////////////////////////////////////////////////////////////
int main( int argc, char **argv ) {
    
    testing::InitGoogleTest(&argc, argv);
    auto nRetCode = RUN_ALL_TESTS();
    getchar();
    return nRetCode;
}

/////////////////////////////////////////////////////////////////////////////
TEST( TestProd, Test1 ) {

    kblas::KMat<double,1,2> m1;
    kblas::KVec<double,2> v1;

    m1(0,0) = 1;  m1(0,1) = 2;
    v1(0) = 3;  v1(1) = 4;

    auto v2 = kblas::prod(m1,v1);

    EXPECT_NEAR( 1*3+2*4, v2(0), 1E-30 );

#if 0
    boost::numeric::ublas::matrix<double> m10(1,2);
    boost::numeric::ublas::vector<double> v10(2);
    m10(0,0) = 1;   m10(0,1) = 2;
    v10(0) = 3; v10(1) = 4;
    auto v12 = prod(m10,v10);
    std::cout << v12 << std::endl;
#endif
}

/////////////////////////////////////////////////////////////////////////////
TEST( TestProd, Test2 ) {

    kblas::KMat<double,2,1> m1;
    kblas::KVec<double,2> v1;

    m1(0,0) = 1;  m1(1,0) = 2;
    v1(0) = 3;  v1(1) = 4;

    auto v2 = prod(v1,m1);

    EXPECT_NEAR( 1*3+2*4, v2(0), 1E-30 );

#if 0
    boost::numeric::ublas::matrix<double> m10(2,1);
    boost::numeric::ublas::vector<double> v10(2);
    m10(0,0) = 1;   m10(1,0) = 2;
    v10(0) = 3; v10(1) = 4;
    auto v12 = prod(v10,m10);
    std::cout << v12 << std::endl;
#endif
}

/////////////////////////////////////////////////////////////////////////////
TEST( TestProd, Test3 ) {

    kblas::KMat<double,2,1> m1;
    kblas::KMat<double,1,2> m2;

    m1(0,0) = 1;    m1(1,0) = 2;
    m2(0,0) = 3;    m2(0,1) = 4;

    auto m3 = prod(m2,m1);
    EXPECT_NEAR( 1*3+2*4, m3(0,0), 1E-30 );

    auto m4 = prod(m1,m2);
    EXPECT_NEAR( 1*3, m4(0,0), 1E-30 );
    EXPECT_NEAR( 1*4, m4(0,1), 1E-30 );
    EXPECT_NEAR( 2*3, m4(1,0), 1E-30 );
    EXPECT_NEAR( 2*4, m4(1,1), 1E-30 );

#if 0
    boost::numeric::ublas::matrix<double> m10(2,1);
    boost::numeric::ublas::matrix<double> m11(1,2);
    m10(0,0) = 1;   m10(1,0) = 2;
    m11(0,0) = 3;   m11(0,1) = 4;
    auto m12 = prod(m10,m11);
    std::cout << m12 << std::endl;
#endif
}

/////////////////////////////////////////////////////////////////////////////
TEST( TestProd, Test4 ) {

    kblas::KMat<double,2,1> m1;
    kblas::KMat<double,1,3> m2;

    m1(0,0) = 1;    m1(1,0) = 2;
    m2(0,0) = 3;    m2(0,1) = 4;    m2(0,2) = 5;

    auto m3 = prod(m1,m2);
    EXPECT_NEAR( 1*3, m3(0,0), 1E-30 );
    EXPECT_NEAR( 1*4, m3(0,1), 1E-30 );
    EXPECT_NEAR( 1*5, m3(0,2), 1E-30 );
    EXPECT_NEAR( 2*3, m3(1,0), 1E-30 );
    EXPECT_NEAR( 2*4, m3(1,1), 1E-30 );
    EXPECT_NEAR( 2*5, m3(1,2), 1E-30 );

#if 0
    boost::numeric::ublas::matrix<double> m10(2,1);
    boost::numeric::ublas::matrix<double> m11(1,3);
    m10(0,0) = 1;   m10(1,0) = 2;
    m11(0,0) = 3;   m11(0,1) = 4;   m11(0,2) = 5;
    auto m12 = prod(m10,m11);
    std::cout << m12 << std::endl;
#endif
}

/////////////////////////////////////////////////////////////////////////////
TEST( TestTrans, Test1 ) {

    kblas::KMat<double, 2, 2> m1, m2;

    m1(0,0) = 1;    m1(0,1) = 2;    
    m1(1,0) = 3;    m1(1,1) = 4;

    m2(0,0) = 5;    m2(0,1) = 6;
    m2(1,0) = 7;    m2(1,1) = 8;

    auto m3 = prod(trans(m1), m2);

    EXPECT_NEAR(1*5+3*7, m3(0,0), 1E-30 );
    EXPECT_NEAR(1*6+3*8, m3(0,1), 1E-30 );
    EXPECT_NEAR(2*5+4*7, m3(1,0), 1E-30 );
    EXPECT_NEAR(2*6+4*8, m3(1,1), 1E-30 );

#if 0
    boost::numeric::ublas::matrix<double> m10(2,2);
    boost::numeric::ublas::matrix<double> m11(2,2);
    m10(0,0) = 1;   m10(0,1) = 2;
    m10(1,0) = 3;   m10(1,1) = 4;  

    m11(0,0) = 5;   m11(0,1) = 6;
    m11(1,0) = 7;   m11(1,1) = 8;

    auto m12 = prod(trans(m10),m11);
    std::cout << m12 << std::endl;
#endif

}

/////////////////////////////////////////////////////////////////////////////
TEST( TestMMt, Test1 ) {

    kblas::KMat<double, 2, 2> m1, m2;

    m1(0,0) = 1;    m1(0,1) = 2;    
    m1(1,0) = 3;    m1(1,1) = 4;

    m2(0,0) = 5;    m2(0,1) = 6;
    m2(1,0) = 7;    m2(1,1) = 8;

    auto m3 = prod(m1, trans(m2) );

    EXPECT_NEAR(1*5+2*6, m3(0,0), 1E-30 );
    EXPECT_NEAR(1*7+2*8, m3(0,1), 1E-30 );
    EXPECT_NEAR(3*5+4*6, m3(1,0), 1E-30 );
    EXPECT_NEAR(3*7+4*8, m3(1,1), 1E-30 );

#if 0
    boost::numeric::ublas::matrix<double> m10(2,2);
    boost::numeric::ublas::matrix<double> m11(2,2);
    m10(0,0) = 1;   m10(0,1) = 2;
    m10(1,0) = 3;   m10(1,1) = 4;  

    m11(0,0) = 5;   m11(0,1) = 6;
    m11(1,0) = 7;   m11(1,1) = 8;

    auto m12 = prod(trans(m10),m11);
    std::cout << m12 << std::endl;
#endif

}

/////////////////////////////////////////////////////////////////////////////
TEST( TestMMt, Test2 ) {

    kblas::KMat<double, 2, 1> m1;
    kblas::KMat<double, 2, 1> m2;

    m1(0,0) = 1;    m1(1,0) = 2;    
    m2(0,0) = 5;    m2(1,0) = 6;

    auto m3 = prod(m1, trans(m2) );

    EXPECT_NEAR(1*5, m3(0,0), 1E-30 );
    EXPECT_NEAR(1*6, m3(0,1), 1E-30 );
    EXPECT_NEAR(2*5, m3(1,0), 1E-30 );
    EXPECT_NEAR(2*6, m3(1,1), 1E-30 );

#if 0
    boost::numeric::ublas::matrix<double> m10(2,2);
    boost::numeric::ublas::matrix<double> m11(2,2);
    m10(0,0) = 1;   m10(0,1) = 2;
    m10(1,0) = 3;   m10(1,1) = 4;  

    m11(0,0) = 5;   m11(0,1) = 6;
    m11(1,0) = 7;   m11(1,1) = 8;

    auto m12 = prod(trans(m10),m11);
    std::cout << m12 << std::endl;
#endif

}

/////////////////////////////////////////////////////////////////////////////
TEST( TestMtMt, Test1 ) {

    kblas::KMat<double, 2, 2> m1, m2;

    m1(0,0) = 1;    m1(0,1) = 2;    
    m1(1,0) = 3;    m1(1,1) = 4;

    m2(0,0) = 5;    m2(0,1) = 6;
    m2(1,0) = 7;    m2(1,1) = 8;

    auto m3 = prod(trans(m1), trans(m2) );

    EXPECT_NEAR(1*5+3*6, m3(0,0), 1E-30 );
    EXPECT_NEAR(1*7+3*8, m3(0,1), 1E-30 );
    EXPECT_NEAR(2*5+4*6, m3(1,0), 1E-30 );
    EXPECT_NEAR(2*7+4*8, m3(1,1), 1E-30 );

#if 0
    boost::numeric::ublas::matrix<double> m10(2,2);
    boost::numeric::ublas::matrix<double> m11(2,2);
    m10(0,0) = 1;   m10(0,1) = 2;
    m10(1,0) = 3;   m10(1,1) = 4;  

    m11(0,0) = 5;   m11(0,1) = 6;
    m11(1,0) = 7;   m11(1,1) = 8;

    auto m12 = prod(trans(m10),m11);
    std::cout << m12 << std::endl;
#endif

}

/////////////////////////////////////////////////////////////////////////////
TEST( TestMtMt, Test2 ) {

    kblas::KMat<double, 2, 1> m1;
    kblas::KMat<double, 1, 2> m2;

    m1(0,0) = 1;    m1(1,0) = 2;    
    m2(0,0) = 5;    m2(0,1) = 6;

    auto m3 = prod(trans(m1), trans(m2) );

    EXPECT_NEAR(1*5+2*6, m3(0,0), 1E-30 );

#if 0
    boost::numeric::ublas::matrix<double> m10(2,2);
    boost::numeric::ublas::matrix<double> m11(2,2);
    m10(0,0) = 1;   m10(0,1) = 2;
    m10(1,0) = 3;   m10(1,1) = 4;  

    m11(0,0) = 5;   m11(0,1) = 6;
    m11(1,0) = 7;   m11(1,1) = 8;

    auto m12 = prod(trans(m10),m11);
    std::cout << m12 << std::endl;
#endif

}


/////////////////////////////////////////////////////////////////////////////
TEST( TestMtV, Test1 ) {

    kblas::KMat<double, 2, 1> m1;
    kblas::KVec<double, 2> v1;

    m1(0,0) = 1;    m1(0,1) = 2;    
    v1(0) = 5;    v1(1) = 6;

    auto v3 = prod(trans(m1), v1 );

    EXPECT_NEAR(1*5+2*6, v3(0), 1E-30 );

#if 0
    boost::numeric::ublas::matrix<double> m10(2,2);
    boost::numeric::ublas::matrix<double> m11(2,2);
    m10(0,0) = 1;   m10(0,1) = 2;
    m10(1,0) = 3;   m10(1,1) = 4;  

    m11(0,0) = 5;   m11(0,1) = 6;
    m11(1,0) = 7;   m11(1,1) = 8;

    auto m12 = prod(trans(m10),m11);
    std::cout << m12 << std::endl;
#endif

}

/////////////////////////////////////////////////////////////////////////////
TEST( TestAddMM, Test1 ) {

    kblas::KMat<double,1,2> m1, m2;
    m1(0,0) = 1;    m1(0,1) = 2;
    m2(0,0) = 3;    m2(0,1) = 4;

    m1 += m2;

    EXPECT_NEAR( 4, m1(0,0), 1E-10 );
    EXPECT_NEAR( 6, m1(0,1), 1E-10 );
}

