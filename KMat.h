/////////////////////////////////////////////////////////////////////////////
/** @file
    @brief
    @author oniprog
*/
/////////////////////////////////////////////////////////////////////////////

#pragma once

namespace kblas {

// ベクトルクラス 
// T 型
// N ベクトルサイズ
template<class T, int N>
class KVec {
public:
    static const int SIZE = N;
public:
    T & operator()(int i) {
        return m_v[i];
    }
    const T operator()(int i) const {
        return m_v[i];
    }
private:
    T   m_v[N];
};

template<class T, int N, int M>
class KMatTrans;

// 行列クラス
template<class T, int N, int M>
// T 型
// N 行サイズ
// M 列サイズ
class KMat{
public:
    static const int SIZE_X = M;
    static const int SIZE_Y = N;
public:

    KMat() {}

    KMat( const KMatTrans<T,N,M> &m1 ) {
        for( int i=0; i<N; ++i ) for( int j=0; j<M; ++j ) {
            (*this)(i,j) = m1(j,i);
        }
    }

    T & operator()(int i, int j) {
        return m_v[i*M+j];
    }
    const T operator()(int i, int j) const {
        return m_v[i*M+j];
    }

private:
    T   m_v[N*M];

    friend class KMatTrans<T,M,N>;
};

// 転置行列クラス
template<class T, int N, int M>
class KMatTrans {
public:
    static const int SIZE_X = N;
    static const int SIZE_Y = M;
public:

    KMatTrans() {}

    KMatTrans( const KMat<T,M,N> &m ) {
        for(int i=0; i<N*M; ++i) m_v[i] = m.m_v[i];
    }
    KMatTrans( KMat<T,M,N> &&m ) {
        m_v = std::move(m.m_v);
    }

    T & operator()(int i, int j) {
        return m_v[j*N+i];
    }
    const T operator()(int i, int j) const {
        return m_v[j*N+i];
    }

private:
    T   m_v[N*M];
};

namespace Detail {

    ///////////////////////////////////////////////////////////////////////////////////
    // M V の計算
    // かけ算補助クラス
    template<class T, int M, int N, int j, int i>
    struct MultL1 {
        static T f( const KMat<T,M,N> &m1, const KVec<T,N> &v1 ) {
            return MultL1<T, M, N, j, i-1>::f(m1, v1 ) + m1(j, i) * v1(i);
        }
    };

    template<class T, int M, int N, int j>
    struct MultL1<T, M, N, j, -1> {
        static T f( const KMat<T,M,N> &m1, const KVec<T,N> &v1 ) {
            return T();
        }
    };

    // かけ算クラス
    template<class T, int M, int N, int j>
    struct MultL {
        static void f( KVec<T,N> &vr, const KMat<T,M,N> &m1, const KVec<T,N> &v1 ) {
            vr(j) = MultL1<T, M, N, j, N-1>::f( m1, v1 );
            MultL<T,M,N,j-1>::f( vr, m1, v1 );
        }
    };

    template<class T, int M, int N>
    struct MultL<T, M, N, -1> {
        static void f( KVec<T,N> &vr, const KMat<T,M,N> &m1, const KVec<T,N> &v1 ) {
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////
    // Mt V の計算
    namespace MtV {
        
        // かけ算補助クラス
        template<class T, int M, int N, int j, int i>
        struct MultL1 {
            static T f( const KMatTrans<T,M,N> &m1, const KVec<T,N> &v1 ) {
                return MultL1<T, M, N, j, i-1>::f(m1, v1 ) + m1(j, i) * v1(i);
            }
        };

        template<class T, int M, int N, int j>
        struct MultL1<T, M, N, j, -1> {
            static T f( const KMatTrans<T,M,N> &m1, const KVec<T,N> &v1 ) {
                return T();
            }
        };

        // かけ算クラス
        template<class T, int M, int N, int j>
        struct MultL {
            static void f( KVec<T,N> &vr, const KMatTrans<T,M,N> &m1, const KVec<T,N> &v1 ) {
                vr(j) = MultL1<T, M, N, j, N-1>::f( m1, v1 );
                MultL<T,M,N,j-1>::f( vr, m1, v1 );
            }
        };

        template<class T, int M, int N>
        struct MultL<T, M, N, -1> {
            static void f( KVec<T,N> &vr, const KMatTrans<T,M,N> &m1, const KVec<T,N> &v1 ) {
            }
        };
    }

    ///////////////////////////////////////////////////////////////////////////////////
    // V M の計算
    // かけ算補助クラス
    template<class T, int M, int N, int j, int i>
    struct MultR1 {
        static T f( const KVec<T,M> &v1, const KMat<T,M,N> &m1 ) {
            return MultR1<T, M, N, j, i-1>::f( v1, m1 ) + v1(i) * m1(i, j);
        }
    };

    template<class T, int M, int N, int j>
    struct MultR1<T, M, N, j, -1> {
        static T f( const KVec<T,M> &v1, const KMat<T,M,N> &m1 ) {
            return T();
        }
    };

    // かけ算クラス
    template<class T, int M, int N, int j>
    struct MultR {
        static void f( KVec<T,N> &vr, const KVec<T,M> &v1, const KMat<T,M,N> &m1 ) {
            vr(j) = MultR1<T, M, N, j, M-1>::f( v1, m1 );
            MultR<T,M,N,j-1>::f( vr, v1, m1 );
        }
    };

    template<class T, int M, int N>
    struct MultR<T, M, N, -1> {
        static void f( KVec<T,N> &vr, const KVec<T,M> &v1, const KMat<T,M,N> &m1 ) {
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////
    // V Mt の計算
    // かけ算補助クラス
    namespace VMt {
        template<class T, int M, int N, int j, int i>
        struct MultR1 {
            static T f( const KVec<T,M> &v1, const KMat<T,M,N> &m1 ) {
                return MultR1<T, M, N, j, i-1>::f( v1, m1 ) + v1(i) * m1(i, j);
            }
        };

        template<class T, int M, int N, int j>
        struct MultR1<T, M, N, j, -1> {
            static T f( const KVec<T,M> &v1, const KMat<T,M,N> &m1 ) {
                return T();
            }
        };

        // かけ算クラス
        template<class T, int M, int N, int j>
        struct MultR {
            static void f( KVec<T,N> &vr, const KVec<T,M> &v1, const KMat<T,M,N> &m1 ) {
                vr(j) = MultR1<T, M, N, j, M-1>::f( v1, m1 );
                MultR<T,M,N,j-1>::f( vr, v1, m1 );
            }
        };

        template<class T, int M, int N>
        struct MultR<T, M, N, -1> {
            static void f( KVec<T,N> &vr, const KVec<T,M> &v1, const KMat<T,M,N> &m1 ) {
            }
        };
    }

    ///////////////////////////////////////////////////////////////////////////////////
    // M x M 
    // 行列同士の積 Sub2
    template<class T, int M, int N, int O, int i, int j, int k>
    struct MultM3 {
        static T f( const KMat<T,M,N> &m1, const KMat<T,N,O> &m2 ) {
            return MultM3<T,M,N,O,i,j,k-1>::f(m1, m2) + m1(i,k) * m2(k,j);
        }
    };

    template<class T, int M, int N, int O, int i, int j>
    struct MultM3<T,M,N,O,i,j,-1> {
        static T f( const KMat<T,M,N> &m1, const KMat<T,N,O> &m2 ) {
            return T();
        }
    };

    /// 行列同士の積 Sub1
    template<class T, int M, int N, int O, int i, int j>
    struct MultM2 {
        static void f( KMat<T,M,O> &rm, const KMat<T,M,N> &m1, const KMat<T,N,O> &m2 ) {
            rm(i,j) = MultM3<T, M, N, O, i, j, N-1>::f(m1, m2);
            MultM2<T,M,N,O,i,j-1>::f(rm, m1, m2);
        }
    };

    template<class T, int M, int N, int O, int i>
    struct MultM2<T,M,N,O,i,-1> {
        static void f( KMat<T,M,O> &rm, const KMat<T,M,N> &m1, const KMat<T,N,O> &m2 ) {
        }
    };

    /// 行列同士の積
    template<class T, int M, int N, int O, int i>
    struct MultM {
        static void f( KMat<T,M,O> &rm, const KMat<T,M,N> &m1, const KMat<T,N,O> &m2 ) {
            MultM2<T,M,N,O,i,O-1>::f(rm, m1, m2);
            MultM<T,M,N,O,i-1>::f(rm,m1,m2);
        }
    };

    template<class T, int M, int N, int O>
    struct MultM<T,M,N,O,-1> {
        static void f( KMat<T,M,O> &rm, const KMat<T,M,N> &m1, const KMat<T,N,O> &m2 ) {
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////
    // Mt x M  // 上のやつのコピー．受け取る引数が違うだけ
    namespace MtM {

        // 行列同士の積 Sub2
        template<class T, int M, int N, int O, int i, int j, int k>
        struct MultM3 {
            static T f( const KMatTrans<T,M,N> &m1, const KMat<T,N,O> &m2 ) {
                return MultM3<T,M,N,O,i,j,k-1>::f(m1, m2) + m1(i,k) * m2(k,j);
            }
        };

        template<class T, int M, int N, int O, int i, int j>
        struct MultM3<T,M,N,O,i,j,-1> {
            static T f( const KMatTrans<T,M,N> &m1, const KMat<T,N,O> &m2 ) {
                return T();
            }
        };

        /// 行列同士の積 Sub1
        template<class T, int M, int N, int O, int i, int j>
        struct MultM2 {
            static void f( KMat<T,M,O> &rm, const KMatTrans<T,M,N> &m1, const KMat<T,N,O> &m2 ) {
                rm(i,j) = MultM3<T, M, N, O, i, j, N-1>::f(m1, m2);
                MultM2<T,M,N,O,i,j-1>::f(rm, m1, m2);
            }
        };

        template<class T, int M, int N, int O, int i>
        struct MultM2<T,M,N,O,i,-1> {
            static void f( KMat<T,M,O> &rm, const KMatTrans<T,M,N> &m1, const KMat<T,N,O> &m2 ) {
            }
        };

        /// 行列同士の積
        template<class T, int M, int N, int O, int i>
        struct MultM {
            static void f( KMat<T,M,O> &rm, const KMatTrans<T,M,N> &m1, const KMat<T,N,O> &m2 ) {
                MultM2<T,M,N,O,i,O-1>::f(rm, m1, m2);
                MultM<T,M,N,O,i-1>::f(rm,m1,m2);
            }
        };

        template<class T, int M, int N, int O>
        struct MultM<T,M,N,O,-1> {
            static void f( KMat<T,M,O> &rm, const KMatTrans<T,M,N> &m1, const KMat<T,N,O> &m2 ) {
            }
        };
    }

    ///////////////////////////////////////////////////////////////////////////////////
    // M x Mt 
    namespace MMt {

        // 行列同士の積 Sub2
        template<class T, int M, int N, int O, int i, int j, int k>
        struct MultM3 {
            static T f( const KMat<T,M,N> &m1, const KMatTrans<T,N,O> &m2 ) {
                return MultM3<T,M,N,O,i,j,k-1>::f(m1, m2) + m1(i,k) * m2(k,j);
            }
        };

        template<class T, int M, int N, int O, int i, int j>
        struct MultM3<T,M,N,O,i,j,-1> {
            static T f( const KMat<T,M,N> &m1, const KMatTrans<T,N,O> &m2 ) {
                return T();
            }
        };

        /// 行列同士の積 Sub1
        template<class T, int M, int N, int O, int i, int j>
        struct MultM2 {
            static void f( KMat<T,M,O> &rm, const KMat<T,M,N> &m1, const KMatTrans<T,N,O> &m2 ) {
                rm(i,j) = MultM3<T, M, N, O, i, j, N-1>::f(m1, m2);
                MultM2<T,M,N,O,i,j-1>::f(rm, m1, m2);
            }
        };

        template<class T, int M, int N, int O, int i>
        struct MultM2<T,M,N,O,i,-1> {
            static void f( KMat<T,M,O> &rm, const KMat<T,M,N> &m1, const KMatTrans<T,N,O> &m2 ) {
            }
        };

        /// 行列同士の積
        template<class T, int M, int N, int O, int i>
        struct MultM {
            static void f( KMat<T,M,O> &rm, const KMat<T,M,N> &m1, const KMatTrans<T,N,O> &m2 ) {
                MultM2<T,M,N,O,i,O-1>::f(rm, m1, m2);
                MultM<T,M,N,O,i-1>::f(rm,m1,m2);
            }
        };

        template<class T, int M, int N, int O>
        struct MultM<T,M,N,O,-1> {
            static void f( KMat<T,M,O> &rm, const KMat<T,M,N> &m1, const KMatTrans<T,N,O> &m2 ) {
            }
        };
    }
    ///////////////////////////////////////////////////////////////////////////////////
    // Mt x Mt 
    namespace MtMt {

        // 行列同士の積 Sub2
        template<class T, int M, int N, int O, int i, int j, int k>
        struct MultM3 {
            static T f( const KMatTrans<T,M,N> &m1, const KMatTrans<T,N,O> &m2 ) {
                return MultM3<T,M,N,O,i,j,k-1>::f(m1, m2) + m1(i,k) * m2(k,j);
            }
        };

        template<class T, int M, int N, int O, int i, int j>
        struct MultM3<T,M,N,O,i,j,-1> {
            static T f( const KMatTrans<T,M,N> &m1, const KMatTrans<T,N,O> &m2 ) {
                return T();
            }
        };

        /// 行列同士の積 Sub1
        template<class T, int M, int N, int O, int i, int j>
        struct MultM2 {
            static void f( KMat<T,M,O> &rm, const KMatTrans<T,M,N> &m1, const KMatTrans<T,N,O> &m2 ) {
                rm(i,j) = MultM3<T, M, N, O, i, j, N-1>::f(m1, m2);
                MultM2<T,M,N,O,i,j-1>::f(rm, m1, m2);
            }
        };

        template<class T, int M, int N, int O, int i>
        struct MultM2<T,M,N,O,i,-1> {
            static void f( KMat<T,M,O> &rm, const KMatTrans<T,M,N> &m1, const KMatTrans<T,N,O> &m2 ) {
            }
        };

        /// 行列同士の積
        template<class T, int M, int N, int O, int i>
        struct MultM {
            static void f( KMat<T,M,O> &rm, const KMatTrans<T,M,N> &m1, const KMatTrans<T,N,O> &m2 ) {
                MultM2<T,M,N,O,i,O-1>::f(rm, m1, m2);
                MultM<T,M,N,O,i-1>::f(rm,m1,m2);
            }
        };

        template<class T, int M, int N, int O>
        struct MultM<T,M,N,O,-1> {
            static void f( KMat<T,M,O> &rm, const KMatTrans<T,M,N> &m1, const KMatTrans<T,N,O> &m2 ) {
            }
        };
    }
}


///////////////////////////////////////////////////////////////////////////////////
// M V の積
template<class T, int M, int N>
KVec<T, N>   prod( const KMat<T, M, N > &m1, const KVec<T,N> &v1 ) { 
    KVec<T,N> rv;
    Detail::MultL<T, M, N, M-1>::f(rv, m1, v1);
    return rv;
}

///////////////////////////////////////////////////////////////////////////////////
// Mt V の積
template<class T, int M, int N>
KVec<T, N>   prod( const KMatTrans<T, M, N > &m1, const KVec<T,N> &v1 ) { 
    KVec<T,N> rv;
    Detail::MtV::MultL<T, M, N, M-1>::f(rv, m1, v1);
    return rv;
}

///////////////////////////////////////////////////////////////////////////////////
// V Mの積
template<class T, int M, int N>
KVec<T, N>   prod( const KVec<T,M> &v1, const KMat<T, M, N > &m1 ) { 
    KVec<T,N> rv;
    Detail::MultR<T, M, N, N-1>::f(rv, v1, m1 );
    return rv;
}

///////////////////////////////////////////////////////////////////////////////////
// V Mtの積
template<class T, int M, int N>
KVec<T, N>   prod( const KVec<T,M> &v1, const KMatTrans<T, M, N > &m1 ) { 
    KVec<T,N> rv;
    Detail::VMt::MultR<T, M, N, N-1>::f(rv, v1, m1 );
    return rv;
}

///////////////////////////////////////////////////////////////////////////////////
// M Mの積
template<class T, int M, int N, int O>
KMat<T, M, O> prod( const KMat<T, M, N> &m1, const KMat<T, N, O> &m2 ) {
    KMat<T,M,O> rm;
    Detail::MultM<T, M, N, O, M-1>::f(rm, m1, m2 );
    return rm;
}

///////////////////////////////////////////////////////////////////////////////////
// Mt Mの積
template<class T, int M, int N, int O>
KMat<T, M, O> prod( const KMatTrans<T, M, N> &m1, const KMat<T, N, O> &m2 ) {
    KMat<T,M,O> rm;
    Detail::MtM::MultM<T, M, N, O, M-1>::f(rm, m1, m2 );
    return rm;
}

///////////////////////////////////////////////////////////////////////////////////
// M Mtの積
template<class T, int M, int N, int O>
KMat<T, M, O> prod( const KMat<T, M, N> &m1, const KMatTrans<T, N, O> &m2 ) {
    KMat<T,M,O> rm;
    Detail::MMt::MultM<T, M, N, O, M-1>::f(rm, m1, m2 );
    return rm;
}

///////////////////////////////////////////////////////////////////////////////////
// Mt Mtの積
template<class T, int M, int N, int O>
KMat<T, M, O> prod( const KMatTrans<T, M, N> &m1, const KMatTrans<T, N, O> &m2 ) {
    KMat<T,M,O> rm;
    Detail::MtMt::MultM<T, M, N, O, M-1>::f(rm, m1, m2 );
    return rm;
}

///////////////////////////////////////////////////////////////////////////////////
/// trans
template<class T, int M, int N>
KMatTrans<T, N, M> trans( const KMat<T, M, N> &m1 ) {
    return KMatTrans<T,N,M>( m1 );
}

///////////////////////////////////////////////////////////////////////////////////
/// trans
template<class T, int M, int N>
KMatTrans<T, N, M> trans( KMat<T, M, N> &&m1 ) {
    return KMatTrans<T,N,M>( std::move(m1) );
}


///////////////////////////////////////////////////////////////////////////////////
/// trans
template<class T, int M, int N>
KMat<T, N, M> trans( const KMatTrans<T, M, N> &m1 ) {
    return KMat<T,N,M>( m1 );
}

///////////////////////////////////////////////////////////////////////////////////
/// trans
template<class T, int M, int N>
KMat<T, N, M> trans( KMatTrans<T, M, N> &&m1 ) {
    return KMat<T,N,M>( std::move(m1) );
}

} // namespace kblas

/////////////////////////////////////////////////////////////////////////////
