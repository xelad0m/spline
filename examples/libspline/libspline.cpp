/* Python ctypes interface wrapper for spline.hpp library*/

#include "spline.hpp"

typedef double T;
using bc = Spline<T>::bc_type;

extern "C" {

    Spline<T> * spline(T * xarr, T * yarr, int size, int lbc, T lbc_val, int rbc, T rbc_val, T tol) {
        std::vector<T> x(size), y(size);
        for (int i=0; i < size; ++i) {
            x[i] = xarr[i];
            y[i] = yarr[i];
        }
        return new Spline<T>(x, y, (bc)lbc, lbc_val, (bc)rbc, rbc_val, tol);
    }

    void delete_spline(Spline<T> * s) { delete s; }
    
    T eval(Spline<T> * s, T x) { return (*s)(x); }

    void eval_vec(Spline<T> * s, T * x, T * res, int size) { 
        for (int i=0; i < size; ++i) {
            res[i] = (*s)(x[i]);
        }
    }

    T deriv(Spline<T> * s, T x) { return s->deriv(x); }
    T deriv2(Spline<T> * s, T x) { return s->deriv2(x); }
    T deriv3(Spline<T> * s, T x) { return s->deriv3(x); }
    T integrate(Spline<T> * s, T x1, T x2)  { return s->integrate(x1, x2); }

}