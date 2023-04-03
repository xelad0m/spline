/* spline.hpp
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * SPDX-License-Identifier: GPL-2.0-or-later
 */


#include <vector>
#include <algorithm>    // upper_bound
#include <cmath>        // abs
#include <cassert>

template<typename T = double>
class Spline
{
public:

    enum bc_type { clamped = 1, natural = 2, notaknot = 3 };

    Spline() : bc_left_type(natural), bc_left_value(0.0), bc_right_type(natural), bc_right_value(0.0), tol(1e-8) {};
    
    Spline(bc_type t, T bc_value = 0.0, T tol = 1e-8) 
        : bc_left_type(t), bc_left_value(bc_value), bc_right_type(t), bc_right_value(bc_value), tol(tol) {};
    Spline(bc_type left_type, T bc_left_value, bc_type right_type, T bc_right_value, T tol = 1e-8) 
        : bc_left_type(left_type), bc_left_value(bc_left_value), bc_right_type(right_type), bc_right_value(bc_right_value), tol(tol) {};

    Spline(const std::vector<T> & x, const std::vector<T> & y, T tol = 1e-8) 
        : bc_left_type(natural), bc_left_value(0.0), bc_right_type(natural), bc_right_value(0.0), tol(tol), x(x)
    { fit(y); };

    Spline(const std::vector<T> & x, const std::vector<T> & y, 
            bc_type t, T bc_value = 0.0, T tol = 1e-8) 
        : bc_left_type(t), bc_left_value(bc_value), bc_right_type(t), bc_right_value(bc_value), tol(tol), x(x)
    { fit(y); };

    Spline(const std::vector<T> & x, const std::vector<T> & y, 
            bc_type left_type, T bc_left_value,
            bc_type right_type, T bc_right_value,
            T tol = 1e-8)
        : bc_left_type(left_type), bc_left_value(bc_left_value), 
            bc_right_type(right_type), bc_right_value(bc_right_value), tol(tol), x(x)
    { fit(y); };
      
    Spline(const Spline & s) // copy
        : bc_left_type(s.bc_left_type), bc_left_value(s.bc_left_value), 
            bc_right_type(s.bc_right_type), bc_right_value(s.bc_right_value), tol(s.tol),
            x(s.x), a(s.a), b(s.b), c(s.c), d(s.d) {};
    
    Spline (Spline && s) // move
        : bc_left_type(s.bc_left_type), bc_left_value(s.bc_left_value), 
            bc_right_type(s.bc_right_type), bc_right_value(s.bc_right_value), tol(s.tol),
            x(s.x), a(s.a), b(s.b), c(s.c), d(s.d) {};

    T operator() (T val) const
    {   
        int i = di(val);
        T dx = val - x[i];
        return a[i] + (b[i] + (c[i] + d[i]*dx)*dx)*dx;
    };

    int di (T&) const;
    T deriv(T) const;
    T deriv2(T) const;
    T deriv3(T) const;
    T integrate(T, T) const;
    
    std::vector<T> x_() const { return x; };
    std::vector<T> y_() const { return a; };
    std::vector<T> a_() const { return a; };
    std::vector<T> b_() const { return b; };
    std::vector<T> c_() const { return c; };
    std::vector<T> d_() const { return d; };

    bc_type const   bc_left_type;
    T const         bc_left_value;
    bc_type const   bc_right_type;
    T const         bc_right_value;

    T const         tol;

    void fit(const std::vector<T> &, const std::vector<T> &);
private:  
    void fit(const std::vector<T> &);

    std::vector<T>  x;
    std::vector<T>  a;
    std::vector<T>  b;
    std::vector<T>  c;
    std::vector<T>  d;
};

template<typename T>
void Spline<T>::fit(const std::vector<T> & x_new, const std::vector<T> & y_new) 
{
    x = x_new;
    fit(y_new);
}

template<typename T>
void Spline<T>::fit(const std::vector<T> & y) 
{   
    int n = x.size() - 1;
    assert(n+1 >= 3); // you need at least 3 points to build spline
    assert(n+1 == (int)y.size()); // 'x' and 'y' length mismatch
    for (int i = 1; i <= n; ++i)
        assert( x[i] - x[i-1] >= tol ); // you need a strictly increasing 'x' with an interval of at least 'tol' param (default 1e-8)

    a = std::vector<T>(y);
    b = std::vector<T>(n+1, 0);
    c = std::vector<T>(n+1, 0);
    d = std::vector<T>(n+1, 0);

    std::vector<T> l(n+1, 0), u(n+1, 0), z(n+1, 0), F(n+1, 0), h(n, 0);

    for (int i = 0; i < n; ++i) {
        h[i] = x[i+1] - x[i];
    }

    bool zerocorner = ( bc_left_type == bc_type::notaknot && std::abs(h[0] - h[1]) < tol );

    // left boundary condition
    if ( bc_left_type == bc_type::natural ) {
        l[0] = 1;
        z[0] = bc_left_value / 2.0;
        u[0] = 0;
    } else if ( bc_left_type == bc_type::clamped ) {
        F[0] = 3*(a[1]-a[0])/h[0] - 3*bc_left_value;                    // F[0]
        l[0] = 2*h[0];                                                  // C[0]
        u[0] = 0.5;                                                     // B[0] / l[0]
        z[0] = F[0] / l[0];                                             // F[0] / l[0]
    } else if ( bc_left_type == bc_type::notaknot && !zerocorner ) {
        assert(n+1 >= 4); // you need at least 4 points to make not-a-knot bound
        l[0] = h[0]*h[0] - h[1]*h[1];                                   // C[0]
        l[0] = (l[0]==0)?tol:l[0];
        u[0] = ( 2*h[0]*h[0] + 3*h[0]*h[1] + h[1]*h[1] ) / l[0];        // B[0] / l[0]
        z[0] = ( 3*h[0]*((a[2]-a[1])/h[1] - (a[1]-a[0])/h[0]) ) / l[0]; // F[0] / l[0]
    } else if ( bc_left_type == bc_type::notaknot && zerocorner ) {      // zero corner case
        assert(n+1 >= 4); // you need at least 4 points to make not-a-knot bound
        c[1] = (3*h[0]*((a[2]-a[1])/h[1]-(a[1]-a[0])/h[0])) / (2*h[0]*h[0] + 3*h[0]*h[1] + h[1]*h[1]);  // F[0] / B[0]
        l[1] = h[0];                                                                                    // A[1]
        u[1] = h[1] / l[1];                                                                             // B[1] / l[1]
        z[1] = (3*((a[2]-a[1])/h[1] - (a[1]-a[0])/h[0])) - 2*(h[1]+h[0])*c[1];                          // F[1] - C[1]*c[1]
    } else {
    }

    // forward
    if ( bc_left_type == bc_type::notaknot && zerocorner ) {
        for (int i = 1; i < n; ++i) { // adjusted Crout factorization
            F[i] = 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1];           // F[i]
            if ( i == 2 ) {
                F[i] -= h[1]*c[1];                                          // F[2] - A[2]*c[1]
                l[i] = 2*(x[i+1]-x[i-1]);                                   // C[2] - 0 * u[1]
                u[i] = h[i] / l[i];                                         // B[2] / l[1]
                z[i] = F[i] / l[i];                                         // (F[2] - 0 * z[1]) / l[1]
            } else {
                l[i] = 2*(x[i+1]-x[i-1]) - h[i-1]*u[i-1];                   // C[i] - A[i]*u[i-1]
                u[i] = h[i] / l[i];                                         // B[i] / l[i]
                z[i] = (F[i] - h[i-1]*z[i-1]) / l[i];                       // (F[i] - A[i] * z[i-1]) / l[i]
            }
        }
    } else { // regular Crout factorization
        for (int i = 1; i < n; ++i) {
            F[i] = 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1];           // F[i]
            l[i] = 2*(x[i+1] - x[i-1]) - h[i-1] * u[i-1];                   // C[i] - A[i]*u[i-1]
            u[i] = h[i] / l[i];                                             // B[i] / l[i]
            z[i] = (F[i] - h[i-1] * z[i-1]) / l[i];                         // (F[i] - A[i] * z[i-1]) / l[i]
        }
    }

    // right boundary condition 
    if ( bc_right_type == bc_type::natural) {
        l[n] = 0;
        z[n] = 0;
        c[n] = bc_right_value / 2.0;
    } else if ( bc_right_type == bc_type::clamped ) {
        F[n] = 3 * bc_right_value - 3*(a[n]-a[n-1])/h[n-1];             // F[n]
        l[n] = h[n-1]*(2 - u[n-1]);                                     // C[n] - A[n]*u[n-1]
        z[n] = (F[n] - h[n-1]*z[n-1]) / l[n];                           // (F[n] - A[n]*z[n-1])/l[n]
        c[n] = z[n];
    } else if ( bc_right_type == bc_type::notaknot )  {
        assert(n+1 >= 4); // you need at least 4 points to make not-a-knot bound
        F[n] = 3*h[n-1]*((a[n]-a[n-1])/h[n-1] - (a[n-1]-a[n-2] )/h[n-2]);                                       // F[n] * h[n-1]
        l[n] = (h[n-1]*h[n-1] - h[n-2]*h[n-2]) - (h[n-2]*h[n-2] + 3*h[n-2]*h[n-1] + 2*h[n-1]*h[n-1]) * u[n-1];  // C[n] - A[n]*u[n-1]
        z[n] = (F[n] - (h[n-2]*h[n-2] + 3*h[n-2]*h[n-1] + 2*h[n-1]*h[n-1]) * z[n-1]) / l[n];                    // (F[n] - A[n] * z[n-1]) / l[n]
        c[n] = z[n];
    } else {
    }

    // backward
    if ( bc_left_type == bc_type::notaknot && zerocorner ) {  // in case of zero corner
        for (int j = n-1; j >= 0; j--) {
            if (j > 1) c[j] = z[j] - u[j]*c[j+1];
            if (j == 0) {
                c[0] = ( 3*((a[2]-a[1])/h[1] - (a[1]-a[0])/h[0])  - 2*(h[1]+h[0])*c[1] - h[1]*c[2]) / h[0];    // (F[1]-C[1]*c[1]-B[1]*c[2]) / A[1]
            }
            b[j] = (a[j+1]-a[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3;
            d[j] = (c[j+1]-c[j])/(3*h[j]);
        }
    } else { // regular solve
        for (int j = n-1; j >= 0; j--) {
            c[j] = z[j] - u[j]*c[j+1];
            b[j] = (a[j+1]-a[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3;
            d[j] = (c[j+1]-c[j])/(3*h[j]);
        }
    }
}


template<typename T>
int Spline<T>::di (T & val) const
{
    typename std::vector<T>::const_iterator it;
    it = std::upper_bound(x.begin(), --x.end(), val); 
    int i = std::max(int(it - x.begin()) - 1, 0);
    return i;
}

template<typename T>
T Spline<T>::deriv(T val) const 
{
    int i = di(val);
    T dx = val - x[i];
    return b[i] + (2*c[i] + 3*d[i]*dx)*dx;
}

template<typename T>    
T Spline<T>::deriv2(T val) const 
{
    int i = di(val);
    T dx = val - x[i];
    return 2*c[i] + 6*d[i]*dx;
}

template<typename T>    
T Spline<T>::deriv3(T val) const 
{ 
    return 6*d[di(val)]; 
}

template<typename T>    
T Spline<T>::integrate(T x1, T x2) const
{   
    assert(x1 < x2 || x1 >= x[0] || x2 <= x.back());  // limits x1 > x2 or beyond x0..xn
    
    int i = di(x1), j = di(x2);

    T dx = x1 - x[i];
    T s1 = (a[i] + (b[i]/2 + (c[i]/3 + d[i]/4*dx)*dx)*dx)*dx;

    T s2 = 0;
    for (int k=i; k < j; ++k) {
        dx = x[k+1] - x[k];
        s2 += (a[k] + (b[k]/2 + (c[k]/3 + d[k]/4*dx)*dx)*dx)*dx;
    }
    dx = x2 - x[j];
    s2 += (a[j] + (b[j]/2 + (c[j]/3 + d[j]/4*dx)*dx)*dx)*dx;

    return s2 - s1;
}