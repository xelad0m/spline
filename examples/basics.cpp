#include <iostream>
#include <cmath>
#include "spline.hpp"


typedef double T;
using bc = Spline<T>::bc_type;


template<typename Type = T>
void print(const std::vector<Type> & v)
{
    std::cout << "[";
    for (auto & i :  v) 
        std::cout << i << (&i != &v.back() ? ", " :"");
    std::cout << "]\n";
}

template<typename T>
Spline<T> factory(int param) 
{   
    if (param > 0)
        return std::move(Spline<T>(bc::natural, (T)param, bc::natural, (T)param));
    else
        return std::move(Spline<T>(bc::clamped, (T)param, bc::clamped, (T)param));
}

int main(int argc, char** argv) 
{    
    std::vector<T> x;
    std::vector<T> y;

    int n = 3;
    T step = 1.0;
    for (int i = 0; i <= n; ++i) {
        x.push_back(i*step);
        y.push_back(std::exp(i*step));
    }
    std::cout << "Using n = " << n << std::endl;
    std::cout << "Using x = ";
    print(x);
    std::cout << "Using y = ";
    print(y);

    // constructor variants
    Spline<T> s(x, y);
    Spline<T> s1(x, y, 1e-16);
    Spline<T> sc1(x, y, bc::clamped);
    Spline<T> sc2(x, y, bc::clamped, 1);
    Spline<T> sc3(x, y, bc::clamped, 1, bc::clamped, -1);
    Spline<T> sn1(x, y, bc::natural);
    Spline<T> sn2(x, y, bc::natural, 1);
    Spline<T> sn3(x, y, bc::natural, 1, bc::natural, -1);
    Spline<T> sk1(x, y, bc::notaknot);
    Spline<T> sk2(x, y, bc::notaknot, 1);
    Spline<T> sm1(x, y, bc::notaknot, 1, bc::clamped, -1);
    Spline<T> sm2(x, y, bc::natural, 1, bc::clamped, -1);

    // unfitted
    Spline<T> spl(bc::clamped, std::exp(x[0]), bc::clamped, std::exp(x[n]));    

    // fitting/refitting on new data
    spl.fit(x, y);

    // copy variants
    Spline<T> scopy1(spl);
    Spline<T> spline = scopy1;

    // move variants
    Spline<T> s10(factory<T>(1));
    Spline<T> s12 = factory<T>(-1);

    // evaluate
    std::cout << "s(1) = exp(1)            : " << spline(1.0) << " = " << std::exp(1.0) << std::endl;
    std::cout << "s(1.5) -n-> exp(1.5)     : " << spline(1.5) << " -n-> " << std::exp(1.5) << std::endl;

    // deriv
    std::cout << "s'(1) -n-> exp'(1)       : " << spline.deriv(1.0) << " -n-> " << std::exp(1.0) << std::endl;
    std::cout << "s'(1.5) -n-> exp'(1.5)   : " << spline.deriv(1.5) << " -n-> " << std::exp(1.5) << std::endl;
    
    // deriv2
    std::cout << "s''(1) -n-> exp''(1)     : " << spline.deriv2(1.0) << " -n-> " << std::exp(1.0) << std::endl;
    std::cout << "s''(1.5) -n-> exp''(1.5) : " << spline.deriv2(1.5) << " -n-> " << std::exp(1.5) << std::endl;

    // integrate
    T splint = spline.integrate(x[0], x[n]);
    T expint = std::exp(x[n]) - std::exp(x[0]);
    std::cout << "integral(0, 3)           : " << splint << " -n-> " << expint << std::endl;
    
    // params
    std::cout << "Order of parametrized derivative on left:  " << spline.bc_left_type << ", value = " << spline.bc_left_value << std::endl;
    std::cout << "Order of parametrized derivative on right: " << spline.bc_right_type << ", value = " << spline.bc_right_value << std::endl;

    // coefficients
    std::cout << "a = ";
    print(spline.a_());
    std::cout << "b = ";
    print(spline.b_());
    std::cout << "c = ";
    print(spline.c_());
    std::cout << "d = ";
    print(spline.d_());

    std::cout << "\nnot-a-knot\n";

    x = std::vector<T> {0.6, 1.2, 1.8, 3};
    y = std::vector<T> {0.5,   1, 0.8, 2};
    Spline<T> tst(x, y, bc::notaknot);
    // coefficients
    std::cout << "x = ";
    print(tst.x_());
    std::cout << "a = ";
    print(tst.a_());
    std::cout << "b = ";
    print(tst.b_());
    std::cout << "c = ";
    print(tst.c_());
    std::cout << "d = ";
    print(tst.d_());

    return 0;
}