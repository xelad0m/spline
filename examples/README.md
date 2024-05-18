## Examples

1. Basic usage example

```
make basics
./basics
...
```
Evaluating of cubic spline for $y = e^x$ with $x \in \{0, 1, 2, 3\}$.

2. Comparison with [**spline.h**](https://github.com/ttk592/spline) (`tk` in benchmark below) library and [**ALGLIB**](https://www.alglib.net/) implementation of cubic splines. Used bechmark tool from [**spline.h**](https://github.com/ttk592/spline) project

```
./bench 2400 50 1000
                                tk                      alglib                  this
random access:   loops=1e+06,   0.043s ( 102 cycl)      0.058s ( 140 cycl)      0.043s ( 104 cycl)
spline creation: loops=2e+04,   0.109s (1.3e+04 cycl)   0.122s (1.5e+04 cycl)   0.042s (4981 cycl)
grid transform:  loops=2e+04,   0.121s (1.4e+04 cycl)   0.110s (1.3e+04 cycl)   0.052s (6277 cycl)

tk vs ALGLIB accuracy: max difference = 6.66e-15, l2-norm difference = 4.17e-19
this vs ALGLIB accuracy: max difference = 4.44e-15, l2-norm difference = 3.35e-19
```
To build bechmark app you need to have alglib installed kind of this way (package name depends on linux distro, it could be `alglib`, `libalglib`, `libalglib3` and so on)

```
sudo apt-get install libalglib
make bench
```

3. Comparison with [**SciPy**](https://scipy.org/) (`scipy.interpolate`) cubic spline implementations:

- Fortran based `splrep`+`splev`
- Cython based `make_interp_spline`
- Using simple `ctypes` wrapper for this implementation 
- [Results](./tests.ipynb) of this implementation are 3-6x faster than Fortran and Cython implementations in average:

|implementation|build|evaluate| build & evaluate|
|-|:-:|:-:|:-:|
|this|8.9 µs|14.8 µs|23.6 µs
`splrep`+`splev`|10.3 µs|61.7 µs|70.6 µs
`make_interp_spline`|73 µs|61.1 µs|141 µs

