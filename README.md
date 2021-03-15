## gram_savitzky_golay

![CI of gram-savitzy-golay](https://github.com/arntanguy/gram_savitzky_golay/workflows/CI%20of%20gram-savitzy-golay/badge.svg)
[![Documentation](https://img.shields.io/badge/website-online-brightgreen?logo=read-the-docs&style=flat)](https://arntanguy.github.io/gram_savitzky_golay/)
![Package gram_savitzky_golay](https://github.com/arntanguy/gram_savitzky_golay/workflows/Package%20gram_savitzky_golay/badge.svg)

C++ Implementation of Savitzky-Golay filtering based on Gram polynomials, as described in
- [General Least-Squares Smoothing and Differentiation by the Convolution (Savitzky-Golay) Method](http://pubs.acs.org/doi/pdf/10.1021/ac00205a007)

## Installation

### From Ubunu packages

```sh
# Make sure you have required tools
sudo apt install apt-transport-https lsb-release ca-certificates gnupg
# Add our key
sudo apt-key adv --keyserver 'hkp://keyserver.ubuntu.com:80' --recv-key F6D3710D0B5016967A994DFFA650E12EFF6D3EDE
# Add our repository
echo "deb https://dl.bintray.com/arntanguy/ppa-head bionic main" | sudo tee -a /etc/apt/sources.list.d/arntanguy-head.list
# Install packages
sudo apt install libgram-savitzy-golay-dev
```

### From source

Dependencies:
- eigen3
- [SpaceVecAlg](https://github.com/jrl-umi3218/SpaceVecAlg) (optional): If present then you will be able to filter its main types (`PTransformd`, `MotionVecd`, etc)

```sh
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
sudo make install
```

How to include in your projects?
==
This package uses a modern cmake approach and exports its targets. To include in your own project, simply use:

```cmake
cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)

project(Example LANGUAGES CXX)

# Find the state-observation package and all its dependencies (Eigen)
find_package(gram_savitzky_golay REQUIRED)

# Creates a new executable and link it with the `gram_savitzky_golay` target
add_executable(Example example.cpp)
target_link_libraries(Example PUBLIC gram_savitzky_golay::gram_savitzky_golay)
# For SpaceVecAlg types
target_link_libraries(Example PUBLIC gram_savitzky_golay::gram_savitzky_golay_sva)
```

Example
==

```cpp
#include <gram_savitzky_golay/gram_savitzky_golay.h>

// Window size is 2*m+1
const size_t m = 3;
// Polynomial Order
const size_t n = 2;
// Initial Point Smoothing (ie evaluate polynomial at first point in the window)
// Points are defined in range [-m;m]
const size_t t = m;
// Derivation order? 0: no derivation, 1: first derivative, 2: second derivative...
const int d = 0;


// Real-time filter (filtering at latest data point)
gram_sg::SavitzkyGolayFilter filter(m, t, n, d);
// Filter some data
std::vector<double> data = {.1, .7, .9, .7, .8, .5, -.3};
double result = filter.filter(data);


// Real-time derivative filter (filtering at latest data point)
// Use first order derivative
// NOTE that the derivation timestep is assumed to be 1. If this is not the case,
// divide the filter result by the timestep to obtain the correctly scaled derivative
// See Issue #1
d=1;
gram_sg::SavitzkyGolayFilter first_derivative_filter(m, t, n, d);
// Filter some data
std::vector<double> values = {.1, .2, .3, .4, .5, .6, .7};
// Should be =.1
double derivative_result = first_derivative_filter.filter(values);
```


Filtering Rotations
==

```cpp
#include <gram_savitzky_golay/spatial_filters.h>

gram_sg::SavitzkyGolayFilterConfig sg_conf(50, 50, 2, 0);
gram_sg::RotationFilter filter(sg_conf);

filter.reset(Eigen::Matrix3d::Zero());

// Add rotation matrices to the filter
filter.add(Eigen::Matrix3d ...)
filter.add(Eigen::Matrix3d ...)
filter.add(Eigen::Matrix3d ...)

// Filter rotation matrices (weighted average of rotation matrices followed by an orthogonalization)
// See Peter Cork lecture here:
// https://www.cvl.isy.liu.se/education/graduate/geometry2010/lectures/Lecture7b.pdf
const Eigen::Matrix3d res = filter.filter();
```

This header also contains a filter for homogeneous transformations defined as `Eigen::Affine3d`, and a generic filter for eigen vectors.


Filtering Transformations
==

```cpp
#include <gram_savitzky_golay/spatial_filters.h>

gram_sg::SavitzkyGolayFilterConfig sg_conf(50, 50, 2, 0);
gram_sg::TransformFilter<Eigen::Affine3d> filter(sg_conf);  // Or Eigen::Matrix4d

filter.reset(Eigen::Matrix4d::Zero());

// Add homegenous transformation matrices to the filter
filter.add(Eigen::Matrix4d ...)
filter.add(Eigen::Matrix4d ...)
filter.add(Eigen::Matrix4d ...)

const Eigen::Matrix4d res = filter.filter();
```

If you have SpaceVecAlg installed, you can also filter `sva::PTransformd` in the same way.


Filtering Velocities
==

```cpp
#include <gram_savitzky_golay/spatial_filters.h>

gram_sg::SavitzkyGolayFilterConfig sg_conf(50, 50, 2, 0);
// Linear and angular velocity
gram_sg::VelocityFilter<Eigen::Vector6d> filter(sg_conf);

using Vector6d = Eigen::Matrix<double, 6, 1>;

// Add homegenous transformation matrices to the filter
filter.add(Vector6d{...})
filter.add(Vector6d{...})
...

const Eigen::Vector6d res = filter.filter();
```

Note that for standard vector types, this is the same as using `EigenVectorFilter` directly. However, you can also use more complex types such as `sva::MotionVecd`

```cpp
#include <gram_savitzky_golay/sva_filters.h>

gram_sg::SavitzkyGolayFilterConfig sg_conf(50, 50, 2, 0);
// Linear and angular velocity
gram_sg::VelocityFilter<sva::MotionVecd> filter(sg_conf);

// Add homegenous transformation matrices to the filter
filter.add(sva::MotionVecd{ 0.5, 1.2, 0} , {0.1, 0.2, 0.3}}};
filter.add(sva::MotionVecd{ 0.6, 1.3, 0} , {0.1, 0.2, 0.3}}};
filter.add(sva::MotionVecd{ 0.7, 1.4, 0} , {0.1, 0.2, 0.3}}};
...

sva::MotionVecd res = filter.filter();
```
