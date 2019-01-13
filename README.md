# EDGƎ [![Build Status](https://travis-ci.org/3343/edge.svg?branch=master)](https://img.shields.io/coverity/scan/12907.svg)](https://scan.coverity.com/projects/3343-edge) [![codecov](https://codecov.io/gh/3343/edge/branch/master/graph/badge.svg)](https://codecov.io/gh/3343/edge) [![FOSSA Status](https://app.fossa.io/api/projects/git%2Bhttps%3A%2F%2Fgithub.com%2F3343%2Fedge.svg?type=shield)](https://app.fossa.io/projects/git%2Bhttps%3A%2F%2Fgithub.com%2F3343%2Fedge?ref=badge_shield)

## Equations
* *advection (1D)*:

  Solves the one-dimensional advection equation. ```q(x,t)``` is a scalar. The scalar advection speed ```a(x)``` can be set per element, but has to be either positive or negative for the entire domain.

  ```
  q_t + a * q_x = 0
  ```

* *advection (2D)*:

  Solves the two-dimensional advection equation. ```q(x,y,t)``` is a scalar. The scalar advection speeds ```a(x,y)``` and ```b(x,y)``` can be set per element. Each has to be either positive or negative for the entire domain.

  ```
  q_t + a * q_x + b * q_y = 0
  ```

* *advection (3D)*:

  Solves the three-dimensional advection equation. ```q(x,y,z,t)``` is a scalar. The scalar advection speeds ```a(x,y,z)```, ```b(x,y,z)``` and ```c(x,y,z)``` can be set per element. Each has to be either positive or negative for the entire doman.

  ```
  q_t + a * q_x + b * q_y + c q_z = 0
  ```

* *elastics (2D)*:

  Solves the two-dimensional elastic wave equations. The vector of quantities ```q(x,y,t)=(sigma_xx, sigma_yy, sigma_xy, u, v )``` contains the normal stress components ```sigma_xx``` and ```sigma_yy```, the shear stress ```sigma_xy``` and the two particle velocities ```u``` and ```v``` in ```x-``` and ```y-```direction respectively. The Jacobians ```A(x,y)``` and ```B(x,y)``` are allowed to be set per element and summarize the material parameters.

  ```
  q_t + A q_x + B q_y = 0
  ```
  
* *elastics (3D)*:

  Solves the three-dimensional elastic wave equations. The vector of quantities ```q(x,y,z,t)=(sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz, u, v, w )``` contains the normal stress components ```sigma_xx```, ```sigma_yy``` and ```sigma_zz```, the shear stresses ```sigma_xy```, ```sigma_xz``` and ```sigma_yz``` and the three particle velocities ```u```, ```v``` ```w```  in ```x-```, ```y-``` and ```z-```direction respectively. The Jacobians ```A(x,y,z)```, ```B(x,y,z)``` and ```C(x,y,z)``` are allowed to be set per element and summarize the material parameters.

  ```
  q_t + A q_x + B q_y + C q_z = 0
  ```

* *swe (1D)*:

  Solves the one-dimensional Shallow Water Equations (SWE) in conservative form. The conserved quantities ```q(x,t)=(h,hu)``` are the water height ```h``` and the momentum ```hu```. The flux function is nonlinear. Bathymetry is supported.

  ```
  q_t + f(q)_x = 0,

         |         hu           |
  f(q) = |                      |
         | hu^2 + 1/2 * g * h^2 |
  ```

* *swe (2D)*:

  Solves the two-dimensional Shallow Water Equations (SWE) in conservative form. The conserved quantities `q(x,t)=(h,hu,hv)` are the water height `h`, the momentum `hu` in x-direction and the momentum `hv` in y-direction. The flux function is nonlinear. Bathymetry is supported.

  ```
  q_t + f(q)_x + g(q)_y = 0,

         |         hu           |         |          hv          |
         |                      |         |                      |
  f(q) = | hu^2 + 1/2 * g * h^2 |, g(q) = |          huv         |
         |                      |         |                      |
         |         huv          |         | hv^2 + 1/2 * g * h^2 |
  ```

## Elements
* *line (1D)*:

  Line element. Element width ```dx``` is allowed to change in every element.

* *quad4r (2D)*:

  Rectangular, 4-node quadrilaterals. Widths ```dx``` and ```dy``` are allowed to change on a per-row/per-column basis (conforming mesh).

* *tria3 (2D)*:

  3-node triangles. Arbitrary, conforming triangulations of the computational domain are supported.

* *hex8r (3D)*:

  Rectangular, 8-node hexahedrons (bricks). Widths ```dx```, ```dy``` and ```dz``` are allowed to change on a conforming mesh basis.
  
* *tet4 (3D)*:

  4-node tetrahedrons. Arbitrary, conforming tetrahedralization are allowed.

## Feature table

Based on the equations and the element type, the following table shows the implemented features:

| equations |         element types            | CFR | FV | ADER-DG | LIBXSMM |
|-----------|:--------------------------------:|:---:|:--:|:-------:|:-------:|
| advection | line, quad4r, tria3, hex8r, tet4 |  x  |  x |    x    |         |
| elastics  | quad4r, tria3, hex8r, tet4       |  x  |  x |    x    |    x    |
| swe       | line, quad4r, tria3              |  x  |  x |         |         |

## High Performance Support
| Microarchitecture | Machine(s) |
|-------------------|:----------:|
| [Sandy Bridge](https://ark.intel.com/products/codename/64276/Sandy-Bridge-EP) | [Stampede 1](https://portal.tacc.utexas.edu/user-guides/stampede) |
| [Bulldozer](http://products.amd.com/en-us/search/cpu/amd-opteron%E2%84%A2/amd-opteron%E2%84%A2-6200-series-processor) | [Blue Waters](https://bluewaters.ncsa.illinois.edu/hardware-summary) |
| [Haswell](https://ark.intel.com/products/codename/42174/Haswell) | [Comet](http://www.sdsc.edu/support/user_guides/comet.html), [Cori Phase 1](http://www.nersc.gov/users/computational-systems/cori/configuration/) |
| [Knights Landing](https://ark.intel.com/products/codename/48999/Knights-Landing)    | [Stampede 2](https://portal.tacc.utexas.edu/user-guides/stampede2), [Cori Phase 2](http://www.nersc.gov/users/computational-systems/cori/configuration/), [Theta](https://www.alcf.anl.gov/theta) |
| [Skylake](https://ark.intel.com/products/codename/37572/Skylake) | [Amazon Elastic Compute Cloud](https://aws.amazon.com/ec2/), [Google Cloud Platform](https://cloud.google.com/), [Stampede 2](https://portal.tacc.utexas.edu/user-guides/stampede2) |
| [EPYC](http://products.amd.com/en-us/search/cpu/amd-epyc%E2%84%A2)                                         | [packet](https://www.packet.com/) |
| [Knights Mill](https://ark.intel.com/products/codename/57723/Knights-Mill) | - |
