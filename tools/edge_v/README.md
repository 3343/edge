# EDGƎ-V

[![License](https://img.shields.io/badge/license-BSD3-blue.svg)](LICENSE.md) 
<!--[![Travis CI](https://travis-ci.org/hfp/libxsmm.svg?branch=master "Master branch build status")](https://github.com/hfp/libxsmm/wiki/Status) -->

EDGE Velocity (EDGE-V) is a tool used to annotate the meshes with velocity model data. It is developed based on [Unified Community Velocity Model](http://scec.usc.edu/scecpedia/UCVMC)(UCVM), which is the data source, and [Mesh-Oriented dAtABase](http://sigma.mcs.anl.gov/moab-library)(MOAB).

EDGE-V is released as an open-source scientific software under the BSD3 software license.


## System and Software Requirements

Edge-V works with the following operating sysetms and software stacks.

*  CentOS 7 Linux x86_64-linux 
*  GNU gcc compilers version 4.8.5
*  MPI c/c++ compilers - openmpi 1.8.8 or mpich 1.2.6 or intel MPI compilers
*  Autotools build software for Linux
*  Unified Community Velocity Model C-language (UCVMC) library: http://scec.usc.edu/scecpedia/UCVMC/ 
*  Proj.4 projection library: http://trac.osgeo.org/proj/ (provided in UCVMC)
*  Mesh-Oriented dAtABase (MOAB) library: http://sigma.mcs.anl.gov/moab-library/ (contained in submodules of EDGE)


## Build Instructions

### UCVMC

To set up the UCVM C-interface library, please follow the commands below:

```bash
git clone https://github.com/SCECcode/UCVMC.git
cd UCVMC/largefiles
./get_large_files.py
./check_largefiles_md5.py
./stage_large_files.py
cd ..
./ucvm_setup.py
```

During the installation process in `./ucvm_setup.py`, it is recommended to install all the libraries and models to get the full support. Be aware that it will occupy a large amount of storage.
Please refer to [UCVMC repo](https://github.com/SCECcode/UCVMC#ucvmc) for detailed intructions.

### MOAB

Please refer to [EDGE installation guide](https://usr.dial3343.org/chapters/install/edge.html) for compatible setup.

### Edge-V

Finally, to build Edge-V tool, please run the following command:

```bash
PREFIX=$(pwd) make MOABDIR="path_to_MOAB" UCVMDIR="path_to_UCVMC" PROJ4DIR="path_to_Proj_4" 
```
The paths to MOAB library and UCVMC library are required. If `PROJ4DIR` is not provided, the Proj.4 library within UCVMC will be searched and used. By default, the released tool is in `$(pwd)`. Please set `PREFIX` to change it.

One can also set up the dependent libraries in the `Makefile.inc` file, and simply run:
```bash
make
```

**NOTE**: The default c++ compiler to build Edge-V is the MPI compiler because MOAB and its dependency are usually built as MPI version. Be careful when switching to GNU c++ compiler.



## Usage

### Overview

To use the Edge-V tool for annotation, simply run:

```bash
./bin/edge_v -f config_file.log
```

Here `config_file.log` is the path to the annotation configuration file. There is an example configuration in `./example/annotation.conf`:

```bash
# initialization params for UCVMC
ucvm_config=[@path_to_UCVMC]/conf/ucvm.conf
ucvm_model_list=cvmsi

# input mesh file path
mesh_file=./meshes/ucvm_mini.msh

# output file path
node_vm_file=./meshes/ucvm_mini_node.vel
elmt_vm_file=./meshes/ucvm_mini_elmt.vel
h5m_file=./meshes/ucvm_mini_vmtags.h5m
```

All the settings are required. 
* `ucvm_config` : 
    The configuration file for UCVM. There is a reference config file at `$(UCVMDIR)/conf/ucvm.conf`.
    (correspondent to the argument of `-f` option in `ucvm_query`)

* `ucvm_model_list` :
    The sub-model used to provide velocity model data.
    (correspondent to the argument of `-m` option in `ucvm_query`)

**Note**: For details in the UCVMC initial parameters, please refer to [UCVMC wiki](https://github.com/SCECcode/UCVMC/wiki).

* `mesh_file` :
    The input mesh file in Gmsh’s native “MSH” ASCII file format.
* `node_vm_file` :
    The output velocity model data for all the nodes in the mesh in ASCII file format.
    The velocity model data is parameterized as `lambda`, `mu` and `rho`.

* `elmt_vm_file` :
    The output velocity model data for all the tets in the mesh in ASCII file format.
    The velocity model data is parameterized as `lambda`, `mu` and `rho`.

* `h5m_file` :
    The output mesh file with velocity model annotation in H5M file format.
    The velocity model data is parameterized as `lambda`, `mu` and `rho`, recorded as dense, 4-byte tags for only tetrahedron elements in the mesh. Each 4-byte data is interpreted as a single-precision floating-point number for use.

### Example

A basic example is provided with `example/annotation.conf` and `meshes/ucvm_mini.msh`. Change the `ucvm_config` setting to the path of the reference config file in UCVMC (could be found at ${UCVMC_DIR}/conf/ucvm.conf), and run the following command:

```bash
~bash $ ./bin/edge_v -f example/annotation.conf 
Reading Annotation Config File: example/annotation.conf ... Done!
Reading Mesh File: ./meshes/ucvm_mini.msh ... 
 | Number of vertices is 13789
 | Number of elements is 58232
Done!
UCVM Query ... Done!
Write Velocity Model: ./meshes/ucvm_mini_node.vel ... Done!
Writing Velocity Model: ./meshes/ucvm_mini_elmt.vel ... Done!
Writing Annotated Mesh File: ./meshes/ucvm_mini_vmtags.h5m ... Done!
```

It should output the same logging info as above. The velocity model files and the annotation file are generated in `meshes/`. It is also recommend to convert `.h5m` file to `.vtk` for visualization.

```bash
$(MOAD_DIR)/bin/mbconvert -f VTK ./meshes/ucvm_mini_vmtags.h5m ./meshes/ucvm_mini_vmtags.vtk
```
