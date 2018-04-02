# EDGƎ-V

[![License](https://img.shields.io/badge/license-BSD3-blue.svg)](LICENSE.md) 


EDGE Velocity (EDGE-V) is a tool used to annotate the meshes with velocity model data. It is developed based on [Unified Community Velocity Model](http://scec.usc.edu/scecpedia/UCVMC)(UCVM), which is the data source, and [Mesh-Oriented datABase](http://sigma.mcs.anl.gov/moab-library)(MOAB).

EDGE-V is released as an open-source scientific software under the BSD3 software license.


## System and Software Requirements

Edge-V works with the following operating systems and software stacks.

*  CentOS 7 Linux x86_64-linux 
*  GNU gcc compilers version 4.8.5
*  MPI C/C++ compilers - openmpi 1.8.8 or mpich 1.2.6 or Intel MPI compilers
*  Autotools build software for Linux
*  Unified Community Velocity Model C-language (UCVMC) library: https://scec.usc.edu/scecpedia/UCVMC
*  Proj.4 projection library: http://trac.osgeo.org/proj/ (provided in UCVMC)
*  Mesh-Oriented datABase (MOAB) library: http://sigma.mcs.anl.gov/moab-library/ (contained in submodules of EDGE)


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

**NOTE**: The default C++ compiler to build Edge-V is the MPI compiler because MOAB and its dependency are usually built as MPI version. Be careful when switching to GNU C++ compiler.



## Usage

### Overview

To use the Edge-V tool for annotation, simply run:

```bash
./bin/edge_v -f config_file.log
```

Here `config_file.log` is the path to the annotation configuration file. There is an example configuration in `./example/annotation.conf`:

```bash
# initialization params for UCVMC
ucvm_config=./conf/ucvm.conf
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

A basic example is provided with `example/annotation.conf`. First, download the mini mesh:

```bash
mkdir -p meshes
wget https://bitbucket.org/3343/edge_opt/raw/HEAD/tools/edge_v/meshes/ucvm_mini.msh -O ./meshes/ucvm_mini.msh
```

Change the `ucvm_config` setting in `example/annotation.conf` to the path of the reference config file in UCVMC (could be found at `${UCVMC_DIR}/conf/ucvm.conf` ), and run the following command:

```bash
~bash $ ./bin/edge_v -f example/annotation.conf 
Reading Annotation Config File: example/annotation.conf... Done!
Reading Mesh File: ./meshes/ucvm_mini.msh... 
 | Number of vertices is 13789
 | Number of elements is 58232
Done!
UCVM Query... Done!
Write Velocity Model: ./meshes/ucvm_mini_node.vel... Done!
No fault input files... Skipping fault annotation
Writing Velocity Model: ./meshes/ucvm_mini_elmt.vel... Done!
Writing Annotated Mesh File: ./meshes/ucvm_mini_vmtags.h5m... Done!
```

This should output the same logging info as above. The velocity model files and the annotation file are generated in `meshes/`.


## Velocity Model based Mesh Refinement

<img style="float: right;" src="https://scec.usc.edu/scecwiki/images/thumb/6/62/Base_lahabra-win1.png/250px-Base_lahabra-win1.png">

Here is an example of a velocity based mesh refinement process for the [La Habra simulation region](https://scec.usc.edu/scecpedia/La_Habra_Simulation_Region). A script is provided (`mesh_refinement.sh`) to perform iterative mesh refinement in order to generate a final refined mesh that is based on the velocity model (CVM-S4.26.M01) queried from UCVM.

First, we need to download the geo file:
```bash
mkdir -p meshes
wget https://bitbucket.org/3343/edge_opt/raw/HEAD/tools/edge_v/meshes/la_habra_small.geo -O ./meshes/la_habra_small.geo
```

The script uses this `la_habra_small.geo` file to produce the first mesh file (without any optimization), uses this mesh file to run `write_pos` program which generates a background velocity map to be used in the next meshing iteration. The subsequent mesh iteration uses this background pos file and generates a refined mesh.

Open `mesh_refinement.sh` in your favorite text-editor and replace the `PLACEHOLDER` tokens with user-specific values, like:
```
declare -g model='la_habra_small'
confDir='./example/'
mshDir='./meshes/'
writePos='./bin/write_pos'
```

There is a provision in the script to utilize any high clock speed CPU rack (that you have access) to perform the heavy-duty meshing operations and send intermediate files back and forth automatically between the local machine and the remote client. In order to use this provision, the user needs to set the script variable `remoteMsh` to 1 (by default 0). **Do note that this will require you to have a working SSH public key (generated using ssh-keygen) established with the remote client and also gmsh to have been installed on the remote client somewhere (we will need this location)**.
```
declare -g remoteMsh=1
```

If you are using the above-mentioned remote meshing provision, you also need to set your ssh username, remote-client domain, gmsh installation location on the remote client and the mesh directory on the remote client which will store the intermediate files. For example, (change accordingly)
```
declare -g remoteUsr='myuser'
declare -g remoteClient='client.domain.name'
declare -g remoteGmshDir='/home/myuser/path-to/gmsh-3.0.6-Linux64/bin/gmsh'
declare -g remoteMshDir='/home/myuser/path-to/la_habra_small/'
```

Change the `ucvm_config` setting in `example/la_habra_small.conf` to the path of the reference config file in UCVMC (could be found at `${UCVMC_DIR}/conf/ucvm.conf` ), and run the following command:
```bash
./mesh_refinement.sh 2>&1 | tee mesh_refinement.sh.log
```

At the end of the iterations, a final refined mesh file is generated in `meshes/` and this file can now be passed on to `edge_v` (see above sections for instructions on how to use the `edge_v` program) in order to generate the velocity information annotated h5m file of the La Habra simulation region.

This process of mesh refinement can be repeated for any other model for which you have the geo file. Simply use that model-name (name used for geo file) in the `mesh_refinement.sh` script and make sure your other paths are correct.
