#!/bin/sh

# global build group variables
EDGE_CXX=mpiicpc
EDGE_PAR_COMPILE=16
EDGE_ARCH=skx
EDGE_PARALLEL=mpi+omp
EDGE_ELEMENT=tet4
EDGE_EQUATION=elastic
EDGE_BIN_DIR=./bin/

# configs to build, $order_$precision_$cfr
EDGE_CONFIGS="2_64_1 2_32_1 2_64_8 2_32_16 3_64_1 3_32_1 3_64_8 3_32_16 4_64_1 4_32_1 4_64_8 4_32_16 5_64_1 5_32_1 5_64_8 5_32_16 6_64_1 6_32_1 6_64_8 6_32_16 7_64_1 7_32_1 7_64_8 7_32_16"

# remove bin dir
rm -rf ${EDGE_BIN_DIR}
mkdir ${EDGE_BIN_DIR}

for c in ${EDGE_CONFIGS}
do
  # extract config
  EDGE_ORDER=`echo ${c} | awk -F"_" '{print $1}'`
  EDGE_PRECISION=`echo ${c} | awk -F"_" '{print $2}'`
  EDGE_CFR=`echo ${c} | awk -F"_" '{print $3}'`

  # cleanup
  rm -rf build/
  rm -rf .sconf_temp
  rm -rf .sconsign.dblite

  # finally let's build the code
  CXX=${EDGE_CXX} scons equations=${EDGE_EQUATION} order=${EDGE_ORDER} precision=${EDGE_PRECISION} cfr=${EDGE_CFR} element_type=${EDGE_ELEMENT} parallel=${EDGE_PARALLEL} arch=${EDGE_ARCH} xsmm=./libs zlib=./libs hdf5=./libs netcdf=./libs moab=./libs -j ${EDGE_PAR_COMPILE}

  # copy binary
  cp ./build/edge ${EDGE_BIN_DIR}/edge_${EDGE_ARCH}_${EDGE_PARALLEL}_${EDGE_ELEMENT}_${EDGE_EQUATION}_${EDGE_ORDER}_${EDGE_PRECISION}_${EDGE_CFR}
  ln -s edge_${EDGE_ARCH}_${EDGE_PARALLEL}_${EDGE_ELEMENT}_${EDGE_EQUATION}_${EDGE_ORDER}_${EDGE_PRECISION}_${EDGE_CFR} ${EDGE_BIN_DIR}/edge_${EDGE_ORDER}_${EDGE_PRECISION}_${EDGE_CFR}
done

