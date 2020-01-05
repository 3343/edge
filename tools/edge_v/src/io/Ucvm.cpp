/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020, Alexander Breuer
 * Copyright (c) 2018, Regents of the University of California
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Interface to UCVM.
 **/
#include "Ucvm.h"
extern "C" {
#include <ucvm.h>
extern int ucvm_init_flag; // used as workaround for multiple objects of this class
}
#include <proj_api.h>
#include "io/logging.h"

edge_v::io::Ucvm::Ucvm( std::string const & i_config,
                        std::string const & i_models,
                        std::string const & i_crdMode ) {
  // set coordinate mode
  ucvm_ctype_t l_crdMode = UCVM_COORD_GEO_ELEV;
  if( i_crdMode == "UCVM_COORD_GEO_DEPTH" ) {
    l_crdMode = UCVM_COORD_GEO_DEPTH;
  }
  else if( i_crdMode == "UCVM_COORD_GEO_ELEV" ) {
    l_crdMode = UCVM_COORD_GEO_ELEV;
  }
  else EDGE_V_LOG_FATAL;

  // only init UCVM once due to issues with finalize
  if( !ucvm_init_flag ) {
    // init UCVM, add models and set parameters
    int l_err = ucvm_init( i_config.c_str() );
    EDGE_V_CHECK_EQ( l_err, UCVM_CODE_SUCCESS );

    l_err = ucvm_add_model_list( i_models.c_str() );
    EDGE_V_CHECK_EQ( l_err, UCVM_CODE_SUCCESS );

    l_err = ucvm_setparam( UCVM_PARAM_QUERY_MODE, l_crdMode );
    EDGE_V_CHECK_EQ( l_err, UCVM_CODE_SUCCESS );
  }
}

edge_v::io::Ucvm::~Ucvm() {
  // intentionally no ucvm_finalize due to memory errors of the lib
}

void edge_v::io::Ucvm::getVels( std::size_t          i_nPts,
                                double      const    i_trafoSrc[3][3],
                                std::string const  & i_projSrc,
                                std::string const  & i_projDes,
                                std::string const  & i_ucvmType,
                                double            (* i_pts)[3],
                                float              * o_vps,
                                float              * o_vss,
                                float              * o_rhos ) {
  // allocate memory for UCVM points and data
  ucvm_point_t *l_ucvmPts  = new ucvm_point_t[ i_nPts ];
  ucvm_data_t  *l_ucvmData = new ucvm_data_t[  i_nPts ];

  // init UCVM points by applying the source trafo
  for( std::size_t l_pt = 0; l_pt < i_nPts; l_pt++ ) {
    double l_tmp[3] = {0,0,0};
    for( unsigned short l_d1 = 0; l_d1 < 3; l_d1++ )
      for( unsigned short l_d2 = 0; l_d2 < 3; l_d2++ )
        l_tmp[l_d1] += i_trafoSrc[l_d1][l_d2] * i_pts[l_pt][l_d2];

    for( unsigned short l_di = 0; l_di < 3; l_di++ )
      l_ucvmPts[l_pt].coord[l_di] = l_tmp[l_di];
  }

  // init the proj.4 projections
  projPJ l_projSrc = pj_init_plus( i_projSrc.c_str() );
  projPJ l_projDes = pj_init_plus( i_projDes.c_str() );

  // check the size of a UCVM point for the following query
  static_assert( sizeof(ucvm_point_t) == 3 * sizeof(double), "UCVM point size unexpected" );

  // perform projections and convert to degree
  pj_transform(   l_projSrc,
                  l_projDes,
                  i_nPts,
                  3,
                &(l_ucvmPts[0].coord[0]),
                &(l_ucvmPts[0].coord[1]),
                  NULL );

  for( std::size_t l_pt = 0; l_pt < i_nPts; l_pt++ ) {
    l_ucvmPts[l_pt].coord[0] *= RAD_TO_DEG;
    l_ucvmPts[l_pt].coord[1] *= RAD_TO_DEG;
    l_ucvmPts[l_pt].coord[2] = i_pts[l_pt][2];
  }

  // perform the actual UCVM query
  int l_err  = ucvm_query( i_nPts,
                           l_ucvmPts,
                           l_ucvmData );
  EDGE_V_CHECK_EQ( l_err, 0 );

  // convert from UCVM data type to float
  for( std::size_t l_pt = 0; l_pt < i_nPts; l_pt++ ) {
    // get pointer to UCVM properties
    ucvm_prop_t *l_propPtr;
    if( i_ucvmType == "gtl" ) {
      l_propPtr = &(l_ucvmData[l_pt].gtl);
    }
    else if( i_ucvmType == "crust" ) {
      l_propPtr = &(l_ucvmData[l_pt].crust);
    }
    else {
      l_propPtr = &(l_ucvmData[l_pt].cmb);
    }

    // store velocities and density
    o_vps[l_pt]  = l_propPtr->vp;
    o_vss[l_pt]  = l_propPtr->vs;
    o_rhos[l_pt] = l_propPtr->rho;
  }

  // free memory
  delete[] l_ucvmPts;
  delete[] l_ucvmData;
}