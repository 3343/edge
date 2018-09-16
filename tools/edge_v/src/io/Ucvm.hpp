/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
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
#ifndef EDGEV_IO_UCVM_HPP
#define EDGEV_IO_UCVM_HPP

#include <string>
#include <iostream>

extern "C" {
#include "ucvm.h"
}
#include "proj_api.h"

namespace edge_v {
  namespace io {
    class Ucvm;
  }
}

/**
 * @brief Interface to UCVM.
 */
class edge_v::io::Ucvm {
  public:
    /**
     * @brief Constructor, which initializes the UCVM interface.
     * 
     * @param i_config path to the UCVM config.
     * @param i_models UCVM models.
     * @param i_crdMode UCVM coordinate mode.
     */
    Ucvm( std::string const &i_config,
          std::string const &i_models,
          std::string const &i_crdMode ) {
      // coordinate mode
      ucvm_ctype_t l_crdMode = UCVM_COORD_GEO_ELEV;
      if( i_crdMode == "UCVM_COORD_GEO_DEPTH" )
        l_crdMode = UCVM_COORD_GEO_DEPTH;
  
      if( ucvm_init( i_config.c_str() ) != UCVM_CODE_SUCCESS ) {
        std::cerr << "UCVM init failed, aborting" << std::endl;
        exit(1);
      };

      if( ucvm_add_model_list( i_models.c_str() ) != UCVM_CODE_SUCCESS ) {
        std::cerr << "adding UCVM models failed, aborting: " << i_models << std::endl;
        exit(1);
      }

      if( ucvm_setparam( UCVM_PARAM_QUERY_MODE, l_crdMode ) != UCVM_CODE_SUCCESS  ) {
        std::cerr << "setting UCVM query mode failed, aborting: " << i_crdMode << std::endl;
        exit(1);
      }
    }

    /**
     * @brief Get the velocities from UCVM.
     * 
     * @param i_nPts number of points.
     * @param i_trafoSrc trafo, applied to the source coordinates before the proj.4 projection.
     * @param i_projSrc proj.4 string of the source projection.
     * @param i_projDes proh.4 string of the destination projection (used to query UCVM).
     * @param i_ucvmType ucvm type: gtl, crust or cmb.
     * @param o_vps will be set to p-wave velocities of the points.
     * @param o_vss will be set to s-wave velocities of the points.
     * @param o_rhos will be set to densities of the points.
     *
     * @paramt TL_T_ID integral type of the ids. 
     * @paramt TL_T_REAL floating point precision.
     */
    template< typename TL_T_ID,
              typename TL_T_REAL >
    void getVels( TL_T_ID             i_nPts,
                  double      const   i_trafoSrc[3][3],
                  std::string const  &i_projSrc,
                  std::string const  &i_projDes,
                  std::string const  &i_ucvmType,
                  double            (*i_pts)[3],
                  TL_T_REAL          *o_vps,
                  TL_T_REAL          *o_vss,
                  TL_T_REAL          *o_rhos ) {
      // allocate memory for UCVM points and data
      ucvm_point_t *l_ucvmPts  = new ucvm_point_t[ i_nPts ];
      ucvm_data_t  *l_ucvmData = new ucvm_data_t[  i_nPts ];

      // init UCVM points by applying the source trafo
      for( TL_T_ID l_pt = 0; l_pt < i_nPts; l_pt++ ) {
        double l_tmp[3] = {0,0,0};
        for( unsigned short l_d1 = 0; l_d1 < 3; l_d1++ )
          for( unsigned short l_d2 = 0; l_d2 < 3; l_d2++ )
            l_tmp[l_d1] += i_trafoSrc[l_d1][l_d2] * i_pts[l_pt][l_d2];

        for( unsigned short l_di = 0; l_di < 3; l_di++ )
          l_ucvmPts[l_pt].coord[l_di] = l_tmp[l_di];
      }

      // perform the Proj.4 projection
      projPJ l_projSrc = pj_init_plus( i_projSrc.c_str() );
      projPJ l_projDes = pj_init_plus( i_projDes.c_str() );

      // check the size of a UCVM point for the following query
      static_assert( sizeof(ucvm_point_t) == 3 * sizeof(double), "UCVM point size unexpected" );

      pj_transform(   l_projSrc,
                      l_projDes,
                      i_nPts,
                      3,
                    &(l_ucvmPts[0].coord[0]),
                    &(l_ucvmPts[0].coord[1]),
                    &(l_ucvmPts[0].coord[2]) );

      // apply rad to degree
      for( TL_T_ID l_pt = 0; l_pt < i_nPts; l_pt++ ) {
        l_ucvmPts[l_pt].coord[0] *= RAD_TO_DEG;
        l_ucvmPts[l_pt].coord[1] *= RAD_TO_DEG;
      }

      // perform the actual UCVM query
      int l_err  = ucvm_query( i_nPts,
                                      l_ucvmPts,
                                      l_ucvmData );
      if( l_err ) {
        std::cerr << "UCVM query failed, aborting" << std::endl;
        exit(1);
      }

      // convert from UCVM data type for float
      for( TL_T_ID l_pt = 0; l_pt < i_nPts; l_pt++ ) {
        // get pointer to UCVM properties
        ucvm_prop_t *l_propPtr;
        if( i_ucvmType == "gtl" )
          l_propPtr = &(l_ucvmData[l_pt].gtl);
        else if( i_ucvmType == "crust" )
          l_propPtr = &(l_ucvmData[l_pt].crust);
        else
          l_propPtr = &(l_ucvmData[l_pt].cmb);

        // store velocities
        o_vps[l_pt] = l_propPtr->vp;
        o_vss[l_pt] = l_propPtr->vs;
        o_rhos[l_pt] = l_propPtr->rho;
      }

      delete[] l_ucvmPts;
      delete[] l_ucvmData;
    }
};

#endif