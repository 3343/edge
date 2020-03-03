/**
 * @file This file is part of EDGE.
 *
 * @author Rajdeep Konwar (rkonwar AT ucsd.edu)
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
 * Copyright (c) 2017-2018, Regents of the University of California
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
 * Entry point of standalone EDGE-V.
 **/
#include "io/Config.h"
#include "io/Csv.h"
#include "io/Hdf5.h"
#include "io/BgMeshMsh4.h"
#include "models/Constant.h"
#include "models/seismic/Expression.h"
#ifdef PP_HAS_UCVM
#include "models/seismic/Ucvm.h"
#endif
#include "mesh/Refinement.h"
#include "mesh/Mesh.h"
#include "time/Cfl.h"
#include "time/Groups.h"
#include "mesh/Partition.h"
#include "mesh/Communication.h"

#include "io/logging.h"
#ifdef PP_USE_EASYLOGGING
INITIALIZE_EASYLOGGINGPP
#endif

int main( int i_argc, char *i_argv[] ) {
#ifdef PP_USE_EASYLOGGING
  START_EASYLOGGINGPP( i_argc, i_argv );
#endif

  EDGE_V_LOG_INFO << "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVVV   VVVVVVVVVVVVVV            VVVVVVVVVVVVVVV  VVVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVVV   VVVVVVVVVVVVVVV         VVVVVVVVVVVVVVVV   VVVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVV            VVVVV       VVVVV      VVVVVV                       VVVVV";
  EDGE_V_LOG_INFO << "VVVVV            VVVVV        VVVVV    VVVVV                         VVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVV    VVVVV         VVVVV  VVVVV                  VVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVV    VVVVV         VVVVV  VVVVV      VVVVVVVVV   VVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVV            VVVVV         VVVVV  VVVVV      VVVVVVVVV           VVVVV";
  EDGE_V_LOG_INFO << "VVVVV            VVVVV        VVVVV    VVVVV        VVVVVV           VVVVV";
  EDGE_V_LOG_INFO << "VVVVV            VVVVV       VVVVV      VVVVV       VVVVV            VVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVVVV  VVVVVVVVVVVVVVV         VVVVVVVVVVVVVVV   VVVVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVVVV  VVVVVVVVVVVVVV           VVVVVVVVVVVVV    VVVVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV";
  EDGE_V_LOG_INFO << "";
  EDGE_V_LOG_INFO << "              EDGE-V is available from: https://dial3343.org";
  EDGE_V_LOG_INFO << "";

  // parse argumenets
  std::vector< std::string > l_args( i_argv, i_argv + i_argc );
  if( l_args.size() != 3 || (l_args[1] != "-x" && l_args[1] != "--xml") ) {
    EDGE_V_LOG_INFO << "USAGE: " << l_args[0] << " [options]";
    EDGE_V_LOG_INFO;
    EDGE_V_LOG_INFO << "  Options:";
    EDGE_V_LOG_INFO << "   --help, -h Print usage and exit.";
    EDGE_V_LOG_INFO << "   --xml,  -x Location of the XML configuration.";

    for( std::size_t l_ar = 0; l_ar < l_args.size(); l_ar++ ) {
      if( l_args[l_ar] == "-h" || l_args[l_ar] == "--help" ) {
        return EXIT_SUCCESS;
      }
    }

    return EXIT_FAILURE;
  }

  EDGE_V_LOG_INFO << "EDGE version: " << PP_EDGE_VERSION;
  EDGE_V_LOG_INFO << "entity size in bytes: " << sizeof(edge_v::t_idx);
  EDGE_V_LOG_INFO << "parsing xml config";
  std::string l_xml = l_args[2];
  edge_v::io::Config l_config( l_xml );
  EDGE_V_LOG_INFO << "sharing runtime config:";
  EDGE_V_LOG_INFO << "  mesh:";
  EDGE_V_LOG_INFO << "    periodic:                  " << l_config.getPeriodic();
  EDGE_V_LOG_INFO << "    write_element_annotations: " << l_config.getWriteElAn();
  EDGE_V_LOG_INFO << "    n_partitions:              " << l_config.nPartitions();
  EDGE_V_LOG_INFO << "    in:                        " << l_config.getMeshIn();
  EDGE_V_LOG_INFO << "    out:";
  EDGE_V_LOG_INFO << "      base:                    " << l_config.getMeshOutBase();
  EDGE_V_LOG_INFO << "      extension:               " << l_config.getMeshOutExt();
  EDGE_V_LOG_INFO << "  velocity model:";
#ifdef PP_HAS_UCVM
  if( l_config.getVelModUcvmProjSrc() != "" ) {
    EDGE_V_LOG_INFO << "    ucvm:";
    EDGE_V_LOG_INFO << "      source_transformation:";
    EDGE_V_LOG_INFO << "        " << l_config.getVelModUcvmTrafoSrc()[0][0]
                    << " "        << l_config.getVelModUcvmTrafoSrc()[0][1]
                    << " "        << l_config.getVelModUcvmTrafoSrc()[0][2];
    EDGE_V_LOG_INFO << "        " << l_config.getVelModUcvmTrafoSrc()[1][0]
                    << " "        << l_config.getVelModUcvmTrafoSrc()[1][1]
                    << " "        << l_config.getVelModUcvmTrafoSrc()[1][2];
    EDGE_V_LOG_INFO << "        " << l_config.getVelModUcvmTrafoSrc()[2][0]
                    << " "        << l_config.getVelModUcvmTrafoSrc()[2][1]
                    << " "        << l_config.getVelModUcvmTrafoSrc()[2][2];
    EDGE_V_LOG_INFO << "      projections:";
    EDGE_V_LOG_INFO << "        source: " << l_config.getVelModUcvmProjSrc();
    EDGE_V_LOG_INFO << "        destination: " << l_config.getVelModUcvmProjDes();
    EDGE_V_LOG_INFO << "      models: " << l_config.getVelModUcvmModels();
    EDGE_V_LOG_INFO << "      model_type: " << l_config.getVelModUcvmModelType();
    EDGE_V_LOG_INFO << "      coordinate_mode: " << l_config.getVelModUcvmCrdMode();
    EDGE_V_LOG_INFO << "      normalization_rule: " << l_config.getVelModUcvmRule();
    EDGE_V_LOG_INFO << "      lower_to_surface: " << l_config.getLowerToSurf();
  }
#else
  if( false ){}
#endif
  else if( l_config.getVelModSeismicExpr() != "" ) {
    EDGE_V_LOG_INFO << "    seismic expression: ";
    EDGE_V_LOG_INFO << "      " << l_config.getVelModSeismicExpr();
  }
  else {
    EDGE_V_LOG_INFO << "  constant";
  }
  if( l_config.getRefExpr() != "" ) {
    EDGE_V_LOG_INFO << "  refinement:";
    EDGE_V_LOG_INFO << "    expression: " << l_config.getRefExpr();
    EDGE_V_LOG_INFO << "    out: " << l_config.getRefOut();
  }
  EDGE_V_LOG_INFO << "  time:";
  EDGE_V_LOG_INFO << "    #time groups: " << l_config.nTsGroups();
  EDGE_V_LOG_INFO << "    fun dt:       " << l_config.getFunDt();
  EDGE_V_LOG_INFO << "    out:          " << l_config.getTsOut();

  EDGE_V_LOG_INFO << "initializing MOAB";
  std::string l_meshIn = l_config.getMeshIn();
  edge_v::io::Moab l_moab( l_meshIn );
  l_moab.printStats();

  // initialize and set mesh data
  EDGE_V_LOG_INFO << "initializing mesh interface";
  if( l_config.getPeriodic() ) EDGE_V_LOG_INFO << "  assuming periodioc boundaries";
  edge_v::mesh::Mesh* l_mesh = new edge_v::mesh::Mesh( l_moab,
                                                       l_config.getPeriodic() );
  l_mesh->printStats();

  EDGE_V_LOG_INFO << "initializing velocity model";
  edge_v::models::Model *l_velMod = nullptr;
#ifdef PP_HAS_UCVM
  std::string l_tagVp = "edge_v_vp";
  std::string l_tagVs = "edge_v_vs";
  std::string l_tagRho = "edge_v_rho";

  edge_v::io::Ucvm l_ucvmReader( PP_UCVM_CONF,
                                 l_config.getVelModUcvmModels(),
                                 l_config.getVelModUcvmCrdMode(),
                                 l_config.getVelModUcvmRule(),
                                 l_config.getLowerToSurf() );
  if( l_config.getVelModUcvmProjSrc() != "" ) {
    l_velMod = new edge_v::models::seismic::Ucvm( l_ucvmReader,
                                                  l_config.getVelModUcvmTrafoSrc(),
                                                  l_config.getVelModUcvmProjSrc(),
                                                  l_config.getVelModUcvmProjDes(),
                                                  l_config.getVelModUcvmModelType() );
  }
#else
  if( false ){}
#endif
  else if( l_config.getVelModSeismicExpr() != ""  ) {
    l_velMod = new edge_v::models::seismic::Expression( l_config.getVelModSeismicExpr() );
  }
  else {
    l_velMod = new edge_v::models::Constant( 1 );
  }

  if( l_config.getRefExpr() != "" ) {
    EDGE_V_LOG_INFO << "computing mesh refinement";
    l_velMod->init( l_mesh->nVes(),
                    l_mesh->getVeCrds() );

    edge_v::mesh::Refinement l_ref;
    l_ref.init(  l_mesh->nVes(),
                 l_mesh->nEls(),
                 l_mesh->nElVes(),
                 l_mesh->getElVe(),
                 l_mesh->getVeCrds(),
                 l_config.getRefExpr(),
                *l_velMod );

    EDGE_V_LOG_INFO << "writing mesh refinement: " << l_config.getRefOut();
    std::ofstream l_stream( l_config.getRefOut(),
                            std::ios::out );
    edge_v::io::BgMeshMsh4::write( l_mesh->getTypeEl(),
                                   l_mesh->nVes(),
                                   l_mesh->nEls(),
                                   l_mesh->getElVe(),
                                   l_mesh->getVeCrds(),
                                   l_ref.getTargetLengthsVe(),
                                   l_stream );
    l_stream.close();

    if( l_config.getWriteElAn() ) {
      EDGE_V_LOG_INFO << "storing elements' target lengths";
      l_moab.setEnData( l_mesh->getTypeEl(),
                        "edge_v_target_lengths",
                        l_ref.getTargetLengthsEl() );
    }
  }

  // abort if mesh doesn't have output
  if( l_config.getMeshOutBase() == "" ) {
    delete l_mesh;
    delete l_velMod;
    return EXIT_SUCCESS;
  }

  EDGE_V_LOG_INFO << "computing CFL time steps";
  edge_v::time::Cfl *l_cfl = new edge_v::time::Cfl(  l_mesh->getTypeEl(),
                                                     l_mesh->nVes(),
                                                     l_mesh->nEls(),
                                                     l_mesh->getElVe(),
                                                     l_mesh->getVeCrds(),
                                                     l_mesh->getInDiasEl(),
                                                    *l_velMod );
  l_cfl->printStats();

  if( l_config.getWriteElAn() ) {
    EDGE_V_LOG_INFO << "storing elements' cfl time steps";
    std::string l_tagCfl = "edge_v_cfl_time_steps";
    l_moab.setEnData( l_mesh->getTypeEl(),
                      l_tagCfl,
                      l_cfl->getTimeSteps() );
  }

  if( l_config.getTsOut() != "" ) {
    EDGE_V_LOG_INFO << "writing cfl time steps to ascii";
    std::string l_colNames = "ts_cfl";
    double const * l_data[1] = { l_cfl->getTimeSteps() };
    edge_v::io::Csv::write( l_config.getTsOut(),
                            1,
                            l_mesh->nEls(),
                            &l_colNames,
                            5,
                            l_data );
  }

  EDGE_V_LOG_INFO << "computing time step groups";
  unsigned short l_nRates = l_config.nTsGroups()-1;
  double *l_rates = new double[l_nRates];
  for( unsigned short l_tg = 0; l_tg < l_nRates; l_tg++ ) {
    l_rates[l_tg] = 2.0;
  }

  double l_funDt = l_config.getFunDt();
  // search for fundamental dt, if not specified
  double l_speedUp = 0;
  if( l_funDt == 0 ) {
    double l_dt = 1.0;
    while( l_dt > 0.5 ) {
      edge_v::time::Groups l_tsGroups( l_mesh->getTypeEl(),
                                       l_mesh->nEls(),
                                       l_mesh->getElFaEl(),
                                       l_nRates,
                                       l_rates,
                                       l_dt,
                                       l_cfl->getTimeSteps() );
      if( l_tsGroups.getSpeedUp() > l_speedUp ) {
        l_speedUp = l_tsGroups.getSpeedUp();
        l_funDt = l_dt;
        EDGE_V_LOG_INFO << "  found new best fundamental dt / speedup: " << l_funDt << " / " << l_speedUp;
      }
      l_dt -= 0.01;
    }
  }

  edge_v::time::Groups *l_tsGroups = new edge_v::time::Groups(  l_mesh->getTypeEl(),
                                                                l_mesh->nEls(),
                                                                l_mesh->getElFaEl(),
                                                                l_nRates,
                                                                l_rates,
                                                                l_funDt,
                                                                l_cfl->getTimeSteps() );
  l_tsGroups->printStats();

  if( l_config.getWriteElAn() ) {
    EDGE_V_LOG_INFO << "storing elements' time step groups";
    std::string l_tagElTg = "edge_v_element_time_groups";
    l_moab.setEnData( l_mesh->getTypeEl(),
                      l_tagElTg,
                      l_tsGroups->getElTg() );
  }

  EDGE_V_LOG_INFO << "partitioning the mesh";
  edge_v::mesh::Partition l_part( *l_mesh,
                                   l_tsGroups->getElTg() );
  if( l_config.nPartitions() > 1 ) {
    l_part.kWay( l_config.nPartitions() );
  }
  if( l_config.getWriteElAn() ) {
    EDGE_V_LOG_INFO << "storing elements' partitions";
    std::string l_tagElPa = "edge_v_partitions";
    l_moab.setEnData( l_mesh->getTypeEl(),
                      l_tagElPa,
                      l_part.getElPa() );
  }

  EDGE_V_LOG_INFO << "computing and storing elements' priorities";
  std::string l_tagElPr = "edge_v_element_priorities";
  l_moab.setEnData( l_mesh->getTypeEl(),
                    l_tagElPr,
                    l_part.getElPr() );

  EDGE_V_LOG_INFO << "reordering by rank and time group";
  l_moab.reorder( l_mesh->getTypeEl(),
                  l_tagElPr );

  if( !l_config.getWriteElAn() ) {
    EDGE_V_LOG_INFO << "deleting elements' priorities";
    l_moab.deleteTag( l_tagElPr );
  }

  if( l_config.getWriteElAn() ) {
    EDGE_V_LOG_INFO << "storing element ids";
    edge_v::t_idx * l_elIds = new edge_v::t_idx[ l_mesh->nEls() ];
    for( edge_v::t_idx l_el = 0; l_el < l_mesh->nEls(); l_el++ ) {
      l_elIds[l_el] = l_el;
    }
    std::string l_tagElIds = "edge_v_element_ids";
    l_moab.setEnData( l_mesh->getTypeEl(),
                      l_tagElIds,
                      l_elIds );
    delete[] l_elIds;
  }

  EDGE_V_LOG_INFO << "freeing ts-groups, cfl, mesh and velocity model";
  delete l_tsGroups;
  delete l_cfl;
  delete l_mesh;
  delete l_velMod;

  EDGE_V_LOG_INFO << "re-initializing mesh interface with reordered data";
  l_mesh = new edge_v::mesh::Mesh( l_moab,
                                   l_config.getPeriodic() );

  EDGE_V_LOG_INFO << "re-initializing velocity model";
#ifdef PP_HAS_UCVM
  float * l_velP = new float[l_mesh->nEls()];
  float * l_velS = new float[l_mesh->nEls()];
  float * l_rho  = new float[l_mesh->nEls()];
  if( l_config.getVelModUcvmProjSrc() != "" ) {
    l_velMod = new edge_v::models::seismic::Ucvm( l_ucvmReader,
                                                  l_config.getVelModUcvmTrafoSrc(),
                                                  l_config.getVelModUcvmProjSrc(),
                                                  l_config.getVelModUcvmProjDes(),
                                                  l_config.getVelModUcvmModelType() );
    l_velMod->init( l_mesh->nVes(),
                    l_mesh->getVeCrds() );

    // compute elemnent averages of seismic velocities
    l_velMod->getElAve( l_mesh->nElVes(),
                        l_mesh->nEls(),
                        l_mesh->getElVe(),
                        ( (edge_v::models::seismic::Ucvm*) l_velMod)->getVelP(),
                        l_velP );
    l_velMod->getElAve( l_mesh->nElVes(),
                        l_mesh->nEls(),
                        l_mesh->getElVe(),
                        ( (edge_v::models::seismic::Ucvm*) l_velMod)->getVelS(),
                        l_velS );
    l_velMod->getElAve( l_mesh->nElVes(),
                        l_mesh->nEls(),
                        l_mesh->getElVe(),
                        ( (edge_v::models::seismic::Ucvm*) l_velMod)->getRho(),
                        l_rho );
  }
#else
  if( false ){}
#endif
  else if( l_config.getVelModSeismicExpr() != ""  ) {
    l_velMod = new edge_v::models::seismic::Expression( l_config.getVelModSeismicExpr() );
  }
  else {
    l_velMod = new edge_v::models::Constant( 1 );
  }

#ifdef PP_HAS_UCVM
  if( l_config.getWriteElAn() ) {
    EDGE_V_LOG_INFO << "storing vp, vs and rho";
    l_moab.setEnData( l_mesh->getTypeEl(),
                      l_tagVp,
                      l_velP );
    l_moab.setEnData( l_mesh->getTypeEl(),
                      l_tagVs,
                      l_velS );
    l_moab.setEnData( l_mesh->getTypeEl(),
                      l_tagRho,
                      l_rho );
  }
#endif

  EDGE_V_LOG_INFO << "re-initializing CFL interface with reordered data";
  l_cfl = new edge_v::time::Cfl(  l_mesh->getTypeEl(),
                                  l_mesh->nVes(),
                                  l_mesh->nEls(),
                                  l_mesh->getElVe(),
                                  l_mesh->getVeCrds(),
                                  l_mesh->getInDiasEl(),
                                  *l_velMod );

  EDGE_V_LOG_INFO << "re-initializing time step groups with reordered data";
  l_tsGroups = new edge_v::time::Groups( l_mesh->getTypeEl(),
                                          l_mesh->nEls(),
                                          l_mesh->getElFaEl(),
                                          l_nRates,
                                          l_rates,
                                          l_funDt,
                                          l_cfl->getTimeSteps() );

  std::string l_tagNtgElsIn = "edge_v_n_time_group_elements_inner";
  std::string l_tagNtgElsSe = "edge_v_n_time_group_elements_send";
  std::string l_tagRelTs = "edge_v_relative_time_steps";

  if( l_config.nPartitions() == 1 ) {
    std::string l_pathMesh = l_config.getMeshOutBase() + l_config.getMeshOutExt();
    EDGE_V_LOG_INFO << "writing mesh: " << l_pathMesh;
    l_moab.writeMesh( l_pathMesh );

    EDGE_V_LOG_INFO << "adding meta data: " << l_pathMesh;
    edge_v::io::Hdf5 l_hdf( l_pathMesh,
                            false );

    l_hdf.set( l_tagNtgElsIn,
               l_tsGroups->nGroups(),
               l_tsGroups->nGroupEls() );

    // dummy number of send elements
    edge_v::t_idx * l_nTgElsSe = new edge_v::t_idx[ l_tsGroups->nGroups() ];
    for( unsigned short l_tg = 0; l_tg < l_tsGroups->nGroups(); l_tg++ ) l_nTgElsSe[l_tg] = 0;
    l_hdf.set( l_tagNtgElsSe,
               l_tsGroups->nGroups(),
               l_nTgElsSe );
    delete[] l_nTgElsSe;

    l_hdf.set( l_tagRelTs,
               l_tsGroups->nGroups()+1,
               l_tsGroups->getTsIntervals() );
  }

  if( l_config.nPartitions() > 1 ) {
    EDGE_V_LOG_INFO << "deriving communication structure";
    edge_v::mesh::Communication l_comm( l_tsGroups->nGroups(),
                                        l_mesh->getTypeEl(),
                                        l_mesh->nEls(),
                                        l_mesh->getElFaEl(),
                                        l_part.nPas(),
                                        l_part.nPaEls(),
                                        l_tsGroups->getElTg() );

    EDGE_V_LOG_INFO << "writing mesh by partition";
    edge_v::t_idx l_first = 0;
    for( edge_v::t_idx l_pa = 0; l_pa < l_config.nPartitions(); l_pa++ ) {
      // get number of elements in the partition
      edge_v::t_idx l_nPaEls = l_part.nPaEls()[l_pa];

      std::string l_pathMesh = l_config.getMeshOutBase() + "_" + std::to_string(l_pa) + l_config.getMeshOutExt();
      EDGE_V_LOG_INFO << "  writing mesh: " << l_pathMesh;
      l_moab.writeMesh( l_first,
                        l_nPaEls,
                        l_pathMesh );

      EDGE_V_LOG_INFO << "  adding meta data: " << l_pathMesh;
      edge_v::io::Hdf5 l_hdf( l_pathMesh,
                              false );

#ifdef PP_HAS_UCVM
      if( l_config.getVelModUcvmProjSrc() != "" ) {
        l_hdf.set( l_tagVp,
                   l_nPaEls,
                   l_velP + l_first );
        l_hdf.set( l_tagVs,
                   l_nPaEls,
                   l_velS + l_first );
        l_hdf.set( l_tagRho,
                   l_nPaEls,
                   l_rho + l_first );
      }
#endif

      // store relative time steps of the groups
      l_hdf.set( l_tagRelTs,
                 l_tsGroups->nGroups()+1,
                 l_tsGroups->getTsIntervals() );

      // write number of elements per time group partition-local
      l_hdf.set( l_tagNtgElsIn,
                 l_tsGroups->nGroups(),
                 l_comm.nGroupElsIn( l_pa ) );
      l_hdf.set( l_tagNtgElsSe,
                 l_tsGroups->nGroups(),
                 l_comm.nGroupElsSe( l_pa ) );

      // annotate with comm data
      std::string l_tagCoSt = "edge_v_communication_structure";
      edge_v::t_idx l_coSize = l_comm.getStruct( l_pa )[0]*4 + 1;
      l_hdf.set( l_tagCoSt,
                 l_coSize,
                 l_comm.getStruct( l_pa ) );

      std::string l_tagSeEl = "edge_v_send_el";
      l_hdf.set( l_tagSeEl,
                 l_comm.nSeRe( l_pa ),
                 l_comm.getSendEl( l_pa ) );

      std::string l_tagSeFa = "edge_v_send_fa";
      l_hdf.set( l_tagSeFa,
                 l_comm.nSeRe( l_pa ),
                 l_comm.getSendFa( l_pa ) );

      std::string l_tagReEl = "edge_v_recv_el";
      l_hdf.set( l_tagReEl,
                 l_comm.nSeRe( l_pa ),
                 l_comm.getRecvEl( l_pa ) );

      std::string l_tagReFa = "edge_v_recv_fa";
      l_hdf.set( l_tagReFa,
                 l_comm.nSeRe( l_pa ),
                 l_comm.getRecvFa( l_pa ) );

      // write face- and vertex-ids for ghost elements, adjacent to the send-elements
      unsigned short * l_idsAd = new unsigned short[ l_comm.nSeRe( l_pa ) ];
      l_mesh->getFaIdsAd( l_comm.nSeRe( l_pa ),
                          l_first,
                          l_comm.getSendEl( l_pa ),
                          l_comm.getSendFa( l_pa ),
                          l_idsAd );
      std::string l_tagFaIdsAd = "edge_v_send_face_ids";
      l_hdf.set( l_tagFaIdsAd,
                 l_comm.nSeRe( l_pa ),
                 l_idsAd );

      l_mesh->getVeIdsAd( l_comm.nSeRe( l_pa ),
                          l_first,
                          l_comm.getSendEl( l_pa ),
                          l_comm.getSendFa( l_pa ),
                          l_idsAd );
      std::string l_tagVeIdsAd = "edge_v_send_vertex_ids";
      l_hdf.set( l_tagVeIdsAd,
                 l_comm.nSeRe( l_pa ),
                 l_idsAd );
      delete[] l_idsAd;

      l_first += l_nPaEls;
    }
  }

#ifdef PP_HAS_UCVM
  delete[] l_velP;
  delete[] l_velS;
  delete[] l_rho;
#endif
  delete l_tsGroups;
  delete l_cfl;
  delete l_mesh;
  delete l_velMod;
  delete[] l_rates;

  EDGE_V_LOG_INFO << "thank you for using EDGE-V; exiting..";

  return EXIT_SUCCESS;
}