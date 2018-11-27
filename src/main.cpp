/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2015-2018, Regents of the University of California
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
 * This is the main file of EDGE.
 **/

#include "parallel/Mpi.h"
#include "parallel/Shared.h"

#include "io/logging.h"
#ifdef PP_USE_EASYLOGGING
INITIALIZE_EASYLOGGINGPP
#endif

#include <limits>
#include <string>
#include "io/OptionParser.h"
#include "io/Config.h"
#include "dg/Basis.h"
#include "dg/QuadratureEval.hpp"
#include "sc/Memory.hpp"
#include "sc/Init.hpp"
#include "sc/ibnd/SuperCell.hpp"
#include "io/Receivers.h"
#include "io/ReceiversSf.hpp"
#include "io/WaveField.h"
#include "data/Dynamic.h"
#include "data/DataLayout.hpp"
#include "data/SparseEntities.hpp"
#include "data/Internal.hpp"
#include "time/Manager.h"
#include "mesh/common.hpp"
#include "mesh/SparseTypes.hpp"
#include "mesh/setup_dep.inc"
#include "monitor/Timer.hpp"
#include "monitor/instrument.hpp"

// include dependencies of the setups
#if defined PP_T_EQUATIONS_ADVECTION
#include "impl/advection/setup_dep.inc"
#include "impl/advection/fin_dep.inc"
#elif defined PP_T_EQUATIONS_ELASTIC
#include "impl/elastic/setup_dep.inc"
#include "impl/elastic/fin_dep.inc"
#elif defined PP_T_EQUATIONS_SWE
#include "impl/swe/setup_dep.inc"
#else
#error implementation not supported
#endif

int main( int i_argc, char *i_argv[] ) {
  PP_INSTR_FUN("main")

  // disable logging file-IO
  edge::io::logging::config();

  // create a timer
  edge::monitor::Timer l_timer;
  l_timer.start();
  PP_INSTR_REG_DEF(init)
  PP_INSTR_REG_BEG(init,"init")

  // start shared memory parallelization
  edge::parallel::Shared l_shared;
  l_shared.init();

  // start MPI
  edge::parallel::Mpi l_mpi;
  l_mpi.start( i_argc, i_argv );

  // reconfigure the logging interface with rank and thread id
  edge::io::logging::config();

  EDGE_LOG_INFO << "##########################################################################";
  EDGE_LOG_INFO << "##############   ##############            ###############  ##############";
  EDGE_LOG_INFO << "##############   ###############         ################   ##############";
  EDGE_LOG_INFO << "#####            #####       #####      ######                       #####";
  EDGE_LOG_INFO << "#####            #####        #####    #####                         #####";
  EDGE_LOG_INFO << "#############    #####         #####  #####                  #############";
  EDGE_LOG_INFO << "#############    #####         #####  #####      #########   #############";
  EDGE_LOG_INFO << "#####            #####         #####  #####      #########           #####";
  EDGE_LOG_INFO << "#####            #####        #####    #####        ######           #####";
  EDGE_LOG_INFO << "#####            #####       #####      #####       #####            #####";
  EDGE_LOG_INFO << "###############  ###############         ###############   ###############";
  EDGE_LOG_INFO << "###############  ##############           #############    ###############";
  EDGE_LOG_INFO << "##########################################################################";
  EDGE_LOG_INFO << "";
  EDGE_LOG_INFO << "please come in, have a seat, and.. let's go!!";

#ifdef PP_USE_MPI
  EDGE_LOG_INFO << "our mpi settings:";
  EDGE_LOG_INFO << "  standard-version " << edge::parallel::Mpi::m_verStd[0] << "."
                                         << edge::parallel::Mpi::m_verStd[1];
  EDGE_LOG_INFO << "  #ranks: "          << edge::parallel::g_nRanks;
#endif

#ifdef PP_USE_OMP
  l_shared.print();
#endif

  // print memory information
  edge::data::common::printHugePages();
  edge::data::common::printNumaSizes();
  edge::data::common::printMemStats();

  // parse command line options
  EDGE_LOG_INFO << "parsing command line options";
  edge::io::OptionParser l_options( i_argc, i_argv );

  EDGE_LOG_INFO << "parsing xml config";
  edge::io::Config l_config( l_options.getXmlPath() );

  // parse mesh
  EDGE_LOG_INFO << "parsing mesh";
#include "mesh/setup.inc"

  // get the data layout
  EDGE_LOG_INFO << "taking care of data layout now";
#include "data/setup.inc"

  // initialize all elements/faces
  edge::data::Internal l_internal;
  l_internal.initScratch();
  l_internal.initDense(  l_enLayouts[0].nEnts,
                         l_enLayouts[1].nEnts,
                         l_enLayouts[2].nEnts );

  // setup constant data structures for DG
  EDGE_LOG_INFO << "setting up basis and DG-structure";
  edge::dg::Basis l_basis( T_SDISC.ELEMENT, ORDER );
  l_basis.print();

#include "dg/setup_ader.inc"
#include "sc/setup.inc"

  // initialize internal chars and connectivity information
  EDGE_LOG_INFO << "initializing internal chars and connectivity info";
  l_mesh.getVeChars( l_internal.m_vertexChars  );
  l_mesh.getElChars( l_internal.m_elementChars );
  l_mesh.getFaChars( l_internal.m_faceChars    );
  // int_el l_nElVeEl = l_mesh.getNelVeEl();
  // l_internal.m_connect.elVeEl[0] = (int_el*) l_dynMem.allocate( l_nElVeEl * sizeof(int_el) );
  l_mesh.getConnect( l_internal.m_vertexChars,
                     l_internal.m_faceChars,
                     l_internal.m_connect      );

  // enhance entity chars if set in the config
  if( l_config.m_spTypesDoms[0].size() > 0 ) EDGE_LOG_FATAL << "not implemented";
  if( l_config.m_spTypesDoms[1].size() > 0 ) edge::mesh::SparseTypes<
                                               T_SDISC.FACE
                                             >::set(  l_enLayouts[1].nEnts,
                                                      l_internal.m_connect.faVe,
                                                     &l_config.m_spTypesVals[1][0],
                                                      l_config.m_spTypesDoms[1],
                                                      l_internal.m_vertexChars,
                                                      l_internal.m_faceChars );
  if( l_config.m_spTypesDoms[2].size() > 0 ) edge::mesh::SparseTypes<
                                               T_SDISC.ELEMENT
                                             >::set(  l_enLayouts[2].nEnts,
                                                      l_internal.m_connect.elVe,
                                                     &l_config.m_spTypesVals[2][0],
                                                      l_config.m_spTypesDoms[2],
                                                      l_internal.m_vertexChars,
                                                      l_internal.m_elementChars );

  EDGE_VLOG(2) << "  printing neigh relations (loc_fa-nei_fa-nei_ve):";
  if (EDGE_VLOG_IS_ON(2)) {
    edge::mesh::common< T_SDISC.ELEMENT> ::printNeighRel( l_enLayouts[2],
                                                          l_internal.m_connect.fIdElFaEl[0],
                                                          l_internal.m_connect.vIdElFaEl[0] );
  }

//TODO: add to internal
std::vector< int_gid > l_gIdsEl;
l_mesh.getGIdsEl( l_gIdsEl );

  // setup receivers
#include "io/inc/setup_recv.inc"

  // time step statistics
  double l_dT[3];

  EDGE_LOG_INFO << "performing equation-specific setup";
  PP_INSTR_REG_DEF(equSpe)
  PP_INSTR_REG_BEG(equSpe,"eq_spec_setup")
#if defined PP_T_EQUATIONS_ADVECTION
#include "impl/advection/setup.inc"
#elif defined PP_T_EQUATIONS_ELASTIC
#include "impl/elastic/setup.inc"
#elif defined PP_T_EQUATIONS_SWE
#include "impl/swe/setup.inc"
#endif
  PP_INSTR_REG_END(equSpe)

  // determine global time step stats
  double l_dtG[3];
#ifdef PP_USE_MPI
  MPI_Allreduce( l_dT,   l_dtG,   1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
  MPI_Allreduce( l_dT+1, l_dtG+1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  MPI_Allreduce( l_dT+2, l_dtG+2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
#else
  l_dtG[0] = l_dT[0];
  l_dtG[1] = l_dT[1];
  l_dtG[2] = l_dT[2];
#endif
  l_dtG[1] /= edge::parallel::g_nRanks;

  // construct single GTS cluster
  edge::time::TimeGroupStatic l_cluster( std::numeric_limits< int_ts >::max(),
                                         1,
                                         l_internal );

  EDGE_LOG_INFO << "time step stats coming thru (min,ave,max): "
                << l_dtG[0] << ", " << l_dtG[1] << ", " << l_dtG[2];

  // add cluster to time manager
  edge::time::Manager l_time( l_dtG[0], l_shared, l_mpi, l_receivers, l_recvsSf );
  l_time.add( &l_cluster );

  // set up simulation times and synchronization intervals
  double l_simTime = 0;
  double l_endTime = l_config.m_endTime;
  double l_syncInt = l_config.m_waveFieldInt;
         l_syncInt = std::min( l_syncInt, l_config.m_iBndInt );

  if( std::abs(l_syncInt) < TOL.TIME ) l_syncInt = l_endTime;

  // create a wave field writer
  edge::io::WaveField l_writer( l_config.m_waveFieldType,
                                l_config.m_waveFieldFile,
                                l_enLayouts[2],
                                l_mesh.getInMap(),
                                l_internal.m_vertexChars,
                                l_internal.m_elementChars,
                                l_internal.m_connect.elVe,
                                l_internal.m_elementModePrivate1,
                                l_config.m_waveFieldSpType );

  // write setup
  EDGE_LOG_INFO << "reached synchronization point #0: " << l_simTime;
  if( l_config.m_waveFieldInt > TOL.TIME ) {
    EDGE_LOG_INFO << "  writing wave field #0";
    l_writer.write( 0,
                    l_internal.m_globalShared2[0].limSync );
  }

#if defined(PP_T_EQUATIONS_ELASTIC)
  if( l_config.m_iBndInt > TOL.TIME ) {
    EDGE_LOG_INFO << "  writing internal boundary #0";
    l_rupWriter.write( 0,
                       l_enLayouts[l_rupLayoutFa].nEnts,
                       5*N_CRUNS + 4 * (N_DIM-1) * N_CRUNS,
                       5*N_CRUNS + 4 * (N_DIM-1) * N_CRUNS,
                       l_internal.m_globalShared5[0].sfQtNaPtr,
                       (real_base*) l_internal.m_globalShared5[0].sf );
  }
#endif

  // print mem stats
  edge::data::common::printMemStats();

  // print timing info for init
  l_timer.end();
  PP_INSTR_REG_END(init)
  EDGE_LOG_INFO << "initialization phase took us " << l_timer.elapsed() << " seconds";

  PP_INSTR_REG_DEF(comp)
#ifdef PP_USE_MPI
  int l_err = MPI_Barrier( MPI_COMM_WORLD );
  EDGE_CHECK_EQ( l_err, MPI_SUCCESS );
#endif
  PP_INSTR_REG_BEG(comp,"comp")
  l_timer.start();

  // iterate over sync points
  unsigned int l_step    = 0;
  unsigned int l_stepWf  = 0;
  unsigned int l_stepBnd = 0;
  while( l_endTime - l_simTime > TOL.TIME ) {
    // derive time to advance in this step
    double l_stepTime = std::max( 0.0, l_endTime - l_simTime );
           l_stepTime = std::min( l_stepTime, l_syncInt );

    l_time.simulate( l_stepTime );

    // update simulation time
    l_simTime += l_stepTime;

    EDGE_LOG_INFO << "reached synchronization point #" << l_step << ": " << l_simTime;

    // write this sync step
    if( l_simTime + TOL.TIME > (l_stepWf+1)*l_config.m_waveFieldInt ) {
      EDGE_LOG_INFO << "  writing wave field #" << l_stepWf+1;
      l_writer.write( l_stepTime,
                      l_internal.m_globalShared2[0].limSync );
      l_stepWf++;
    }

#if defined(PP_T_EQUATIONS_ELASTIC)
    if( l_simTime + TOL.TIME > (l_stepBnd+1)*l_config.m_iBndInt ) {
      EDGE_LOG_INFO << "  writing internal boundary #" << l_stepBnd+1;
      l_rupWriter.write( 0,
                         l_enLayouts[l_rupLayoutFa].nEnts,
                         2*N_CRUNS + 3 * (N_DIM-1) * N_CRUNS,
                         5*N_CRUNS + 4 * (N_DIM-1) * N_CRUNS,
                         l_internal.m_globalShared5[0].sfQtNaPtr + (N_DIM+2)*N_CRUNS,
                         (real_base*) l_internal.m_globalShared5[0].sf[0][0].muf );
      l_stepBnd++;
    }
#endif

    // increase step and derive next synchronization point
    l_step++;
    if( (l_stepWf+1)*l_config.m_waveFieldInt < (l_stepBnd+1)*l_config.m_iBndInt )
      l_syncInt = (l_stepWf +1)*l_config.m_waveFieldInt - l_simTime;
    else
      l_syncInt = (l_stepBnd+1)*l_config.m_iBndInt      - l_simTime;

    if( std::abs(l_syncInt) < TOL.TIME ) l_syncInt = l_endTime;
  }

  // print time info for compute
  l_timer.end();
  PP_INSTR_REG_END(comp)
  EDGE_LOG_INFO << "that's the duration of the computations ("
                << l_cluster.getUpdatesPer() << " time steps): "
                << l_timer.elapsed() << " seconds";
  PP_INSTR_REG_DEF(fin)
  PP_INSTR_REG_BEG(fin,"fin")
  l_timer.start();

#if defined PP_T_EQUATIONS_ADVECTION
#include "impl/advection/fin.inc"
#elif defined PP_T_EQUATIONS_ELASTIC
#include "impl/elastic/fin.inc"
#endif

  // shutdown internal structure
  l_internal.finalize();

  EDGE_LOG_INFO << "that was fun: EDGE over and out!";

  // stop MPI
  l_mpi.fin();

  // print duration of finalizaiton
  l_timer.end();
  PP_INSTR_REG_END(fin)
  EDGE_LOG_INFO << "finalizing time: " << l_timer.elapsed();
}
