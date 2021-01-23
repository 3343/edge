/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2020-2021, Friedrich Schiller University Jena
 * Copyright (c) 2019-2020, Alexander Breuer
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
#if defined(PP_USE_GASPI)
#include "parallel/Gaspi.h"
#elif defined(PP_USE_MPI)
#include "parallel/MpiRemix.h"
#else
#include "parallel/DistributedDummy.hpp"
#endif

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
#include "setups/Cpu.h"
#include "io/Receivers.h"
#include "io/WaveField.h"
#include "data/Dynamic.h"
#include "data/DataLayout.hpp"
#include "data/SparseEntities.hpp"
#include "data/Internal.hpp"
#include "time/Manager.h"
#include "mesh/common.hpp"
#include "mesh/SparseTypes.hpp"
#include "mesh/EdgeV.h"

#include "monitor/Timer.hpp"
#include "monitor/instrument.hpp"

// include dependencies of the setups
#if defined PP_T_EQUATIONS_ADVECTION
#include "impl/advection/setup_dep.inc"
#include "impl/advection/fin_dep.inc"
#elif defined PP_T_EQUATIONS_SEISMIC
#include "impl/seismic/setup_dep.inc"
#include "impl/seismic/fin_dep.inc"
#elif defined PP_T_EQUATIONS_SWE
#include "impl/swe/setup_dep.inc"
#else
#error implementation not supported
#endif

int main( int i_argc, char *i_argv[] ) {
  PP_INSTR_FUN("main")

  // set CPU flags
  edge::setups::Cpu::setFlushToZero( true );
  edge::setups::Cpu::setDenormalsAreZero( true );

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

  // init distributed memory parallelization
#if defined(PP_USE_GASPI)
  edge::parallel::Gaspi l_distributed( i_argc,
                                       i_argv );
#elif defined(PP_USE_MPI)
  edge::parallel::MpiRemix l_distributed( i_argc,
                                          i_argv );
#else
  edge::parallel::DistributedDummy l_distributed( i_argc,
                                                  i_argv );
#endif

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
  EDGE_LOG_INFO;
  EDGE_LOG_INFO << "               EDGE is available from: https://dial3343.org";
  EDGE_LOG_INFO;

#ifdef PP_USE_MPI
  EDGE_LOG_INFO << "our mpi settings:";
  EDGE_LOG_INFO << "  standard-version " << l_distributed.getVerStr();
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

  // parse mesh and mesh supplement
  EDGE_LOG_INFO << "parsing mesh and supplement";
  std::string l_meshPath = l_config.m_meshInBase;
  l_meshPath += "_" + std::to_string(edge::parallel::g_rank+1);
  l_meshPath += l_config.m_meshInExt;

  std::string l_supplementPath = l_config.m_meshInBase;
  l_supplementPath += "_" + std::to_string(edge::parallel::g_rank+1) + ".h5";

  edge::mesh::EdgeV l_edgeV( l_meshPath,
                             l_supplementPath,
                             l_config.m_periodic );

  // dynamic memory allocates
  edge::data::Dynamic l_dynMem;

  // initialize all elements/faces
  edge::data::Internal l_internal;
  l_internal.initScratch();
  l_internal.initDense(  l_edgeV.nVes(),
                         l_edgeV.nFas(),
                         l_edgeV.nEls() );

  // setup constant data structures for DG
  EDGE_LOG_INFO << "setting up basis and DG-structure";
  edge::dg::Basis l_basis( T_SDISC.ELEMENT, ORDER );
  l_basis.print();

#include "dg/setup_ader.inc"

// initialize internal chars and connectivity information
EDGE_LOG_INFO << "initializing internal chars and connectivity info";
l_internal.m_connect.faVe   = (std::size_t (*)[C_ENT[T_SDISC.ELEMENT].N_FACE_VERTICES]) l_edgeV.getFaVe();
l_internal.m_connect.faEl   = (std::size_t (*)[2])                                      l_edgeV.getFaEl();
l_internal.m_connect.elVe   = (std::size_t (*)[C_ENT[T_SDISC.ELEMENT].N_VERTICES])      l_edgeV.getElVe();
l_internal.m_connect.elFa   = (std::size_t (*)[C_ENT[T_SDISC.ELEMENT].N_FACES])         l_edgeV.getElFa();
l_internal.m_connect.elFaEl = (std::size_t (*)[C_ENT[T_SDISC.ELEMENT].N_FACES])         l_edgeV.getElFaEl();

double const (* l_veCrds)[3] = l_edgeV.getVeCrds();
double const * l_volEl = l_edgeV.getVolumesEl();
double const * l_volFa = l_edgeV.getAreasFa();
double const * l_inDiaEl = l_edgeV.getInDiasEl();
double const (* l_normalsFa)[3] = l_edgeV.getNormalsFa();
double const (* l_tangentsFa)[2][3] = l_edgeV.getTangentsFa();

// fill internal data structures
for( std::size_t l_ve = 0; l_ve < l_edgeV.nVes(); l_ve++ ) {
  for( unsigned short l_di = 0; l_di < 3; l_di++ )
    l_internal.m_vertexChars[l_ve].coords[l_di] = l_veCrds[l_ve][l_di];

  l_internal.m_vertexChars[l_ve].spType = 0;
}

for( std::size_t l_fa = 0; l_fa < l_edgeV.nFas(); l_fa++ ) {
  l_internal.m_faceChars[l_fa].spType = 0;
  l_internal.m_faceChars[l_fa].area = l_volFa[l_fa];

  for( unsigned short l_di = 0; l_di < 3; l_di++ ) {
    l_internal.m_faceChars[l_fa].outNormal[l_di] = l_normalsFa[l_fa][l_di];
    l_internal.m_faceChars[l_fa].tangent0[l_di] = l_tangentsFa[l_fa][0][l_di];
    l_internal.m_faceChars[l_fa].tangent1[l_di] = l_tangentsFa[l_fa][1][l_di];
  }
}

for( std::size_t l_el = 0; l_el < l_edgeV.nEls(); l_el++ ) {
  l_internal.m_elementChars[l_el].spType = 0;
  l_internal.m_elementChars[l_el].volume = l_volEl[l_el];
  l_internal.m_elementChars[l_el].inDia = l_inDiaEl[l_el];
}

l_edgeV.setSpTypes( l_internal.m_vertexChars,
                    l_internal.m_faceChars,
                    l_internal.m_elementChars );

  l_edgeV.setSeVeFaIdsAd( l_internal.m_connect.vIdElFaEl[0],
                          l_internal.m_connect.fIdElFaEl[0] );

  l_edgeV.setLtsTypes( l_internal.m_elementChars );

  // enhance entity chars if set in the config
  if( l_config.m_spTypesDoms[0].size() > 0 ) EDGE_LOG_FATAL << "not implemented";
  if( l_config.m_spTypesDoms[1].size() > 0 ) edge::mesh::SparseTypes<
                                               T_SDISC.FACE
                                             >::set(  l_edgeV.nFas(),
                                                      l_internal.m_connect.faVe,
                                                     &l_config.m_spTypesVals[1][0],
                                                      l_config.m_spTypesDoms[1],
                                                      l_internal.m_vertexChars,
                                                      l_internal.m_faceChars );
  if( l_config.m_spTypesDoms[2].size() > 0 ) edge::mesh::SparseTypes<
                                               T_SDISC.ELEMENT
                                             >::set(  l_edgeV.nEls(),
                                                      l_internal.m_connect.elVe,
                                                     &l_config.m_spTypesVals[2][0],
                                                      l_config.m_spTypesDoms[2],
                                                      l_internal.m_vertexChars,
                                                      l_internal.m_elementChars );

  EDGE_VLOG(2) << "  printing neigh relations (loc_fa-nei_fa-nei_ve):";
  if (EDGE_VLOG_IS_ON(2)) {
    edge::mesh::common< T_SDISC.ELEMENT> ::printNeighRel( l_edgeV.nEls(),
                                                          l_internal.m_connect.fIdElFaEl[0],
                                                          l_internal.m_connect.vIdElFaEl[0] );
  }

  // setup receivers
#include "io/inc/setup_recv.inc"

  // time step statistics
  double l_dT[3];

  EDGE_LOG_INFO << "performing equation-specific setup";
  PP_INSTR_REG_DEF(equSpe)
  PP_INSTR_REG_BEG(equSpe,"eq_spec_setup")
#if defined PP_T_EQUATIONS_ADVECTION
#include "impl/advection/setup.inc"
#elif defined PP_T_EQUATIONS_SEISMIC
#include "impl/seismic/setup.inc"
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

  // assemble LTS groups
  std::vector< edge::time::TimeGroupStatic > l_tgs;

  for( unsigned short l_tg = 0; l_tg < l_edgeV.nTgs(); l_tg++ ) {
    l_tgs.push_back( edge::time::TimeGroupStatic( l_edgeV.nTgs(),
                                                  l_tg,
                                                  l_internal,
                                                  l_distributed.nCommBuffers(),
                                                  l_distributed.getSendPtrs(),
                                                  l_distributed.getRecvPtrs() ) );
  }

  EDGE_LOG_INFO << "time step stats coming thru (min,ave,max): "
                << l_dtG[0] << ", " << l_dtG[1] << ", " << l_dtG[2];

  // create time manager
  EDGE_LOG_INFO << "creating time manager, groups relative to the min. dt:";
  for( unsigned short l_tg = 0; l_tg < l_edgeV.nTgs(); l_tg++ ) {
    EDGE_LOG_INFO << "  [" << l_edgeV.getRelDt()[l_tg] << ", " << l_edgeV.getRelDt()[l_tg+1] << "[";
  }
  edge::time::Manager l_time( l_edgeV.getRelDt()[0]*l_dtG[0],
                              l_shared,
                              l_distributed,
                              l_tgs,
                              l_receivers );

  // set up simulation times and synchronization intervals
  double l_simTime = 0;
  double l_endTime = l_config.m_endTime;
  double l_syncInt = l_config.m_waveFieldInt;
         l_syncInt = std::min( l_syncInt, l_config.m_syncMaxInt );

  if( std::abs(l_syncInt) < TOL.TIME ) l_syncInt = l_endTime;

  // create a wave field writer
  edge::io::WaveField l_writer( l_config.m_waveFieldType,
                                l_config.m_waveFieldFile,
                                l_edgeV.nVes(),
                                l_edgeV.nEls(),
                                l_internal.m_vertexChars,
                                l_internal.m_elementChars,
                                l_internal.m_connect.elVe,
                                l_internal.m_elementModePrivate1,
                                l_config.m_waveFieldSpType );

  // write setup
  EDGE_LOG_INFO << "reached synchronization point #0";
  EDGE_LOG_INFO << "  simulation time: " << l_simTime;
  if( l_config.m_waveFieldInt < l_config.m_endTime ) {
    EDGE_LOG_INFO << "  writing wave field #0";
    l_writer.write( 0 );
  }

  // print mem stats
  edge::data::common::printMemStats();

  // print timing info for init
  l_timer.end();
  PP_INSTR_REG_END(init)
  EDGE_LOG_INFO << "initialization phase took us " << l_timer.elapsed() << " seconds";
  l_timer.reset();

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
  while( l_endTime - l_simTime > TOL.TIME ) {
    // derive time to advance in this step
    double l_stepTime = std::max( 0.0, l_endTime - l_simTime );
           l_stepTime = std::min( l_stepTime, l_syncInt );


    PP_INSTR_REG_DEF( sync )
    PP_INSTR_REG_BEG( sync, "sync" )
#pragma warning push
#pragma warning(disable:68)
    PP_INSTR_PAR_UINT64("sync_id",  (uint64_t) l_step )
#pragma warning pop

    EDGE_LOG_INFO << "progressing simulation by " << l_stepTime;
    l_time.simulate( l_stepTime );
    PP_INSTR_REG_END( sync )

    // update simulation time
    l_simTime += l_stepTime;

    EDGE_LOG_INFO << "reached synchronization point #" << l_step+1;
    EDGE_LOG_INFO << "  simulation time: " << l_simTime;
    l_timer.end();
    EDGE_LOG_INFO << "  estimated remaining time: " << (l_endTime / l_simTime - 1.0) * l_timer.elapsed() << " seconds";
    l_timer.start();

    // write this sync step
    if( l_simTime + TOL.TIME > (l_stepWf+1)*l_config.m_waveFieldInt ) {
      EDGE_LOG_INFO << "  writing wave field #" << l_stepWf+1;
      l_writer.write( l_stepTime );
      l_stepWf++;
    }

    // increase step and derive next synchronization point
    l_step++;
    l_syncInt = (l_stepWf +1)*l_config.m_waveFieldInt - l_simTime;

    l_syncInt = std::min( l_syncInt, l_config.m_syncMaxInt );
    l_syncInt = std::min( l_syncInt, l_endTime-l_simTime );

    if( l_syncInt < TOL.TIME ) l_syncInt = l_endTime;
  }

  // print time info for compute
  l_timer.end();
  PP_INSTR_REG_END(comp)
  EDGE_LOG_INFO << "that's the duration of the computations ("
                << l_tgs[0].getUpdatesPer() << " fundamental time steps): "
                << l_timer.elapsed() << " seconds";
  l_timer.reset();
  PP_INSTR_REG_DEF(fin)
  PP_INSTR_REG_BEG(fin,"fin")
  l_timer.start();

#if defined PP_T_EQUATIONS_ADVECTION
#include "impl/advection/fin.inc"
#elif defined PP_T_EQUATIONS_SEISMIC
#include "impl/seismic/fin.inc"
#endif

  // shutdown internal structure
  l_internal.finalize();

  EDGE_LOG_INFO << "that was fun: EDGE over and out!";

  // stop MPI
  l_distributed.fin();

  // print duration of finalizaiton
  l_timer.end();
  PP_INSTR_REG_END(fin)
  EDGE_LOG_INFO << "finalizing time: " << l_timer.elapsed();
  l_timer.reset();
}