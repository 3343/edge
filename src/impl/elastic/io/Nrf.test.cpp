/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017, Regents of the University of California
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
 * Unit tests for the NetCDF Rupture Format reader.
 **/

#include <catch.hpp>
#define private public
#include "Nrf.h"
#undef private

#include <fstream>

TEST_CASE( "2D Nrf: Init, updateBuf.", "[nrf2d][init]" ) {
  std::string l_ncFile = "cont/unit_tests/elastic/sources/nrf_2d.nc";

  // return silently if file does not exist
  if( !std::ifstream(l_ncFile) ) return;


  // NRF reader with 2-source buffer
  edge::elastic::io::Nrf< 2 > l_nrf1( 2 );

  // init
  l_nrf1.init( l_ncFile );

  // file info
  REQUIRE( l_nrf1.m_bSize == 2 );
  REQUIRE( l_nrf1.m_nrf[0]->path   == l_ncFile );
  REQUIRE( l_nrf1.m_nrf[0]->nSrcsG == 3 );
  REQUIRE( l_nrf1.m_nrf[0]->nSplsG[0] == 0  );
  REQUIRE( l_nrf1.m_nrf[0]->nSplsG[1] == 31 );
  REQUIRE( l_nrf1.m_nrf[0]->gIdFb == 0 );

  // src info
  REQUIRE( l_nrf1.m_nrf[0]->crds[0][0] == Approx( 913609.762166668 ) );
  REQUIRE( l_nrf1.m_nrf[0]->crds[0][1] == Approx( 3795426.22805577 ) );

  REQUIRE( l_nrf1.m_nrf[0]->crds[1][0] == Approx( 913609.593554401 ) );
  REQUIRE( l_nrf1.m_nrf[0]->crds[1][1] == Approx( 3795526.32237143 ) );

  // subfault: 1st source
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].tInit == Approx( 10.5201 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dt    == Approx( 0.1 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].mu    == Approx( 3168000000 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].area  == Approx( 12222 ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[0][0] == Approx(  0.911320106822418  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[0][1] == Approx( -0.063725709730111  ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[1][0] == Approx( -0.0838743060458509 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[1][1] == Approx( -0.995966174736174  ) );

  // subfault: 2nd source
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].tInit == Approx( 10.5255 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dt    == Approx( 0.2 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].mu    == Approx( 3168111111 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].area  == Approx( 11111 ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dir[0][0] == Approx(  0.912988950432977  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dir[0][1] == Approx( -0.0318822766865827 ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dir[1][0] == Approx( -0.0911321737575558 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dir[1][1] == Approx( -0.987689275364108  ) );

  // offsets
  REQUIRE( l_nrf1.m_nrf[0]->srOff[0][0] ==  0 );
  REQUIRE( l_nrf1.m_nrf[0]->srOff[0][1] ==  0 );

  REQUIRE( l_nrf1.m_nrf[0]->srOff[1][0] ==  0 );
  REQUIRE( l_nrf1.m_nrf[0]->srOff[1][1] == 11 );

  REQUIRE( l_nrf1.m_nrf[0]->srOff[2][0] ==  0 ); // ghost
  REQUIRE( l_nrf1.m_nrf[0]->srOff[2][1] == 21 ); // ghost

  // samples
  REQUIRE( l_nrf1.m_nrf[0]->nSrMax[0] == 0 );
  REQUIRE( l_nrf1.m_nrf[0]->nSrMax[1] == (std::size_t) (21*l_nrf1.m_srOvAll) );

  REQUIRE( l_nrf1.m_nrf[0]->sr[1][ 0] == Approx( 0.0       ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][ 1] == Approx( 0.332642  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][ 2] == Approx( 0.289227  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][ 3] == Approx( 0.0753253 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][ 4] == Approx( 0.0669723 ) );
  // lazy: ignore middle part
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][19] == Approx( 0.0177707  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][20] == Approx( 0.00595686 ) );


  // update the buffer with the last entry
  l_nrf1.updateBufAll( 2 );

  // src info
  REQUIRE( l_nrf1.m_nrf[0]->crds[0][0] == Approx( 913593.932784665 ) );
  REQUIRE( l_nrf1.m_nrf[0]->crds[0][1] == Approx( 3795625.7327648  ) );

  // subfault
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].tInit == Approx( 10.4184 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dt    == Approx( 0.3 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].mu    == Approx( 3168222222 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].area  == Approx( 10000 ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[0][0] == Approx(  0.916427924584317  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[0][1] == Approx(  0.178135543283137  ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[1][0] == Approx(  0.100036968135533 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[1][1] == Approx( -0.969011238609667  ) );

  // offsets
  REQUIRE( l_nrf1.m_nrf[0]->srOff[0][0] ==  0 );
  REQUIRE( l_nrf1.m_nrf[0]->srOff[0][1] ==  0 );

  REQUIRE( l_nrf1.m_nrf[0]->srOff[1][0] ==  0 ); // ghost
  REQUIRE( l_nrf1.m_nrf[0]->srOff[1][1] == 10 ); // ghost

  // samples
  REQUIRE( l_nrf1.m_nrf[0]->nSrMax[0] == 0 );
  REQUIRE( l_nrf1.m_nrf[0]->nSrMax[1] == (std::size_t) (21*l_nrf1.m_srOvAll) );

  REQUIRE( l_nrf1.m_nrf[0]->sr[1][0] == Approx( 0.0        ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][1] == Approx( 0.647666   ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][2] == Approx( 0.383135   ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][3] == Approx( 0.130604   ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][4] == Approx( 0.111769   ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][5] == Approx( 0.0877252  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][6] == Approx( 0.0616282  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][7] == Approx( 0.0369034  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][8] == Approx( 0.0167959  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][9] == Approx( 0.00394474 ) );
}

TEST_CASE( "3D Nrf: Init, updateBuf.", "[nrf3d][init]" ) {
  std::string l_ncFile = "cont/unit_tests/elastic/sources/nrf_3d.nc";

  // return silently if file does not exist
  // we assume that the unit tests executable is simply run from a different location
  if( !std::ifstream(l_ncFile) ) return;

  // NRF reader with 2-source buffer
  edge::elastic::io::Nrf< 3 > l_nrf1( 2 );

  // init
  l_nrf1.init( l_ncFile );

  // file info
  REQUIRE( l_nrf1.m_bSize == 2 );
  REQUIRE( l_nrf1.m_nrf[0]->path   == l_ncFile );
  REQUIRE( l_nrf1.m_nrf[0]->nSrcsG == 3 );
  REQUIRE( l_nrf1.m_nrf[0]->nSplsG[0] == 0  );
  REQUIRE( l_nrf1.m_nrf[0]->nSplsG[1] == 31 );
  REQUIRE( l_nrf1.m_nrf[0]->nSplsG[2] == 0  );
  REQUIRE( l_nrf1.m_nrf[0]->gIdFb == 0 );

  // src info
  REQUIRE( l_nrf1.m_nrf[0]->crds[0][0] == Approx( 913609.762166668 ) );
  REQUIRE( l_nrf1.m_nrf[0]->crds[0][1] == Approx( 3795426.22805577 ) );
  REQUIRE( l_nrf1.m_nrf[0]->crds[0][2] == Approx( -50.0 ) );

  REQUIRE( l_nrf1.m_nrf[0]->crds[1][0] == Approx( 913609.593554401 ) );
  REQUIRE( l_nrf1.m_nrf[0]->crds[1][1] == Approx( 3795526.32237143 ) );
  REQUIRE( l_nrf1.m_nrf[0]->crds[1][2] == Approx( -50.0 ) );

  // subfault: 1st source
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].tInit == Approx( 10.5201 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dt    == Approx( 0.1 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].mu    == Approx( 3168000000 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].area  == Approx( 12222 ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[0][0] == Approx(  0.911320106822418  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[0][1] == Approx( -0.063725709730111  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[0][2] == Approx(  0.4067366430758    ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[1][0] == Approx( -0.0838743060458509 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[1][1] == Approx( -0.995966174736174  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[1][2] == Approx(  0.0318822766865828 ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[2][0] == Approx(  0.403064217819573  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[2][1] == Approx( -0.0631697134771603 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[2][2] == Approx( -0.912988950432977  ) );

  // subfault: 2nd source
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].tInit == Approx( 10.5255 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dt    == Approx( 0.2 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].mu    == Approx( 3168111111 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].area  == Approx( 11111 ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dir[0][0] == Approx(  0.912988950432977  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dir[0][1] == Approx( -0.0318822766865827 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dir[0][2] == Approx(  0.4067366430758    ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dir[1][0] == Approx( -0.0911321737575558 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dir[1][1] == Approx( -0.987689275364108  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dir[1][2] == Approx(   0.127140954208103 ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dir[2][0] == Approx(  0.397675877183308  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dir[2][1] == Approx( -0.153145080769852  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[1].dir[2][2] == Approx( -0.904654896047372  ) );

  // offsets
  REQUIRE( l_nrf1.m_nrf[0]->srOff[0][0] ==  0 );
  REQUIRE( l_nrf1.m_nrf[0]->srOff[0][1] ==  0 );
  REQUIRE( l_nrf1.m_nrf[0]->srOff[0][2] ==  0 );

  REQUIRE( l_nrf1.m_nrf[0]->srOff[1][0] ==  0 );
  REQUIRE( l_nrf1.m_nrf[0]->srOff[1][1] == 11 );
  REQUIRE( l_nrf1.m_nrf[0]->srOff[1][2] ==  0 );

  REQUIRE( l_nrf1.m_nrf[0]->srOff[2][0] ==  0 ); // ghost
  REQUIRE( l_nrf1.m_nrf[0]->srOff[2][1] == 21 ); // ghost
  REQUIRE( l_nrf1.m_nrf[0]->srOff[2][2] ==  0 ); // ghost

  // samples
  REQUIRE( l_nrf1.m_nrf[0]->nSrMax[0] == 0 );
  REQUIRE( l_nrf1.m_nrf[0]->nSrMax[1] == (std::size_t) (21*l_nrf1.m_srOvAll) );

  REQUIRE( l_nrf1.m_nrf[0]->sr[1][ 0] == Approx( 0.0       ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][ 1] == Approx( 0.332642  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][ 2] == Approx( 0.289227  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][ 3] == Approx( 0.0753253 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][ 4] == Approx( 0.0669723 ) );
  // lazy: ignore middle part
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][19] == Approx( 0.0177707  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][20] == Approx( 0.00595686 ) );


  // update the buffer with the last entry
  l_nrf1.updateBufAll( 2 );

  // src info
  REQUIRE( l_nrf1.m_nrf[0]->crds[0][0] == Approx( 913593.932784665 ) );
  REQUIRE( l_nrf1.m_nrf[0]->crds[0][1] == Approx( 3795625.7327648  ) );
  REQUIRE( l_nrf1.m_nrf[0]->crds[0][2] == Approx( -50.0 ) );

  // subfault
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].tInit == Approx( 10.4184 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dt    == Approx( 0.3 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].mu    == Approx( 3168222222 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].area  == Approx( 10000 ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[0][0] == Approx(  0.916427924584317  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[0][1] == Approx(  0.178135543283137  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[0][2] == Approx(  0.3583679495453    ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[1][0] == Approx(  0.100036968135533 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[1][1] == Approx( -0.969011238609667  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[1][2] == Approx(  0.225853546472949 ) );

  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[2][0] == Approx(  0.38749511487028  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[2][1] == Approx( -0.171128453704753 ) );
  REQUIRE( l_nrf1.m_nrf[0]->sub[0].dir[2][2] == Approx( -0.905849097965157  ) );

  // offsets
  REQUIRE( l_nrf1.m_nrf[0]->srOff[0][0] ==  0 );
  REQUIRE( l_nrf1.m_nrf[0]->srOff[0][1] ==  0 );
  REQUIRE( l_nrf1.m_nrf[0]->srOff[0][2] ==  0 );

  REQUIRE( l_nrf1.m_nrf[0]->srOff[1][0] ==  0 ); // ghost
  REQUIRE( l_nrf1.m_nrf[0]->srOff[1][1] == 10 ); // ghost
  REQUIRE( l_nrf1.m_nrf[0]->srOff[1][2] ==  0 ); // ghost

  // samples
  REQUIRE( l_nrf1.m_nrf[0]->nSrMax[0] == 0 );
  REQUIRE( l_nrf1.m_nrf[0]->nSrMax[1] == (std::size_t) (21*l_nrf1.m_srOvAll) );

  REQUIRE( l_nrf1.m_nrf[0]->sr[1][0] == Approx( 0.0        ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][1] == Approx( 0.647666   ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][2] == Approx( 0.383135   ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][3] == Approx( 0.130604   ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][4] == Approx( 0.111769   ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][5] == Approx( 0.0877252  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][6] == Approx( 0.0616282  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][7] == Approx( 0.0369034  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][8] == Approx( 0.0167959  ) );
  REQUIRE( l_nrf1.m_nrf[0]->sr[1][9] == Approx( 0.00394474 ) );
}

TEST_CASE( "2D Nrf: srcCrds", "[nrf2d][srcCrds]" ) {
  std::string l_ncFile = "cont/unit_tests/elastic/sources/nrf_2d.nc";

  // return silently if file does not exist
  if( !std::ifstream(l_ncFile) ) return;

  // NRF reader with 2-source buffer
  edge::elastic::io::Nrf< 2 > l_nrf1( 2 );

  // init
  l_nrf1.init( l_ncFile );
  l_nrf1.init( l_ncFile );

  // check basic quantities
  REQUIRE( l_nrf1.nIn() == 2 );
  REQUIRE( l_nrf1.nSrcsG( 1 ) == 3 );

  // source coordinates
  double l_srcCrds[3][2];

  // get source coordinates without arguments
  l_nrf1.getSrcCrds( 1, l_srcCrds );

  REQUIRE( l_srcCrds[0][0] == Approx( 913609.762166668 ) );
  REQUIRE( l_srcCrds[0][1] == Approx( 3795426.22805577 ) );

  REQUIRE( l_srcCrds[1][0] == Approx( 913609.593554401 ) );
  REQUIRE( l_srcCrds[1][1] == Approx( 3795526.32237143 ) );

  REQUIRE( l_srcCrds[2][0] == Approx( 913593.932784665 ) );
  REQUIRE( l_srcCrds[2][1] == Approx( 3795625.7327648  ) );

  // get source coordinates of #2 only
  l_nrf1.getSrcCrds( 1, l_srcCrds, 1, 1 );
  REQUIRE( l_srcCrds[0][0] == Approx( 913609.593554401 ) );
  REQUIRE( l_srcCrds[0][1] == Approx( 3795526.32237143 ) );

  // get source coordinats of #2 and #3
  l_nrf1.getSrcCrds( 1, l_srcCrds, 1, 2 );

  REQUIRE( l_srcCrds[0][0] == Approx( 913609.593554401 ) );
  REQUIRE( l_srcCrds[0][1] == Approx( 3795526.32237143 ) );

  REQUIRE( l_srcCrds[1][0] == Approx( 913593.932784665 ) );
  REQUIRE( l_srcCrds[1][1] == Approx( 3795625.7327648  ) );
}

TEST_CASE( "3D Nrf: srcCrds", "[nrf3d][srcCrds]" ) {
  std::string l_ncFile = "cont/unit_tests/elastic/sources/nrf_3d.nc";

  // return silently if file does not exist
  // we assume that the unit tests executable is simply run from a different location
  if( !std::ifstream(l_ncFile) ) return;

  // NRF reader with 2-source buffer
  edge::elastic::io::Nrf< 3 > l_nrf1( 2 );

  // init
  l_nrf1.init( l_ncFile );
  l_nrf1.init( l_ncFile );

  // check basic quantities
  REQUIRE( l_nrf1.nIn() == 2 );
  REQUIRE( l_nrf1.nSrcsG( 1 ) == 3 );

  // source coordinates
  double l_srcCrds[3][3];

  // get source coordinates without arguments
  l_nrf1.getSrcCrds( 1, l_srcCrds );

  REQUIRE( l_srcCrds[0][0] == Approx( 913609.762166668 ) );
  REQUIRE( l_srcCrds[0][1] == Approx( 3795426.22805577 ) );
  REQUIRE( l_srcCrds[0][2] == Approx( -50.0            ) );

  REQUIRE( l_srcCrds[1][0] == Approx( 913609.593554401 ) );
  REQUIRE( l_srcCrds[1][1] == Approx( 3795526.32237143 ) );
  REQUIRE( l_srcCrds[1][2] == Approx( -50.0            ) );

  REQUIRE( l_srcCrds[2][0] == Approx( 913593.932784665 ) );
  REQUIRE( l_srcCrds[2][1] == Approx( 3795625.7327648  ) );
  REQUIRE( l_srcCrds[2][2] == Approx( -50.0            ) );

  // get source coordinates of #2 only
  l_nrf1.getSrcCrds( 1, l_srcCrds, 1, 1 );
  REQUIRE( l_srcCrds[0][0] == Approx( 913609.593554401 ) );
  REQUIRE( l_srcCrds[0][1] == Approx( 3795526.32237143 ) );
  REQUIRE( l_srcCrds[0][2] == Approx( -50.0            ) );

  // get source coordinats of #2 and #3
  l_nrf1.getSrcCrds( 1, l_srcCrds, 1, 2 );

  REQUIRE( l_srcCrds[0][0] == Approx( 913609.593554401 ) );
  REQUIRE( l_srcCrds[0][1] == Approx( 3795526.32237143 ) );
  REQUIRE( l_srcCrds[0][2] == Approx( -50.0            ) );

  REQUIRE( l_srcCrds[1][0] == Approx( 913593.932784665 ) );
  REQUIRE( l_srcCrds[1][1] == Approx( 3795625.7327648  ) );
  REQUIRE( l_srcCrds[1][2] == Approx( -50.0            ) );
}

TEST_CASE( "2D Nrf: getOffSetsG", "[nrf2d][getOffSetsG]" ) {
  std::string l_ncFile = "cont/unit_tests/elastic/sources/nrf_2d.nc";

  // return silently if file does not exist
  if( !std::ifstream(l_ncFile) ) return;

  // NRF reader with 2-source buffer
  edge::elastic::io::Nrf< 2 > l_nrf1( 2 );

  // init
  l_nrf1.init( l_ncFile );
  l_nrf1.init( l_ncFile );

  // slip rate offsets
  unsigned int l_off[4];

  // all offsets, normal, first kinematic source
  l_nrf1.getOffSetsG( 0,
                      0,
                      l_off );
  for( unsigned short l_di = 0; l_di < 2; l_di++ ) REQUIRE( l_off[l_di] ==  0 );

  // all offsets, normal, second kinematic source
  l_nrf1.getOffSetsG( 1,
                      0,
                      l_off );
  for( unsigned short l_di = 0; l_di < 2; l_di++ ) REQUIRE( l_off[l_di] ==  0 );

  // all offsets, first along-fault, first kinematic source
  l_nrf1.getOffSetsG( 0,
                      1,
                      l_off );
  REQUIRE( l_off[0] ==  0 );
  REQUIRE( l_off[1] == 11 );
  REQUIRE( l_off[2] == 21 );

  // all offsets, first along-fault, second kinematic source
  l_nrf1.getOffSetsG( 1,
                      1,
                      l_off );
  REQUIRE( l_off[0] ==  0 );
  REQUIRE( l_off[1] == 11 );
  REQUIRE( l_off[2] == 21 );

  // second offset only
  l_off[1] = 33; l_off[2] = 33;

  // normal, first kinematic source
  l_nrf1.getOffSetsG( 0,
                      0,
                      l_off,
                      1, 1 );
  REQUIRE( l_off[0] == 0 );
  REQUIRE( l_off[1] == 33 );
  REQUIRE( l_off[2] == 33 );

  // first along-fault, first kinematic source
  l_nrf1.getOffSetsG( 0,
                      1,
                      l_off,
                      1, 1 );
  REQUIRE( l_off[0] == 11 );
  REQUIRE( l_off[1] == 33 );
  REQUIRE( l_off[2] == 33 );
}

TEST_CASE( "3D Nrf: getOffSetsG", "[nrf3d][getOffSetsG]" ) {
  std::string l_ncFile = "cont/unit_tests/elastic/sources/nrf_3d.nc";

  // return silently if file does not exist
  // we assume that the unit tests executable is simply run from a different location
  if( !std::ifstream(l_ncFile) ) return;

  // NRF reader with 2-source buffer
  edge::elastic::io::Nrf< 3 > l_nrf1( 2 );

  // init
  l_nrf1.init( l_ncFile );
  l_nrf1.init( l_ncFile );

  // slip rate offsets
  unsigned int l_off[4];

  // all offsets, normal, first kinematic source
  l_nrf1.getOffSetsG( 0,
                      0,
                      l_off );
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) REQUIRE( l_off[l_di] ==  0 );

  // all offsets, normal, second kinematic source
  l_nrf1.getOffSetsG( 1,
                      0,
                      l_off );
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) REQUIRE( l_off[l_di] ==  0 );

  // all offsets, first along-fault, first kinematic source
  l_nrf1.getOffSetsG( 0,
                      1,
                      l_off );
  REQUIRE( l_off[0] ==  0 );
  REQUIRE( l_off[1] == 11 );
  REQUIRE( l_off[2] == 21 );

  // all offsets, first along-fault, second kinematic source
  l_nrf1.getOffSetsG( 1,
                      1,
                      l_off );
  REQUIRE( l_off[0] ==  0 );
  REQUIRE( l_off[1] == 11 );
  REQUIRE( l_off[2] == 21 );

  // all offsets, second along-fault, first kinematic source
  l_nrf1.getOffSetsG( 0,
                      2,
                      l_off );
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) REQUIRE( l_off[l_di] ==  0 );

  // all offsets, second along-fault, second kinematic source
  l_nrf1.getOffSetsG( 1,
                      2,
                      l_off );
  for( unsigned short l_di = 0; l_di < 3; l_di++ ) REQUIRE( l_off[l_di] ==  0 );

  // second offset only
  l_off[1] = 33; l_off[2] = 33;

  // normal, first kinematic source
  l_nrf1.getOffSetsG( 0,
                      0,
                      l_off,
                      1, 1 );
  REQUIRE( l_off[0] == 0 );
  REQUIRE( l_off[1] == 33 );
  REQUIRE( l_off[2] == 33 );

  // first along-fault, first kinematic source
  l_nrf1.getOffSetsG( 0,
                      1,
                      l_off,
                      1, 1 );
  REQUIRE( l_off[0] == 11 );
  REQUIRE( l_off[1] == 33 );
  REQUIRE( l_off[2] == 33 );
}

TEST_CASE( "2D Nrf: Gets of source data", "[nrf2d][gets]" ) {
  std::string l_ncFile = "cont/unit_tests/elastic/sources/nrf_2d.nc";

  // return silently if file does not exist
  if( !std::ifstream(l_ncFile) ) return;

  // NRF reader with 2-source buffer
  edge::elastic::io::Nrf< 2 > l_nrf1( 2 );

  // init
  l_nrf1.init( l_ncFile );
  l_nrf1.init( l_ncFile );
  l_nrf1.init( l_ncFile );

  for( unsigned short l_is = 0; l_is < 3; l_is++ ) {
    REQUIRE( l_nrf1.getOnSet(l_is, 0) == Approx(10.5201) );
    REQUIRE( l_nrf1.getOnSet(l_is, 1) == Approx(10.5255) );
    REQUIRE( l_nrf1.getOnSet(l_is, 2) == Approx(10.4184) );

    REQUIRE( l_nrf1.getDt(l_is, 0) == Approx(0.1) );
    REQUIRE( l_nrf1.getDt(l_is, 1) == Approx(0.2) );
    REQUIRE( l_nrf1.getDt(l_is, 2) == Approx(0.3) );

    REQUIRE( l_nrf1.getMu(l_is, 0) == Approx(3168000000) );
    REQUIRE( l_nrf1.getMu(l_is, 1) == Approx(3168111111) );
    REQUIRE( l_nrf1.getMu(l_is, 2) == Approx(3168222222) );

    REQUIRE( l_nrf1.getA(l_is, 0) == Approx(12222) );
    REQUIRE( l_nrf1.getA(l_is, 1) == Approx(11111) );
    REQUIRE( l_nrf1.getA(l_is, 2) == Approx(10000) );

    // slip-diretions
    double l_sds[2][2];

    l_nrf1.getSds(l_is, 0, l_sds);
    REQUIRE( l_sds[0][0] == Approx(  0.911320106822418  ) );
    REQUIRE( l_sds[0][1] == Approx( -0.063725709730111  ) );

    REQUIRE( l_sds[1][0] == Approx( -0.0838743060458509 ) );
    REQUIRE( l_sds[1][1] == Approx( -0.995966174736174  ) );


    l_nrf1.getSds(l_is, 1, l_sds);
    REQUIRE( l_sds[0][0] == Approx(  0.912988950432977  ) );
    REQUIRE( l_sds[0][1]== Approx(  -0.0318822766865827 ) );

    REQUIRE( l_sds[1][0] == Approx( -0.0911321737575558 ) );
    REQUIRE( l_sds[1][1] == Approx( -0.987689275364108  ) );

    l_nrf1.getSds(l_is, 2, l_sds);
    REQUIRE( l_sds[0][0] == Approx(  0.916427924584317  ) );
    REQUIRE( l_sds[0][1] == Approx(  0.178135543283137  ) );

    REQUIRE( l_sds[1][0] == Approx(  0.100036968135533 ) );
    REQUIRE( l_sds[1][1] == Approx( -0.969011238609667 ) );

    // slip-rate samples
    double l_srs[15];
    for( unsigned short l_sr = 0; l_sr < 15; l_sr++ ) l_srs[l_sr] = -1.0;

    // inactive
    l_nrf1.getSrs(l_is, 0, 0, l_srs);
    for( unsigned short l_sr = 0; l_sr < 15; l_sr++ ) REQUIRE( l_srs[l_sr] == Approx(-1.0) );

    // active
    l_nrf1.getSrs(l_is, 0, 1, l_srs);
    REQUIRE( l_srs[ 0] == Approx( 0.0        ) );
    REQUIRE( l_srs[ 1] == Approx( 0.332642   ) );
    REQUIRE( l_srs[ 2] == Approx( 0.289227   ) );
    REQUIRE( l_srs[ 3] == Approx( 0.0753253  ) );
    REQUIRE( l_srs[ 4] == Approx( 0.0669723  ) );
    REQUIRE( l_srs[ 5] == Approx( 0.0557519  ) );
    REQUIRE( l_srs[ 6] == Approx( 0.0428712  ) );
    REQUIRE( l_srs[ 7] == Approx( 0.0297156  ) );
    REQUIRE( l_srs[ 8] == Approx( 0.0177005  ) );
    REQUIRE( l_srs[ 9] == Approx( 0.0081184  ) );
    REQUIRE( l_srs[10] == Approx( 0.00199999 ) );
    for( unsigned short l_sr = 11; l_sr < 15; l_sr++ ) REQUIRE( l_srs[l_sr] == Approx(-1.0) );

    for( unsigned short l_sr = 0; l_sr < 15; l_sr++ ) l_srs[l_sr] = -1.0;
    l_nrf1.getSrs(l_is, 1, 1, l_srs);
    REQUIRE( l_srs[ 0] == Approx( 0.0        ) );
    REQUIRE( l_srs[ 1] == Approx( 0.497319   ) );
    REQUIRE( l_srs[ 2] == Approx( 0.347524   ) );
    REQUIRE( l_srs[ 3] == Approx( 0.104901   ) );
    REQUIRE( l_srs[ 4] == Approx( 0.0912566  ) );
    REQUIRE( l_srs[ 5] == Approx( 0.0734732  ) );
    REQUIRE( l_srs[ 6] == Approx( 0.0537081  ) );
    REQUIRE( l_srs[ 7] == Approx( 0.0343584  ) );
    REQUIRE( l_srs[ 8] == Approx( 0.0177707  ) );
    REQUIRE( l_srs[ 9] == Approx( 0.00595686 ) );
    for( unsigned short l_sr = 10; l_sr < 15; l_sr++ ) REQUIRE( l_srs[l_sr] == Approx(-1.0) );

    for( unsigned short l_sr = 0; l_sr < 15; l_sr++ ) l_srs[l_sr] = -1.0;
    l_nrf1.getSrs(l_is, 2, 1, l_srs);
    REQUIRE( l_srs[ 0] == Approx( 0.0        ) );
    REQUIRE( l_srs[ 1] == Approx( 0.647666   ) );
    REQUIRE( l_srs[ 2] == Approx( 0.383135   ) );
    REQUIRE( l_srs[ 3] == Approx( 0.130604   ) );
    REQUIRE( l_srs[ 4] == Approx( 0.111769   ) );
    REQUIRE( l_srs[ 5] == Approx( 0.0877252  ) );
    REQUIRE( l_srs[ 6] == Approx( 0.0616282  ) );
    REQUIRE( l_srs[ 7] == Approx( 0.0369034  ) );
    REQUIRE( l_srs[ 8] == Approx( 0.0167959  ) );
    REQUIRE( l_srs[ 9] == Approx( 0.00394474 ) );
    for( unsigned short l_sr = 10; l_sr < 15; l_sr++ ) REQUIRE( l_srs[l_sr] == Approx(-1.0) );
  }
}

TEST_CASE( "3D Nrf: Gets of source data", "[nrf3d][gets]" ) {
  std::string l_ncFile = "cont/unit_tests/elastic/sources/nrf_3d.nc";

  // return silently if file does not exist
  // we assume that the unit tests executable is simply run from a different location
  if( !std::ifstream(l_ncFile) ) return;

  // NRF reader with 2-source buffer
  edge::elastic::io::Nrf< 3 > l_nrf1( 2 );

  // init
  l_nrf1.init( l_ncFile );
  l_nrf1.init( l_ncFile );
  l_nrf1.init( l_ncFile );

  for( unsigned short l_is = 0; l_is < 3; l_is++ ) {
    REQUIRE( l_nrf1.getOnSet(l_is, 0) == Approx(10.5201) );
    REQUIRE( l_nrf1.getOnSet(l_is, 1) == Approx(10.5255) );
    REQUIRE( l_nrf1.getOnSet(l_is, 2) == Approx(10.4184) );

    REQUIRE( l_nrf1.getDt(l_is, 0) == Approx(0.1) );
    REQUIRE( l_nrf1.getDt(l_is, 1) == Approx(0.2) );
    REQUIRE( l_nrf1.getDt(l_is, 2) == Approx(0.3) );

    REQUIRE( l_nrf1.getMu(l_is, 0) == Approx(3168000000) );
    REQUIRE( l_nrf1.getMu(l_is, 1) == Approx(3168111111) );
    REQUIRE( l_nrf1.getMu(l_is, 2) == Approx(3168222222) );

    REQUIRE( l_nrf1.getA(l_is, 0) == Approx(12222) );
    REQUIRE( l_nrf1.getA(l_is, 1) == Approx(11111) );
    REQUIRE( l_nrf1.getA(l_is, 2) == Approx(10000) );

    // slip-diretions
    double l_sds[3][3];

    l_nrf1.getSds(l_is, 0, l_sds);
    REQUIRE( l_sds[0][0] == Approx(  0.911320106822418  ) );
    REQUIRE( l_sds[0][1] == Approx( -0.063725709730111  ) );
    REQUIRE( l_sds[0][2] == Approx(  0.4067366430758    ) );

    REQUIRE( l_sds[1][0] == Approx( -0.0838743060458509 ) );
    REQUIRE( l_sds[1][1] == Approx( -0.995966174736174  ) );
    REQUIRE( l_sds[1][2] == Approx(  0.0318822766865828 ) );

    REQUIRE( l_sds[2][0] == Approx(  0.403064217819573  ) );
    REQUIRE( l_sds[2][1] == Approx( -0.0631697134771603 ) );
    REQUIRE( l_sds[2][2] == Approx( -0.912988950432977  ) );

    l_nrf1.getSds(l_is, 1, l_sds);
    REQUIRE( l_sds[0][0] == Approx(  0.912988950432977  ) );
    REQUIRE( l_sds[0][1]== Approx(  -0.0318822766865827 ) );
    REQUIRE( l_sds[0][2] == Approx(  0.4067366430758    ) );

    REQUIRE( l_sds[1][0] == Approx( -0.0911321737575558 ) );
    REQUIRE( l_sds[1][1] == Approx( -0.987689275364108  ) );
    REQUIRE( l_sds[1][2] == Approx(   0.127140954208103 ) );

    REQUIRE( l_sds[2][0] == Approx(  0.397675877183308  ) );
    REQUIRE( l_sds[2][1] == Approx( -0.153145080769852  ) );
    REQUIRE( l_sds[2][2] == Approx( -0.904654896047372  ) );

    l_nrf1.getSds(l_is, 2, l_sds);
    REQUIRE( l_sds[0][0] == Approx(  0.916427924584317  ) );
    REQUIRE( l_sds[0][1] == Approx(  0.178135543283137  ) );
    REQUIRE( l_sds[0][2] == Approx(  0.3583679495453    ) );

    REQUIRE( l_sds[1][0] == Approx(  0.100036968135533 ) );
    REQUIRE( l_sds[1][1] == Approx( -0.969011238609667 ) );
    REQUIRE( l_sds[1][2] == Approx(  0.225853546472949 ) );

    REQUIRE( l_sds[2][0] == Approx(  0.38749511487028  ) );
    REQUIRE( l_sds[2][1] == Approx( -0.171128453704753 ) );
    REQUIRE( l_sds[2][2] == Approx( -0.905849097965157 ) );

    // slip-rate samples
    double l_srs[15];
    for( unsigned short l_sr = 0; l_sr < 15; l_sr++ ) l_srs[l_sr] = -1.0;

    // inactive
    l_nrf1.getSrs(l_is, 0, 0, l_srs);
    for( unsigned short l_sr = 0; l_sr < 15; l_sr++ ) REQUIRE( l_srs[l_sr] == Approx(-1.0) );
    l_nrf1.getSrs(l_is, 0, 2, l_srs);
    for( unsigned short l_sr = 0; l_sr < 15; l_sr++ ) REQUIRE( l_srs[l_sr] == Approx(-1.0) );

    // active
    l_nrf1.getSrs(l_is, 0, 1, l_srs);
    REQUIRE( l_srs[ 0] == Approx( 0.0        ) );
    REQUIRE( l_srs[ 1] == Approx( 0.332642   ) );
    REQUIRE( l_srs[ 2] == Approx( 0.289227   ) );
    REQUIRE( l_srs[ 3] == Approx( 0.0753253  ) );
    REQUIRE( l_srs[ 4] == Approx( 0.0669723  ) );
    REQUIRE( l_srs[ 5] == Approx( 0.0557519  ) );
    REQUIRE( l_srs[ 6] == Approx( 0.0428712  ) );
    REQUIRE( l_srs[ 7] == Approx( 0.0297156  ) );
    REQUIRE( l_srs[ 8] == Approx( 0.0177005  ) );
    REQUIRE( l_srs[ 9] == Approx( 0.0081184  ) );
    REQUIRE( l_srs[10] == Approx( 0.00199999 ) );
    for( unsigned short l_sr = 11; l_sr < 15; l_sr++ ) REQUIRE( l_srs[l_sr] == Approx(-1.0) );

    for( unsigned short l_sr = 0; l_sr < 15; l_sr++ ) l_srs[l_sr] = -1.0;
    l_nrf1.getSrs(l_is, 1, 1, l_srs);
    REQUIRE( l_srs[ 0] == Approx( 0.0        ) );
    REQUIRE( l_srs[ 1] == Approx( 0.497319   ) );
    REQUIRE( l_srs[ 2] == Approx( 0.347524   ) );
    REQUIRE( l_srs[ 3] == Approx( 0.104901   ) );
    REQUIRE( l_srs[ 4] == Approx( 0.0912566  ) );
    REQUIRE( l_srs[ 5] == Approx( 0.0734732  ) );
    REQUIRE( l_srs[ 6] == Approx( 0.0537081  ) );
    REQUIRE( l_srs[ 7] == Approx( 0.0343584  ) );
    REQUIRE( l_srs[ 8] == Approx( 0.0177707  ) );
    REQUIRE( l_srs[ 9] == Approx( 0.00595686 ) );
    for( unsigned short l_sr = 10; l_sr < 15; l_sr++ ) REQUIRE( l_srs[l_sr] == Approx(-1.0) );

    for( unsigned short l_sr = 0; l_sr < 15; l_sr++ ) l_srs[l_sr] = -1.0;
    l_nrf1.getSrs(l_is, 2, 1, l_srs);
    REQUIRE( l_srs[ 0] == Approx( 0.0        ) );
    REQUIRE( l_srs[ 1] == Approx( 0.647666   ) );
    REQUIRE( l_srs[ 2] == Approx( 0.383135   ) );
    REQUIRE( l_srs[ 3] == Approx( 0.130604   ) );
    REQUIRE( l_srs[ 4] == Approx( 0.111769   ) );
    REQUIRE( l_srs[ 5] == Approx( 0.0877252  ) );
    REQUIRE( l_srs[ 6] == Approx( 0.0616282  ) );
    REQUIRE( l_srs[ 7] == Approx( 0.0369034  ) );
    REQUIRE( l_srs[ 8] == Approx( 0.0167959  ) );
    REQUIRE( l_srs[ 9] == Approx( 0.00394474 ) );
    for( unsigned short l_sr = 10; l_sr < 15; l_sr++ ) REQUIRE( l_srs[l_sr] == Approx(-1.0) );
  }
}
