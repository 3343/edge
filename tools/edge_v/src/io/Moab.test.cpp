/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019-2020, Alexander Breuer
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
 * Tests the MOAB interface.
 **/
#include <catch.hpp>
#define private public
#include "Moab.h"
#undef private

namespace edge_v {
  namespace test {
    extern std::string g_files;
  }
}

TEST_CASE( "Tests element type query.", "[moab][getElTy]" ) {
  // only continue if the unit test files are available
  if( edge_v::test::g_files != "" ) {
    // path to the mesh file
    std::string l_path = edge_v::test::g_files + "/tria3.msh";

    // construct the moab-reader
    edge_v::io::Moab l_moab( l_path );

    REQUIRE( l_moab.getElType() == edge_v::t_entityType::TRIA3 );
  }
}

TEST_CASE( "Tests the derivation of elements adjacent to their faces.", "[moab][getFaEl]" ) {
  // only continue if the unit test files are available
  if( edge_v::test::g_files != "" ) {
    // path to the mesh file
    std::string l_path = edge_v::test::g_files + "/tria3.msh";

    // construct the moab-reader
    edge_v::io::Moab l_moab( l_path );

    edge_v::t_idx l_faEl[48 * 2];
    l_moab.getFaEl( edge_v::TRIA3,
                    l_faEl );

    REQUIRE( l_faEl[0*2+0] == 5 );
    REQUIRE( l_faEl[0*2+1] == std::numeric_limits< edge_v::t_idx > ::max() );

    REQUIRE( l_faEl[12*2+0] == 4 );
    REQUIRE( l_faEl[12*2+1] == 5 );

    REQUIRE( l_faEl[15*2+0] == 1 );
    REQUIRE( l_faEl[15*2+1] == 3 );
  }
}

TEST_CASE( "Tests the derivation of faces adjacent to elements.", "[moab][getElFa]" ) {
  // only continue if the unit test files are available
  if( edge_v::test::g_files != "" ) {
    // path to the mesh file
    std::string l_path = edge_v::test::g_files + "/tria3.msh";

    // construct the moab-reader
    edge_v::io::Moab l_moab( l_path );

    edge_v::t_idx l_elFa[28 * 3];
    l_moab.getElFa( edge_v::TRIA3,
                    l_elFa );

    REQUIRE( l_elFa[0*3 + 0] ==  9 );

    for( unsigned short l_el = 0; l_el < 23; l_el++ ) {
      for( unsigned short l_fa = 0; l_fa < 3-1; l_fa++ ) {
        REQUIRE( l_elFa[l_el*3 + l_fa] <  l_elFa[l_el*3 + l_fa+1] );
      }
    }
  }
}

TEST_CASE( "Set and get data for entities.", "[moab][setGetEnData]" ) {
  // only continue if the unit test files are available
  if( edge_v::test::g_files != "" ) {
    // path to the mesh file
    std::string l_path = edge_v::test::g_files + "/tria3.msh";

    // construct the moab-reader
    edge_v::io::Moab l_moab( l_path );

    // set data
    edge_v::t_idx l_nTria3 = l_moab.nEnsByType( edge_v::TRIA3 );
    unsigned short * l_setData = new unsigned short[l_nTria3];
    for( edge_v::t_idx l_en = 0; l_en < l_nTria3; l_en++ ) {
      l_setData[l_en] = l_en % 25 / 2;
    }

    std::string l_tagName = "some_tag_name";
    l_moab.setEnData( edge_v::TRIA3,
                      l_tagName,
                      l_setData );


    // get data
    unsigned short * l_getData = new unsigned short[l_nTria3];

    l_moab.getEnData( edge_v::TRIA3,
                      l_tagName,
                      l_getData );

    // check the results
    for( edge_v::t_idx l_en = 0; l_en < l_nTria3; l_en++ ) {
      REQUIRE( l_getData[l_en] == l_en % 25 / 2 );
    }

    delete[] l_setData;
    delete[] l_getData;
  }
}

TEST_CASE( "Gets the material set data for faces.", "[moab][getMatSetFa]" ) {
  // only continue if the unit test files are available
  if( edge_v::test::g_files != "" ) {
    // path to the mesh file
    std::string l_path = edge_v::test::g_files + "/tria3.msh";

    // construct the moab-reader
    edge_v::io::Moab l_moab( l_path );

    // check that the MATERIAL_SET tag is present
    std::vector< std::string > l_tagNames;
    l_moab.getTagNames( l_tagNames );

    bool l_foundMs = false;
    for( unsigned short l_tn = 0; l_tn < l_tagNames.size(); l_tn++ ) {
      if( l_tagNames[l_tn] == "MATERIAL_SET" ) {
        l_foundMs = true;
      }
    }
    REQUIRE( l_foundMs );

    // get the material set values
    edge_v::t_idx l_nLine = l_moab.nEnsByType( edge_v::LINE );
    int * l_matSetData = new int[l_nLine];

    l_moab.getEnDataFromSet( edge_v::LINE,
                             "MATERIAL_SET",
                             l_matSetData );

    // check the results
    for( unsigned short l_li = 0; l_li < l_nLine; l_li++ ) {
      if( l_li < 12 ) {
        REQUIRE( l_matSetData[l_li] == 106 );
      }
      else {
        REQUIRE( l_matSetData[l_li] == std::numeric_limits< int >::max() );
      }

    }

    delete[] l_matSetData;
  }
}

TEST_CASE( "Gets the material set data for elements.", "[moab][getMatSetEl]" ) {
  // only continue if the unit test files are available
  if( edge_v::test::g_files != "" ) {
    // path to the mesh file
    std::string l_path = edge_v::test::g_files + "/tria3.msh";

    // construct the moab-reader
    edge_v::io::Moab l_moab( l_path );

    // get the material set values
    edge_v::t_idx l_nTria3 = l_moab.nEnsByType( edge_v::TRIA3 );
    int * l_matSetData = new int[l_nTria3];

    l_moab.getEnDataFromSet( edge_v::TRIA3,
                             "MATERIAL_SET",
                             l_matSetData );

    // check the results
    for( unsigned short l_tr = 0; l_tr < l_nTria3; l_tr++ ) {
      REQUIRE( l_matSetData[l_tr] == 1 );
    }

    delete[] l_matSetData;
  }
}

TEST_CASE( "Mesh reordering.", "[moab][reorder]" ) {
  // only continue if the unit test files are available
  if( edge_v::test::g_files != "" ) {
    // path to the mesh file
    std::string l_path = edge_v::test::g_files + "/tria3.msh";

    // construct the moab-reader
    edge_v::io::Moab l_moab( l_path );

    REQUIRE( l_moab.getElType() == edge_v::t_entityType::TRIA3 );

    // get the vertices of the triangles before the reordering
    edge_v::t_idx l_nTria3 = l_moab.nEnsByType( edge_v::TRIA3 );
    edge_v::t_idx *l_elVe0 = new edge_v::t_idx[l_nTria3*3];
    l_moab.getEnVe( edge_v::TRIA3,
                    l_elVe0 );

    // specify priorities (lower is higher)
    unsigned short * l_prio = new unsigned short[l_nTria3];
    REQUIRE( l_nTria3 > 10 );
    l_prio[0] = 2;
    l_prio[1] = 2;
    l_prio[2] = 1;
    l_prio[3] = 3;
    l_prio[4] = 0;
    l_prio[5] = 2;
    l_prio[6] = 3;
    l_prio[7] = 9;
    for( edge_v::t_idx l_el = 8; l_el < l_nTria3; l_el++ ) {
      l_prio[l_el] = 4;
    }

    // assign tag
    std::string l_tagName = "priority";
    l_moab.setEnData( edge_v::TRIA3,
                      l_tagName,
                      l_prio );

    // perform the reordering
    l_moab.reorder( edge_v::TRIA3,
                    l_tagName );

    // get the vertices of the triangles after the reordering
    edge_v::t_idx *l_elVe1 = new edge_v::t_idx[l_nTria3*3];
    l_moab.getEnVe( edge_v::TRIA3,
                    l_elVe1 );

    for( unsigned short l_ve = 0; l_ve < 3; l_ve++ ) {
      REQUIRE( l_elVe0[4*3 + l_ve] == l_elVe1[0*3 + l_ve] );
      REQUIRE( l_elVe0[2*3 + l_ve] == l_elVe1[1*3 + l_ve] );
      REQUIRE( l_elVe0[0*3 + l_ve] == l_elVe1[2*3 + l_ve] );
      REQUIRE( l_elVe0[1*3 + l_ve] == l_elVe1[3*3 + l_ve] );
      REQUIRE( l_elVe0[5*3 + l_ve] == l_elVe1[4*3 + l_ve] );
      REQUIRE( l_elVe0[3*3 + l_ve] == l_elVe1[5*3 + l_ve] );
      REQUIRE( l_elVe0[6*3 + l_ve] == l_elVe1[6*3 + l_ve] );
      for( edge_v::t_idx l_el = 7; l_el < l_nTria3; l_el++ ) {
        REQUIRE( l_elVe0[(7+1)*3 + l_ve] == l_elVe1[7*3 + l_ve] );
      }
       REQUIRE( l_elVe0[7*3 + l_ve] == l_elVe1[(l_nTria3-1)*3 + l_ve] );
    }

    // free memory
    delete[] l_elVe0;
    delete[] l_elVe1;
    delete[] l_prio;
  }
}