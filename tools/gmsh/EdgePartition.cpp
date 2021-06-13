/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section LICENSE
 * Copyright (c) 2021, Friedrich Schiller University Jena
 * Copyright (c) 2020, Alexander Breuer
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
 * Gmsh plugin which partitions meshes for EDGE.
 **/
#include "EdgePartition.h"
#include "gmsh.h"
#include "GModel.h"
#include "meshPartition.h"
#include "Context.h"
#include <cassert>
#include <limits>

StringXNumber EdgePartitionOptions_Number[] = {
  { GMSH_FULLRC,
    "nPaEls",
    NULL,
    0. }
};

extern "C" {
  GMSH_Plugin *GMSH_RegisterEdgePartitionPlugin() {
    return new GMSH_EdgePartitionPlugin();
  }
}

std::string GMSH_EdgePartitionPlugin::getName() const {
  return "EdgePartition";
}

std::string GMSH_EdgePartitionPlugin::getShortHelp() const {
  return "Support for EDGE's custom partitioning.";
}

std::string GMSH_EdgePartitionPlugin::getHelp() const {
  return "Plugin(EdgePartition) incorporates custom, externally computed partitioning for the solver EDGE (available from https://dial3343.org) into Gmsh.";
}

int GMSH_EdgePartitionPlugin::getNbOptions() const {
  return 1;
}

StringXNumber* GMSH_EdgePartitionPlugin::getOption( int i_opt ) {
  return &EdgePartitionOptions_Number[i_opt];
}

void GMSH_EdgePartitionPlugin::run() {
  // get the number of dimensions
  int l_nDis = gmsh::model::getDimension();
  assert( l_nDis >= 2 );

  // query gmsh for the partitions of the elements given by EDGE
  std::vector< std::vector< double > > l_nPaEls;
  std::vector< int > l_nEls;
  std::vector< std::string > l_tys;
  double l_view = EdgePartitionOptions_Number[0].def;
  gmsh::view::getListData( l_view,
                           l_tys,
                           l_nEls,
                           l_nPaEls );
  assert( l_tys.size()  == 1 );
  assert( l_nEls.size() == 1 );
  assert( l_nPaEls.size() == 1 );
  assert( l_nPaEls[0].size() >= 1 );

  // get entities/regions of the Gmsh model
  GModel* l_model = GModel::current();
  std::vector< GEntity* > l_regions;
  l_model->getEntities( l_regions );

  // assemble a list of explicitly defined faces, consisting of their vertex ids and pointers
  Msg::Info("Assembling list of explicitly stored faces");
  struct Face {
    MElement* ptr;
    std::size_t ves[4] = { std::numeric_limits< std::size_t >::max(),
                           std::numeric_limits< std::size_t >::max(),
                           std::numeric_limits< std::size_t >::max(),
                           std::numeric_limits< std::size_t >::max() };
  };
  std::vector< Face > l_mapFas;

  for( std::size_t l_re = 0; l_re < l_regions.size(); l_re++ ) {
    GEntity* l_region = l_regions[l_re];

    for( std::size_t l_en = 0; l_en < l_region->getNumMeshElements(); l_en++ ) {
      MElement* l_ptr = l_region->getMeshElement(l_en);

      if( l_ptr->getDim() == l_nDis-1 ) {
        l_mapFas.resize( l_mapFas.size()+1 );
        l_mapFas.back().ptr = l_ptr;

        std::vector< MVertex* > l_enVes;
        l_ptr->getVertices( l_enVes );
        assert( l_enVes.size() < 4 );

        // store the (sorted) ids of the vertices
        for( std::size_t l_ve = 0; l_ve < l_enVes.size(); l_ve++ ) {
          l_mapFas.back().ves[ l_ve ] = l_enVes[l_ve]->getNum();
        }
        std::sort( l_mapFas.back().ves,
                   l_mapFas.back().ves+4 );
      }
    }
  }

  // sort the faces by their vertices
  Msg::Info("Sorting the faces by their vertex ids");
  std::sort( l_mapFas.begin(),
             l_mapFas.end(),
             []( Face const & l_f0, Face const & l_f1 ) -> bool {
               return   std::tie( l_f0.ves[0], l_f0.ves[1], l_f0.ves[2], l_f0.ves[3] )
                      < std::tie( l_f1.ves[0], l_f1.ves[1], l_f1.ves[2], l_f1.ves[3] ); } );


  // derive a mapping from the elements and the explicitly defined faces to the partitions
  Msg::Info("Deriving mapping from elements to partitions and faces to partitions");
  std::vector< std::pair< MElement *, int > > l_enPa;
  std::size_t l_el = 0;
  std::size_t l_pa = 0;
  std::size_t l_paNext = l_nPaEls[0][0];

  for( std::size_t l_re = 0; l_re < l_regions.size(); l_re++ ) {
    GEntity* l_region = l_regions[l_re];

    for( std::size_t l_en = 0; l_en < l_region->getNumMeshElements(); l_en++ ) {
      MElement* l_entity = l_region->getMeshElement( l_en );

      if( l_entity->getDim() == l_nDis ) {
        // increase partition counter based on element-id if required (stored ascending by assumption)
        if( l_el == l_paNext ) {
          l_pa++;
          if( l_pa < l_nPaEls[0].size() ) {
            l_paNext += l_nPaEls[0][l_pa];
          }
          else {
            l_paNext = std::numeric_limits< std::size_t >::max();
          }
        }
        l_el++;

        // store element to partition mapping
        l_enPa.push_back( std::pair< MElement*, int >( l_entity, l_pa+1) );

        // assign faces of the elements to respective partitions if explicitly stored
        // this assumes a unique mapping (no internal faces at partition boundaries)
        unsigned short l_nEf = 0;
        if( l_nDis == 3 ) {
          l_nEf = l_entity->getNumFaces();
        }
        else {
          l_nEf = l_entity->getNumEdges();
        }
        for( unsigned short l_ef = 0; l_ef < l_nEf; l_ef++ ) {
          std::vector< MVertex* > l_faVes;
          if( l_nDis == 3 ) {
            l_entity->getFaceVertices( l_ef,
                                       l_faVes );
          }
          else {
            l_entity->getEdgeVertices( l_ef,
                                       l_faVes );
          }
          assert( l_faVes.size() <= 4 );

          Face l_fa;
          for( std::size_t l_ve = 0; l_ve < l_faVes.size(); l_ve++ ) {
            l_fa.ves[l_ve] = l_faVes[l_ve]->getNum();
          }
          // sort vertices of the face
          std::sort( l_fa.ves,
                     l_fa.ves+4 );

          auto l_lb = std::lower_bound( l_mapFas.begin(),
                                        l_mapFas.end(),
                                        l_fa,
                                        []( Face const & l_f0, Face const & l_f1 ) -> bool {
                                          return   std::tie( l_f0.ves[0], l_f0.ves[1], l_f0.ves[2], l_f0.ves[3] )
                                                 < std::tie( l_f1.ves[0], l_f1.ves[1], l_f1.ves[2], l_f1.ves[3] ); } );

          // only save mapping if the faces is part of explicitly stored faces 
          if( l_lb != l_mapFas.end() ) {
            // only save mapping if the faces is present (no a "lower" bound)
            if(    std::tie(    l_fa.ves[0],    l_fa.ves[1],    l_fa.ves[2],    l_fa.ves[3] )
                == std::tie( (*l_lb).ves[0], (*l_lb).ves[1], (*l_lb).ves[2], (*l_lb).ves[3] ) ) {
              l_enPa.push_back( std::pair< MElement *, int >( (*l_lb).ptr, l_pa+1) );
            }
          }
        }
      }
    }
  }

  // partition the mesh
  Msg::Info("Partitioning the mesh through external functions");
  int l_err = PartitionUsingThisSplit( l_model,
                                       l_nPaEls[0].size(),
                                       l_enPa );
  if( !l_err ) {
    opt_mesh_color_carousel(0, GMSH_SET | GMSH_GUI, 3.);
    CTX::instance()->mesh.changed = ENT_ALL;
  }
  else {
    assert( false );
  }
}