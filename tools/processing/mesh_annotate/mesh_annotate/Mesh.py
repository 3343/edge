##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2018, Regents of the University of California
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
# Mesh interface.
##
import pymoab.core
import pymoab.types
import pymoab.rng
import numpy

class Mesh:
  # MOAB core
  m_moab = 0
  # root meshset
  m_root = 0

  m_types = {
    "vertex": pymoab.types.MBVERTEX,
    "edge":   pymoab.types.MBEDGE,
    "quad" :  pymoab.types.MBQUAD,
    "tria":   pymoab.types.MBTRI,
    "tet":    pymoab.types.MBTET,
    "hex":    pymoab.types.MBHEX
  }

  ##
  # Sets the tag "global_id" for the given entities in ascending order (first: 0).
  #
  # @param i_enHas entity handles.
  ##
  def setIds( self, i_enHas ):
    l_gIdTag = self.m_moab.tag_get_handle( "global_id",
                                           1,
                                           pymoab.types.MB_TYPE_INTEGER,
                                           pymoab.types.MB_TAG_DENSE,
                                           True )

    l_data = numpy.array( range(len(i_enHas)), dtype = 'int')

    self.m_moab.tag_set_data( l_gIdTag,
                              i_enHas,
                              l_data )

  ##
  # Wipes the tag "global_id".
  ##
  def deleteIds( self ):
    l_gidTag = self.m_moab.tag_get_handle( "global_id" )
    self.m_moab.tag_delete( l_gidTag )

  ##
  # Gets the tag "global_id" for the given entities.
  #
  # @param i_enHas entity handles.
  ##
  def getIds( self, i_enHas ):
    l_gIdTag = self.m_moab.tag_get_handle( "global_id" )
    l_gIds = self.m_moab.tag_get_data( l_gIdTag, i_enHas )

    # convert to plain list
    l_gIds = [ l_en[0] for l_en in l_gIds ]

    return l_gIds

  ##
  # Reads the given mesh.
  ##
  def __init__( self, i_path ):
    self.m_moab = pymoab.core.Core()

    self.m_moab.load_file( i_path )
    self.m_root = self.m_moab.get_root_set()

    # determine number of dimensions
    for l_di in [3,2]:
      if( len(self.m_moab.get_entities_by_dimension( self.m_root, l_di )) > 0 ):
        self.m_nDis = l_di

    l_ves = self.m_moab.get_entities_by_dimension( self.m_root, 0 )
    l_els = self.m_moab.get_entities_by_dimension( self.m_root, self.m_nDis )
    # create if not existent and get faces
    l_fas = self.m_moab.get_adjacencies( l_els, self.m_nDis-1, True, pymoab.types.UNION )

    # set ids
    for l_en in [l_ves, l_els, l_fas]:
      self.setIds( l_en )

  ##
  # Gets the number of entities.
  #
  # @param i_type type for which the number of entities is determined.
  # @return number of entities.
  ##
  def nEn( self, i_type ):
    l_moabType = self.m_types[i_type]
    l_ens = self.m_moab.get_entities_by_type( self.m_root, l_moabType )
    return len(l_ens)

  ##
  # Gets the connectivity information: entity -> vertex.
  #
  # @param i_ens entities to get the connectivity information for.
  # @return flat array with vertex ids.
  ##
  def conn( self,
            i_ens ):
    l_ves = self.m_moab.get_connectivity( i_ens )
    l_ves = self.getIds( l_ves )

    return l_ves

  ##
  # Determines to-entities adjacent to the from-entities (no bridge).
  #
  # @param i_nDisFrom from dimension.
  # @param i_nDisTo to dimension.
  # @return flat array with vertex ids
  ##
  def enAd( self,
            i_nDisFrom,
            i_nDisTo ):
    l_ens = self.m_moab.get_entities_by_dimension( self.m_root, i_nDisFrom )

    if( i_nDisTo == 0):
      l_enAd = self.conn( l_ens )
    else:
      l_enAd = []
      for l_en in l_ens:
        l_enAdHan = self.m_moab.get_adjacencies( l_en, i_nDisTo, op_type=pymoab.types.UNION )
        l_enAd = l_enAd + self.getIds( l_enAdHan )
    return l_enAd

  ##
  # Determines to-entities, adjacent to the from-entities, using a bridge.
  #
  # @param i_nDisFrom from dimension.
  # @param i_nDisTo to dimension.
  # @param i_nDisBridge bridge dimension.
  # @return adjacency information.
  ##
  def enAdBr( self,
              i_nDisFrom,
              i_nDisTo,
              i_nDisBridge ):

    l_ens = self.m_moab.get_entities_by_dimension( self.m_root, i_nDisFrom )

    l_enAdBr = []
    for l_en in l_ens:
      l_enBr = self.m_moab.get_adjacencies( l_en, i_nDisBridge, op_type=pymoab.types.UNION )

      l_enAd = []
      for l_br in l_enBr:
        l_enTo = self.m_moab.get_adjacencies( l_br, i_nDisTo, op_type=pymoab.types.UNION )
        l_enTo.erase( l_en )

        if( len(l_enTo) > 0 ):
          l_enAd = l_enAd + self.getIds( l_enTo )
        else:
          l_enAd = l_enAd + [-1]

      l_enAdBr = l_enAdBr + [l_enAd]

    return l_enAdBr

  ##
  # Gets the vertex coordinates
  #
  # @return vertex coordinates.
  ##
  def veCrds( self ):
    l_ves = self.m_moab.get_entities_by_dimension( self.m_root, 0 )
    l_crds = self.m_moab.get_coords( l_ves )

    # reshape
    l_crds = numpy.reshape(l_crds, (-1, 3))

    return l_crds

  ##
  # Gets the volume of the given entities.
  #
  # @param i_type type of the entities.
  # @return volumes of the entities.
  ##
  def vol( self,
           i_type ):
    l_moabType = self.m_types[i_type]

    l_veCrds = self.veCrds()
    l_ens = self.m_moab.get_entities_by_type( self.m_root, l_moabType )

    # get adjacent vertices
    l_enVe = []
    for l_en in l_ens:
      l_enVeHan = self.m_moab.get_adjacencies( l_en, 0, op_type=pymoab.types.UNION )
      l_enVe = l_enVe + [ self.getIds( l_enVeHan ) ]

    # compute volumes
    l_vol = []
    for l_en in l_enVe:
      if( i_type == 'line' ):
        l_t = l_veCrds[ l_en[1] ] - l_veCrds[ l_en[0] ]
        l_vol = l_vol + [ numpy.linalg.norm( l_t ) ]
      elif( i_type == 'tria3' ):
        l_mat = numpy.matrix( [ l_veCrds[ l_en[0] ],
                                l_veCrds[ l_en[1] ], 
                                l_veCrds[ l_en[2] ] ] )
        l_mat[:, 2 ] = 1

        l_vol = l_vol + [ 0.5 * abs( numpy.linalg.det( l_mat ) ) ]
      else: assert( False )

    return l_vol
  
  ##
  # Derives the centroids for the given entity type.
  #
  # @return list containing the centroids.
  ##
  def cen( self,
           i_type ):
    l_moabType = self.m_types[i_type]

    l_veCrds = self.veCrds()
    l_veCrds = numpy.array( l_veCrds )
    l_ens = self.m_moab.get_entities_by_type( self.m_root, l_moabType )

    # compute centroids
    l_cens = []
    for l_en in l_ens:
      l_enVeHan = self.m_moab.get_adjacencies( l_en, 0, op_type=pymoab.types.UNION )
      
      l_ves = l_veCrds[ self.getIds( l_enVeHan ) ]
      l_cen = numpy.zeros(3)
      for l_ve in l_ves:
        l_cen = l_cen + l_ve
      l_cen = l_cen / len(l_ves)

      l_cens = l_cens + [l_cen]

    return l_cens

  ##
  # Gets the face normals.
  #
  # @param i_faType type of the face.
  # @return face normals.
  ##
  def faNor( self,
             i_faType ):
    l_moabFaType = self.m_types[i_faType]

    if( i_faType == 'line' ): l_nDis = 2
    else:                     l_nDis = 3

    # get vertex coordinates
    l_veCrds = self.veCrds()

    # get faces
    l_fas = self.m_moab.get_entities_by_type( self.m_root, l_moabFaType )

    # determine normals
    l_nor = []

    for l_fa in l_fas:
      # get face's vertices
      l_faVeHan = self.m_moab.get_adjacencies( l_fa, 0, op_type=pymoab.types.UNION )
      l_faVe = self.getIds( l_faVeHan )

      # get coordinates
      l_faVeCrds = l_veCrds[ l_faVe ]
      
      # compute tangent
      if( i_faType == 'line' ):
        l_faTan = l_faVeCrds[1] - l_faVeCrds[0]

      # compute normal
      if( i_faType == 'line' ):
        l_faNor = numpy.array( [ l_faTan[1], -l_faTan[0], 0 ] )
        l_faNor = l_faNor / numpy.linalg.norm( l_faNor )

      # get adjacent elements
      l_faEl = self.m_moab.get_adjacencies( l_fa, l_nDis, op_type=pymoab.types.UNION )
      assert( len(l_faEl) > 0 )

      # get vertices of first element
      l_elVeHan = self.m_moab.get_adjacencies( l_faEl[0], 0, op_type=pymoab.types.UNION )
      assert( len(l_elVeHan) > 2 )

      # determine additional points, defining "left"
      l_elVeSubHan = pymoab.rng.subtract( l_elVeHan, l_faVeHan )
      assert( len(l_elVeSubHan) > 0 )
      l_elVeSub = self.getIds( l_elVeSubHan )

      # get normal point
      l_nPt = l_veCrds[ l_elVeSub[0] ]

      # determine right-direction
      l_dirR = l_nPt - l_faVeCrds[0]

      l_dp = numpy.dot( l_faNor, l_dirR )
      # change direction, if normal points right
      if( l_dp > 0 ):
        l_faNor = -l_faNor

      l_nor = l_nor + [l_faNor]

    return l_nor
  
  ##
  # Gets the types of the given entities.
  #
  # @param i_type type of the entities, e.g., 'line' or 'tria3'.
  # @return data type of the entities (as integers).
  ##
  def type( self,
            i_type ):
    l_moabType = self.m_types[i_type]

    # get material tag
    l_tagMat = self.m_moab.tag_get_handle( "MATERIAL_SET" );

    # get meshsets
    l_mss = self.m_moab.get_entities_by_type( self.m_root, pymoab.types.MBENTITYSET )

    # get entities of the meshsets (missing `contains_entities` in pympab)
    l_msEns = []
    for l_ms in l_mss:
      l_msEns = l_msEns + [ self.m_moab.get_entities_by_handle( l_ms ) ]

    # get data of the meshsets
    l_msData = self.m_moab.tag_get_data( l_tagMat, l_mss, flat=True )

    # get entities
    l_ens = self.m_moab.get_entities_by_type( self.m_root, l_moabType )

    # data of the entities
    l_enData = numpy.zeros( len(l_ens), dtype=int )

    # iterate over entities and set type
    for l_en in range( len(l_ens) ):
      for l_ms in range( len(l_mss) ):
        if( l_msData[l_ms] != -1 and l_ens[l_en] in l_msEns[l_ms] ):
          assert( l_enData[l_en] == 0 )
          l_enData[l_en] = l_msData[l_ms]
  
    return l_enData

  ##
  # Stores the given data for the entities.
  #
  # @param i_type type of the entities.
  # @param i_tag tag name (string).
  # @param i_data data, which will be stored for the entities.
  ##
  def store( self,
             i_type,
             i_tag,
             i_data ):

    # get tag handle
    l_tag = self.m_moab.tag_get_handle( i_tag,
                                        numpy.size(i_data[0]),
                                        pymoab.types.MB_TYPE_DOUBLE,
                                        pymoab.types.MB_TAG_DENSE,
                                        True )

    # get entities
    l_ens = self.m_moab.get_entities_by_type( self.m_root,
                                              self.m_types[i_type] )

    # set data
    self.m_moab.tag_set_data( l_tag, l_ens, i_data )

  ##
  # Writes the mesh to disk.
  #
  # @param i_path path to the output mesh.
  ##
  def write( self,
             i_path ):
    self.m_moab.write_file( i_path )