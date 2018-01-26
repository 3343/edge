#!/usr/bin/env python3
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2017, Regents of the University of California
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
# This is the entry point of EDGEpre.
##
import logging
import argparse
import os
import sympy
from sympy.core.cache import *
import edge_pre.io.Config
import edge_pre.int.Matrices
import edge_pre.sc.ops.Project
import edge_pre.io.ArrStr
import edge_pre.io.Matrices
import edge_pre.dg.basis.Mod
import edge_pre.dg.basis.Line
import edge_pre.sc.grid.Line
import edge_pre.dg.basis.Quad
import edge_pre.sc.grid.Quad
import edge_pre.dg.basis.Tria
import edge_pre.sc.grid.Tria
import edge_pre.dg.basis.Hex
import edge_pre.sc.grid.Hex
import edge_pre.types.Tet
import edge_pre.dg.basis.Tet


# set up logger
logging.basicConfig( level=logging.DEBUG,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# welcome our users
logging.info( '##########################################################################' )
logging.info( '##############   ##############            ###############  ##############' )
logging.info( '##############   ###############         ################   ##############' )
logging.info( '#####            #####       #####      ######                       #####' )
logging.info( '#####            #####        #####    #####                         #####' )
logging.info( '#############    #####         #####  #####                  #############' )
logging.info( '#############    #####         #####  #####      #########   #############' )
logging.info( '#####            #####         #####  #####      #########           #####' )
logging.info( '#####            #####        #####    #####        ######           #####' )
logging.info( '#####            #####       #####      #####       #####            #####' )
logging.info( '###############  ###############         ###############   ###############' )
logging.info( '###############  ##############           #############    ###############' )
logging.info( '#######################################################################pre' )
logging.info( 'you reached EDGEpre, the pre-processing engine of EDGE' )

# parse command line options
l_parser = argparse.ArgumentParser( description='Generation of pre-processed data structures for EDGE.' )

l_parser.add_argument( '-x', '--xml',
                       dest     = 'xml',
                       required = True,
                       type     = str,
                       help     = 'XML configuration of EDGEpre')
l_args = vars(l_parser.parse_args())

# parse XML-config
l_conf = edge_pre.io.Config.Config( l_args['xml'] )

# iterate over element types
for l_ty in l_conf.m_types:
  logging.info( 'processing element type: '+ l_ty )

  # iterate over polynomial degrees
  for l_de in l_conf.m_degs:
    logging.info( '  polynomial degree: '+ str(l_de) )

    # create dirs
    l_dirExt = '/' + l_ty + '/' + str(l_de) + '/'
    for l_out in l_conf.m_out.keys():
      l_dir = l_conf.m_out[l_out] + l_dirExt
      if not os.path.exists(l_dir):
        os.makedirs(l_dir)

    #
    # clear sympy cache to avoid memory issues
    #
    clear_cache()

    #
    # line elements
    #
    if( l_ty == 'line' ):
      # element type
      l_elTy = edge_pre.types.Line.Line( l_de )

      # set pseudo data
      l_symsFa = None
      l_basisFa = None
      l_intFa = None

      # get basis
      l_symsEl, l_basisEl = edge_pre.dg.basis.Line.gen( l_de )

      # get element integration interval
      l_intEl = l_elTy.intEl( l_symsEl[0] )

      # get face trafos
      l_faToFa = l_elTy.faToFa( l_symsFa )
      l_faToEl = l_elTy.faToEl( l_symsFa )

      # set dummy sub-cell reordering
      l_scDgAd = [ [0] ]

      # get sub-cell integration intervals
      l_intScIn, l_intScSend, l_intScSurf = edge_pre.sc.grid.Line.intSc( l_de, l_symsEl )

      # get sub-face integration matrices
      l_sfInt = []
      for l_fa in [0,1]:
        l_sfInt = l_sfInt + [ edge_pre.int.Matrices.subs( [(l_symsEl[0], l_elTy.ves[l_fa][0])],
                                                          l_basisEl ) ]

      # sub-cell adjacency
      l_scSfScIn, l_scSfScSend, l_scSfScRecv = edge_pre.sc.grid.Line.scSfSc( l_de )

      # type of sub-cell's faces
      l_scTySfIn, l_scTySfSend =  edge_pre.sc.grid.Line.scTySf( l_de )

      # plot basis
      l_path = l_conf.m_out['plots'] + l_dirExt + l_ty + '_' + str(l_de) + '_basis.pdf'
      edge_pre.dg.basis.Line.plot( l_path,
                                   l_symsEl,
                                   l_basisEl )

      # vertices of sub-cells
      l_veCrds = edge_pre.sc.grid.Line.svs( l_de )
      l_scSvIn, l_scSvSend, l_scSvRecv = edge_pre.sc.grid.Line.scSv( l_de )
      l_path = l_conf.m_out['plots'] + l_dirExt + l_ty + '_' + str(l_de) + '_subgrid.pdf'
      edge_pre.sc.grid.Line.plot( l_path,
                                  l_veCrds,
                                  l_scSvIn,
                                  l_scSvSend)

      # get scatter matrix for DG surface sub-cells
      l_scatterSurf = []
      # local face
      for l_fa in l_intScSurf:
        l_scatterSurf = l_scatterSurf + [ edge_pre.sc.ops.Project.scatter( l_basisEl, l_fa ) ]
      # remote face
      for l_fa in l_intScSurf:
        l_scatterSurf = l_scatterSurf + [ edge_pre.sc.ops.Project.scatter( l_basisEl, list( reversed(l_fa) ) ) ]

    #
    # quadrilaterals
    #
    if( l_ty == 'quad4r' ):
      # element type
      l_elTy = edge_pre.types.Quad.Quad( l_de )

      # get face basis
      l_symsFa, l_basisFa = edge_pre.dg.basis.Line.gen( l_de )

      # get element basis
      l_symsEl, l_basisEl = edge_pre.dg.basis.Quad.gen( l_de )

      # unify bases
      l_symsEl, l_basisEl = edge_pre.dg.basis.Mod.unify( l_symsFa, l_basisFa, l_symsEl, l_basisEl )

      # get face integration interval
      l_intFa = edge_pre.types.Line.Line( l_de ).intEl( l_symsFa )

      # get element integration interval
      l_intEl = l_elTy.intEl( l_symsEl )

      # get face trafos
      l_faToFa = l_elTy.faToFa( l_symsFa )
      l_faToEl = l_elTy.faToEl( l_symsFa )

      # get sub-cell reorderings
      l_scDgAd = edge_pre.sc.grid.Quad.scDgAd( l_de )

      # get sub-cell integration intervals
      l_intScIn, l_intScSend, l_intScSurf = edge_pre.sc.grid.Quad.intSc( l_de, l_symsEl )

      # get sub-face integration matrices
      l_subsSfDg, l_intSfDg = edge_pre.sc.grid.Quad.intSfDg( l_de, [sympy.symbols('chi_0')], l_symsEl )

      l_sfInt = []
      for l_fa in range(4):
        l_basisElSubs = edge_pre.int.Matrices.subsAll( l_subsSfDg[l_fa], l_basisEl )
        l_sfInt = l_sfInt + [ edge_pre.int.Matrices.intL( l_intSfDg[l_fa], l_basisElSubs ) ]

      # vertices of sub-cells
      l_veCrds = edge_pre.sc.grid.Quad.svs( l_de )
      l_scSvIn, l_scSvSend, l_scSvRecv = edge_pre.sc.grid.Quad.scSv( l_de )

      # sub-cell adjacency
      l_scSfScIn, l_scSfScSend, l_scSfScRecv = edge_pre.sc.grid.Quad.scSfSc( l_de )

      # type of sub-cell's faces
      l_scTySfIn, l_scTySfSend =  edge_pre.sc.grid.Quad.scTySf( l_de )

      # plot basis
      l_path = l_conf.m_out['plots'] + l_dirExt + l_ty + '_' + str(l_de) + '_basis.pdf'
      edge_pre.dg.basis.Quad.plot( l_path,
                                   l_symsEl,
                                   l_basisEl )

      # plot sub-grid
      l_path = l_conf.m_out['plots'] + l_dirExt + l_ty + '_' + str(l_de) + '_subgrid.pdf'
      edge_pre.sc.grid.Quad.plot( l_path,
                                  l_veCrds,
                                  l_scSvIn,
                                  l_scSvSend )

      # get scatter matrix for DG surface sub-cells
      l_scatterSurf = []
      # local face
      for l_fa in l_intScSurf:
        l_scatterSurf = l_scatterSurf + [ edge_pre.sc.ops.Project.scatter( l_basisEl, l_fa ) ]
      # remote face
      for l_fa in l_intScSurf:
        l_scatterSurf = l_scatterSurf + [ edge_pre.sc.ops.Project.scatter( l_basisEl, list( reversed(l_fa) ) ) ]

    #
    # triangles
    #
    if( l_ty == 'tria3' ):
      # element type
      l_elTy = edge_pre.types.Tria.Tria( l_de )

      # get face basis
      l_symsFa, l_basisFa = edge_pre.dg.basis.Line.gen( l_de )

      # get basis
      l_symsEl, l_basisEl = edge_pre.dg.basis.Tria.gen( l_de )

      # unify bases
      l_symsEl, l_basisEl = edge_pre.dg.basis.Mod.unify( l_symsFa, l_basisFa, l_symsEl, l_basisEl )

      # get face integration interval
      l_intFa = edge_pre.types.Line.Line( l_de ).intEl( l_symsFa )

      # get element integration interval
      l_intEl = l_elTy.intEl( l_symsEl )

      # get face trafos
      l_faToFa = l_elTy.faToFa( l_symsFa )
      l_faToEl = l_elTy.faToEl( l_symsFa )

      # get sub-cell reorderings
      l_scDgAd = edge_pre.sc.grid.Tria.scDgAd( l_de )

      # get sub-cell integration intervals
      l_intScIn, l_intScSend, l_intScSurf = edge_pre.sc.grid.Tria.intSc( l_de, l_symsEl )

      # get sub-face integration matrices
      l_subsSfDg, l_intSfDg = edge_pre.sc.grid.Tria.intSfDg( l_de, [sympy.symbols('chi_0')], l_symsEl )

      l_sfInt = []
      for l_fa in range(3):
        l_basisElSubs = edge_pre.int.Matrices.subsAll( l_subsSfDg[l_fa], l_basisEl )
        l_sfInt = l_sfInt + [ edge_pre.int.Matrices.intL( l_intSfDg[l_fa], l_basisElSubs ) ]

      # vertices of sub-cells
      l_veCrds = edge_pre.sc.grid.Tria.svs( l_de )
      l_scSvIn, l_scSvSend, l_scSvRecv = edge_pre.sc.grid.Tria.scSv( l_de )

      # sub-cell adjacency
      l_scSfScIn, l_scSfScSend, l_scSfScRecv = edge_pre.sc.grid.Tria.scSfSc( l_de )

      # type of sub-cell's faces
      l_scTySfIn, l_scTySfSend =  edge_pre.sc.grid.Tria.scTySf( l_de )

      # plot basis
      l_path = l_conf.m_out['plots'] + l_dirExt + l_ty + '_' + str(l_de) + '_basis.pdf'
      edge_pre.dg.basis.Tria.plot( l_path,
                                   l_symsEl,
                                   l_basisEl )

      # plot sub-grid
      l_path = l_conf.m_out['plots'] + l_dirExt + l_ty + '_' + str(l_de) + '_subgrid.pdf'
      edge_pre.sc.grid.Tria.plot( l_path,
                                  l_veCrds,
                                  l_scSvIn,
                                  l_scSvSend)

      # get scatter matrix for DG surface sub-cells
      l_scatterSurf = []
      # local face
      for l_fa in l_intScSurf:
        l_scatterSurf = l_scatterSurf + [ edge_pre.sc.ops.Project.scatter( l_basisEl, l_fa ) ]
      # remote face
      for l_fa in l_intScSurf:
         # derive reordered list
        l_faRe = [ l_fa[l_sc] for l_sc in l_scDgAd[0] ]
        l_scatterSurf = l_scatterSurf + [ edge_pre.sc.ops.Project.scatter( l_basisEl, l_faRe ) ]

    #
    # hexes
    #
    if( l_ty == 'hex8r' ):
      # element type
      l_elTy = edge_pre.types.Hex.Hex( l_de )

      # get face basis
      l_symsFa, l_basisFa = edge_pre.dg.basis.Quad.gen( l_de )

      # get element basis
      l_symsEl, l_basisEl = edge_pre.dg.basis.Hex.gen( l_de )

      # unify bases
      l_symsEl, l_basisEl = edge_pre.dg.basis.Mod.unify( l_symsFa, l_basisFa, l_symsEl, l_basisEl )

      # get face integration intervals
      l_intFa = edge_pre.types.Quad.Quad( l_de ).intEl( l_symsFa )

      # get element integration interval
      l_intEl = l_elTy.intEl( l_symsEl )

      # get face trafos
      l_faToFa = l_elTy.faToFa( l_symsFa )
      l_faToEl = l_elTy.faToEl( l_symsFa )

      # get sub-cell reorderings
      l_scDgAd = edge_pre.sc.grid.Hex.scDgAd( l_de )

      # get sub-cell integration intervals
      l_intScIn, l_intScSend, l_intScSurf = edge_pre.sc.grid.Hex.intSc( l_de, l_symsEl )

      # get sub-face integration matrices
      l_subsSfDg, l_intSfDg = edge_pre.sc.grid.Hex.intSfDg( l_de, [sympy.symbols('chi_0'), sympy.symbols('chi_1')], l_symsEl )

      l_sfInt = []
      for l_fa in range(6):
        l_basisElSubs = edge_pre.int.Matrices.subsAll( l_subsSfDg[l_fa], l_basisEl )
        l_sfInt = l_sfInt + [ edge_pre.int.Matrices.intL( l_intSfDg[l_fa], l_basisElSubs ) ]

      # vertices of sub-cells
      l_veCrds = edge_pre.sc.grid.Hex.svs( l_de )
      l_scSvIn, l_scSvSend, l_scSvRecv = edge_pre.sc.grid.Hex.scSv( l_de )

      # sub-cell adjacency
      l_scSfScIn, l_scSfScSend, l_scSfScRecv = edge_pre.sc.grid.Hex.scSfSc( l_de )

      # type of sub-cell's faces
      l_scTySfIn, l_scTySfSend =  edge_pre.sc.grid.Hex.scTySf( l_de )

      # plot sub-grid
      l_path = l_conf.m_out['plots'] + l_dirExt + l_ty + '_' + str(l_de) + '_subgrid.pdf'
      edge_pre.sc.grid.Hex.plot( l_path,
                                 l_de )

      # get scatter matrix for DG surface sub-cells
      l_scatterSurf = []
      # local face
      for l_fa in l_intScSurf:
        l_scatterSurf = l_scatterSurf + [ edge_pre.sc.ops.Project.scatter( l_basisEl, l_fa ) ]
      # remote face
      for l_ve in range(4):
        for l_fa in l_intScSurf:
          # derive reordered list
          l_faRe = [ l_fa[l_sc] for l_sc in l_scDgAd[l_ve] ]
          l_scatterSurf = l_scatterSurf + [ edge_pre.sc.ops.Project.scatter( l_basisEl, l_faRe ) ]

    #
    # tets
    #
    if( l_ty == 'tet4' ):
      # element type
      l_elTy = edge_pre.types.Tet.Tet( l_de )

      # get face basis
      l_symsFa, l_basisFa = edge_pre.dg.basis.Tria.gen( l_de )

      # get element basis
      l_symsEl, l_basisEl = edge_pre.dg.basis.Tet.gen( l_de )

      # unify bases
      l_symsEl, l_basisEl = edge_pre.dg.basis.Mod.unify( l_symsFa, l_basisFa, l_symsEl, l_basisEl )

      # get face integration interval
      l_intFa = edge_pre.types.Tria.Tria( l_de ).intEl( l_symsFa )

      # get element integration interval
      l_intEl = l_elTy.intEl( l_symsEl )

      # get face trafos
      l_faToFa = l_elTy.faToFa( l_symsFa )
      l_faToEl = l_elTy.faToEl( l_symsFa )

    #
    # Generic DG structures
    #
    # get mass matrix
    l_mass = edge_pre.int.Matrices.mass( l_intEl, l_basisEl )

    # get stiffness matrices
    l_stiff = edge_pre.int.Matrices.stiff( l_symsEl,
                                           l_basisEl,
                                           l_intEl )

    # stiffness time matrices, transposed stiffness matrices multiplied by inverse mass matrix (after transposing)
    l_stiffT = []
    for l_di in range(len(l_stiff)):
      l_stiffT = l_stiffT + [ l_stiff[l_di].transpose()*l_mass.inv() ]

    # stifness volume matrices, stiffness matrices multiplied by inverse mass matrix
    l_stiffV = []
    for l_di in range(len(l_stiff)):
      l_stiffV = l_stiffV + [ l_stiff[l_di]*l_mass.inv() ]

    # get flux matrices
    l_flux = edge_pre.int.Matrices.flux( l_basisFa,
                                         l_basisEl,
                                         l_symsFa,
                                         l_symsEl,
                                         l_intFa,
                                         l_faToFa,
                                         l_faToEl )

    # derive single step projection and surface int for local contribution
    l_fluxL = [ l_mat*l_flux[1] for l_mat in l_flux[0] ] # premultiply first projection with diagonal face mass matrix

    # derive single step projection and surface int for neigh contribution
    l_fluxN = []
    # Remark: vertices are slowest dim, since our adjacency ordering prevents all but the first case for non-periodic meshes
    for l_ve in range( len(l_flux[2]) ):
      for l_fa in range( len(l_flux[0]) ):
        l_fluxN = l_fluxN + [ l_flux[0][l_fa] * l_flux[2][l_ve] ]

    # pre-multiply back basis-change by inverse mass matrix
    l_fluxT = [ l_mat.transpose()*l_mass.inv() for l_mat in l_flux[0] ]

    # save mass matrix
    l_path = l_conf.m_out['dg'] + l_dirExt + l_ty + '_' + str(l_de) + '_mass.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write( edge_pre.io.ArrStr.float2d( l_mass.tolist() ) )

    # save stiffness matrices
    l_path = l_conf.m_out['dg'] + l_dirExt + l_ty + '_' + str(l_de) + '_stiffV.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write( edge_pre.io.ArrStr.float3d( [ l_ma.tolist() for l_ma in l_stiffV ] ) )

    l_path = l_conf.m_out['dg'] + l_dirExt + l_ty + '_' + str(l_de) + '_stiffT.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write( edge_pre.io.ArrStr.float3d( [ l_ma.tolist() for l_ma in l_stiffT ] ) )

    # save flux matrices
    l_path = l_conf.m_out['dg'] + l_dirExt + l_ty + '_' + str(l_de) + '_fluxL.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write( edge_pre.io.ArrStr.float3d( [ l_ma.tolist() for l_ma in l_fluxL ] ) )

    l_path = l_conf.m_out['dg'] + l_dirExt + l_ty + '_' + str(l_de) + '_fluxN.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write( edge_pre.io.ArrStr.float3d( [ l_ma.tolist() for l_ma in l_fluxN ] ) )

    l_path = l_conf.m_out['dg'] + l_dirExt + l_ty + '_' + str(l_de) + '_fluxT.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write( edge_pre.io.ArrStr.float3d( [ l_ma.tolist() for l_ma in l_fluxT ] ) )

    # plot mass matrix
    l_path = l_conf.m_out['plots'] + l_dirExt + l_ty + '_' + str(l_de) + '_mass.pdf'
    edge_pre.io.Matrices.sparsity( l_mass,
                                   l_path )

    # plot stiffness matrices (volume kernel)
    for l_ma in range( len(l_stiff) ):
      l_path = l_conf.m_out['plots'] + l_dirExt + l_ty + '_' + str(l_de) + '_stiffV'+str(l_ma)+'.pdf'
      edge_pre.io.Matrices.sparsity( l_stiffV[l_ma],
                                     l_path )

    # plot stiffness matrices (time kernel)
    for l_ma in range( len(l_stiff) ):
      l_path = l_conf.m_out['plots'] + l_dirExt + l_ty + '_' + str(l_de) + '_stiffT'+str(l_ma)+'.pdf'
      edge_pre.io.Matrices.sparsity( l_stiffT[l_ma],
                                     l_path )

    # plot flux matrices
    for l_ma in range( len(l_fluxL) ):
      l_path = l_conf.m_out['plots'] + l_dirExt + l_ty + '_' + str(l_de) + '_fluxL'+str(l_ma)+'.pdf'
      edge_pre.io.Matrices.sparsity( l_fluxL[l_ma],
                                     l_path )

    for l_ma in range( len(l_fluxN) ):
      l_path = l_conf.m_out['plots'] + l_dirExt + l_ty + '_' + str(l_de) + '_fluxN'+str(l_ma)+'.pdf'
      edge_pre.io.Matrices.sparsity( l_fluxN[l_ma],
                                     l_path )

    for l_ma in range( len(l_fluxT) ):
      l_path = l_conf.m_out['plots'] + l_dirExt + l_ty + '_' + str(l_de) + '_fluxT'+str(l_ma)+'.pdf'
      edge_pre.io.Matrices.sparsity( l_fluxT[l_ma],
                                     l_path )

    # TODO: Add sc-support for tets
    if( l_ty == 'tet4' ):
      continue

    #
    # Generic sub-cell structures
    #
    # get scatter matrix
    l_scatter = edge_pre.sc.ops.Project.scatter( l_basisEl, l_intScIn+l_intScSend )

    # get gather matrix
    l_gather = edge_pre.sc.ops.Project.gather( l_basisEl, l_intScIn+l_intScSend )

    # scale face int with inverse mass matrix
    for l_fa in range(len(l_sfInt)):
     l_sfInt[l_fa] = l_sfInt[l_fa] * l_mass.inv()

    # save vertex cords
    l_path = l_conf.m_out['subcell'] + l_dirExt + l_ty + '_' + str(l_de) + '_svcrds.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write(  edge_pre.io.ArrStr.float2d( l_veCrds ) )

    # save sCsV
    l_path = l_conf.m_out['subcell'] + l_dirExt + l_ty + '_' + str(l_de) + '_scsv.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write(  edge_pre.io.ArrStr.int2d( l_scSvIn+l_scSvSend+l_scSvRecv ) )

    # save sub-cell scSfSc
    l_path = l_conf.m_out['subcell'] + l_dirExt + l_ty + '_' + str(l_de) + '_scsfsc.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write(  edge_pre.io.ArrStr.int2d( l_scSfScIn+l_scSfScSend+l_scSfScRecv ) )

    # save type of sub-cells' faces
    l_path = l_conf.m_out['subcell'] + l_dirExt + l_ty + '_' + str(l_de) + '_sctysf.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write(  edge_pre.io.ArrStr.int2d( l_scTySfIn+l_scTySfSend ) )

    # save sub-cell reordering
    l_path = l_conf.m_out['subcell'] + l_dirExt + l_ty + '_' + str(l_de) + '_scdgad.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write(  edge_pre.io.ArrStr.int2d( l_scDgAd ) )

    # save scatter matrix
    l_path = l_conf.m_out['subcell'] + l_dirExt + l_ty + '_' + str(l_de) + '_scatter.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write( edge_pre.io.ArrStr.float2d( l_scatter.tolist() ) )

    # save scatter surf matrix
    l_path = l_conf.m_out['subcell'] + l_dirExt + l_ty + '_' + str(l_de) + '_scattersurf.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write( edge_pre.io.ArrStr.float3d( [ l_ss.tolist() for l_ss in l_scatterSurf ] ) )

    # save gather matrix
    l_path = l_conf.m_out['subcell'] + l_dirExt + l_ty + '_' + str(l_de) + '_gather.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write( edge_pre.io.ArrStr.float2d( l_gather.tolist() ) )

    # save sub-face integration matrices
    l_path = l_conf.m_out['subcell'] + l_dirExt + l_ty + '_' + str(l_de) + '_sfint.csv'
    with open( l_path, 'w' ) as l_fi:
      l_fi.write( edge_pre.io.ArrStr.float3d( [ l_sf.tolist() for l_sf in l_sfInt ] ) )
