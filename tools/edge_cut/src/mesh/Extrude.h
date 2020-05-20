/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (breuer AT mytum.de)
 *
 * @section LICENSE
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
 * Extrudes a given 2D rectangular surface mesh.
 **/
#include <cstddef>
#include <vector>
#include <tuple>
#include <limits>
#include <ostream>

namespace edge_cut {
  namespace mesh {
    class Extrude;
  }
}

class edge_cut::mesh::Extrude {
  private:
    //! point in space
    typedef std::tuple< double, double, double > point;

    //! epsilon used for zero comparisons
    double m_eps = -1;

    //! target depth
    double m_zTarget = 0;

    //! number of coarsening levels
    std::size_t m_nLevels = 0;

    //! number of steps in depth
    std::size_t m_nSteps = 0;

    //! minimum coordinates of the surface vertices
    double m_veMinSurf[3] = {0};

    //! maximum coordinates of the surface vertices
    double m_veMaxSurf[3] = {0};

    //! minimum delta in x- and y-direction
    double m_dMin[2] = { std::numeric_limits< double >::max(),
                         std::numeric_limits< double >::max() };

    //! refinement ratio
    unsigned short m_refRatio = 1;

    //! number of points in x-direction; 0: fine, 1: coarse
    std::size_t m_nx[2] = { std::numeric_limits< std::size_t >::max(),
                            std::numeric_limits< std::size_t >::max() };

    //! number of points in y-direction; 0: fine, 1: coarse
    std::size_t m_ny[2] = { std::numeric_limits< std::size_t >::max(),
                            std::numeric_limits< std::size_t >::max() };

    //! vertex coordinates for the surface mesh
    std::vector< point > m_veCrdsSurf;

    //! front vertex coordinates
    std::vector< point > m_veCrdsFront;

    //! back vertex coordinates
    std::vector< point > m_veCrdsBack;

    //! left vertex coordinates
    std::vector< point > m_veCrdsLeft;

    //! right vertex coordinates
    std::vector< point > m_veCrdsRight;

    /**
     * Adds a side to the given output coordinates.
     *
     * @param i_nSteps umber of steps in depth.
     * @param i_zTarget target depth.
     * @param i_veCrdsBnd vertex coordinates of the one-dimensional boundary.
     * @param o_veCrds vector of coordinates to which the vertices are added.
     **/
    static void addSide( std::size_t                  i_nSteps,
                         double                       i_zTarget,
                         std::vector< point > const & i_veCrdsBnd,
                         std::vector< point >       & o_veCrds );

    /**
     * Writes the given triangles with a bottom left to top right diagonal.
     *
     * @param i_off offset of the vertices.
     * @param i_strideX stride in x-direction.
     * @param i_strideY stride in y-direction.
     * @param o_tria stream to which the triangles are written.
     *
     * @return number of written triangles.
     **/
    static unsigned short writeTriasBlTr( std::size_t    i_off,
                                          std::size_t    i_strideX,
                                          std::size_t    i_strideY,
                                          std::ostream & o_tria );

    /**
     * Writes the given triangles with a bottom right to top left diagonal.
     *
     * @param i_off offset of the vertices.
     * @param i_strideX stride in x-direction.
     * @param i_strideY stride in y-direction.
     * @param o_tria stream to which the triangles are written.
     *
     * @return number of written triangles.
     **/
    static unsigned short writeTriasBrTl( std::size_t    i_off,
                                          std::size_t    i_strideX,
                                          std::size_t    i_strideY,
                                          std::ostream & o_tria );

    /**
     * Writes the given triangles which splits the four quads at the center right.
     *
     * @param i_off offset of the vertices.
     * @param i_strideX stride in x-direction.
     * @param i_strideY stride in y-direction.
     * @param o_tria stream to which the triangles are written.
     *
     * @return number of written triangles.
     **/
    static unsigned short writeTriasCr( std::size_t    i_off,
                                        std::size_t    i_strideX,
                                        std::size_t    i_strideY,
                                        std::ostream & o_tria );

    /**
     * Writes the given triangles which splits the four quads at the center left.
     *
     * @param i_off offset of the vertices.
     * @param i_strideX stride in x-direction.
     * @param i_strideY stride in y-direction.
     * @param o_tria stream to which the triangles are written.
     *
     * @return number of written triangles.
     **/
    static unsigned short writeTriasCl( std::size_t    i_off,
                                        std::size_t    i_strideX,
                                        std::size_t    i_strideY,
                                        std::ostream & o_tria );

    /**
     * Writes the given triangles which splits the four quads at the center top.
     *
     * @param i_off offset of the vertices.
     * @param i_strideX stride in x-direction.
     * @param i_strideY stride in y-direction.
     * @param o_tria stream to which the triangles are written.
     *
     * @return number of written triangles.
     **/
    static unsigned short writeTriasCt( std::size_t    i_off,
                                        std::size_t    i_strideX,
                                        std::size_t    i_strideY,
                                        std::ostream & o_tria );

    /**
     * Writes the given triangles which splits the four quads at the center bottom.
     *
     * @param i_off offset of the vertices.
     * @param i_strideX stride in x-direction.
     * @param i_strideY stride in y-direction.
     * @param o_tria stream to which the triangles are written.
     *
     * @return number of written triangles.
     **/
    static unsigned short writeTriasCb( std::size_t    i_off,
                                        std::size_t    i_strideX,
                                        std::size_t    i_strideY,
                                        std::ostream & o_tria );

    /**
     * Writes the triangulation for the given number of points in OFF-format.
     * Assumes y as fast dimension.
     *
     * @param i_nx number of points in x-direction.
     * @param i_ny number of points in y-direction.
     * @param o_tria stream to which the triangulation is written.
     * @param i_nLevels number of coarsening levels.
     *
     * @return number of written triangles.
     **/
    static std::size_t writeTriaOff( std::size_t     i_nx,
                                     std::size_t     i_ny,
                                     std::ostream & o_tria,
                                     unsigned short  i_nLevels = 0 );

    /**
     * Writes the given points in OFF-format.
     *
     * @param i_pts points.
     * @param o_pts stream to which the points are written.
     **/
    static void writePtsOff( std::vector< point > const & i_pts,
                             std::ostream               & o_pts );

  public:
    /**
     * Constructor.
     *
     * @param i_nVes number of vertices.
     * @param i_veCrds vertex coordinates.
     * @param i_zTarget target z of the extrusion; has to be larger than differences in heigh-variation within the surface mesh.
     * @param i_nLevels number of levels used for mesh coarsening.
     * @param i_epsilon epsilon used for zero-comparisons.
     **/
    Extrude( std::size_t          i_nVes,
             double      const (* i_veCrds)[3],
             double               i_zTarget,
             unsigned short       i_nLevels=0,
             double               l_epsilon=1E-3 );

    /**
     * Writes the specified stream in OFF-format to the given streams.
     *
     * @param i_left left side panel.
     * @param i_right right side panel.
     * @param i_front front side panel.
     * @param i_back back side panel.
     * @param i_bottom bottom side panel.
     * @param i_top top which is the meshed surface.
     **/
    void writeOff( std::ostream & io_left,
                   std::ostream & io_right,
                   std::ostream & io_front,
                   std::ostream & io_back,
                   std::ostream & io_bottom,
                   std::ostream & io_top ) const;
};