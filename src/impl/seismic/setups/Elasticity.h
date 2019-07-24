/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2019, Alexander Breuer
 * Copyright (c) 2019, Regents of the University of California
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
 * Elasticity.
 **/
#ifndef EDGE_SEISMIC_SETUPS_ELASTICITY_H
#define EDGE_SEISMIC_SETUPS_ELASTICITY_H

namespace edge {
  namespace seismic {
    namespace setups {
      class Elasticity;
    }
  }
}

class edge::seismic::setups::Elasticity {
  public:
    /**
     * Computes the elastic part of the two-dimensional star matrices.
     *
     * @param i_rho density.
     * @param i_lam lame parameter lambda.
     * @param i_mu shear modulus.
     * @param i_jacInv inverse jacobian matrix of the element.
     * @param o_starE will be set to assembled star matrices.
     **/
    static void star( double       i_rho,
                      double       i_lam,
                      double       i_mu, 
                      double const i_jacInv[2][2],
                      double       o_starE[2][5 * 5] );

    /**
     * Computes the elastic part of the two-dimensional star matrices.
     * Only non-zeros are set (compressed sparse rows)
     *
     * @param i_rho density.
     * @param i_lam lame parameter lambda.
     * @param i_mu shear modulus.
     * @param i_jacInv inverse jacobian matrix of the element.
     * @param o_starE will be set to non-zeros of assembled star matrices.
     **/
    static void star( double       i_rho,
                      double       i_lam,
                      double       i_mu, 
                      double const i_jacInv[2][2],
                      double       o_starE[2][10] );

    /**
     * Computes the elastic part of the three-dimensional star matrices.
     *
     * @param i_rho density.
     * @param i_lam lame parameter lambda.
     * @param i_mu shear modulus.
     * @param i_jacInv inverse jacobian matrix of the element.
     * @param o_starE will be set to assembled star matrix.
     **/
    static void star( double       i_rho,
                      double       i_lam,
                      double       i_mu, 
                      double const i_jacInv[3][3],
                      double       o_starE[3][9 * 9] );

    /**
     * Computes the elastic part of the three-dimensional star matrices.
     * Only non-zeros are set (compressed sparse rows)
     *
     * @param i_rho density.
     * @param i_lam lame parameter lambda.
     * @param i_mu shear modulus.
     * @param i_jacInv inverse jacobian matrix of the element.
     * @param o_starE will be set to non-zeros of assembled star matrices.
     **/
    static void star( double       i_rho,
                      double       i_lam,
                      double       i_mu, 
                      double const i_jacInv[3][3],
                      double       o_starE[3][24] );
};

#endif