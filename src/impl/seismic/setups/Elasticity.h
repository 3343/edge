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
  private:
    /**
     * Computes the left side's contribution to the two-dimensional Riemann solution in face-aligned coordinates.
     *
     * @param i_rhoL left element's density rho.
     * @param i_rhoR right element's density rho.
     * @param i_lamL left element's Lame parameter lambda.
     * @param i_lamR right element's Lame parameter lambda.
     * @param i_muL left element's Lame parameter mu.
     * @param i_muR right element's Lame parameter mu.
     * @param o_flMid will be set to matrix for the left element's contribution.
     **/
    static void setupFlMid( double i_rhoL,
                            double i_rhoR,
                            double i_lamL,
                            double i_lamR,
                            double i_muL,
                            double i_muR,
                            double o_flMid[5][5] );

    /**
     * Computes the right side's contribution to the two-dimensional Riemann solution in face-aligned coordinates.
     *
     * @param i_rhoL left element's density rho.
     * @param i_rhoR right element's density rho.
     * @param i_lamL left element's Lame parameter lambda.
     * @param i_lamR right element's Lame parameter lambda.
     * @param i_muL left element's Lame parameter mu.
     * @param i_muR right element's Lame parameter mu.
     * @param o_frMid will be set to matrix for the right element's contribution.
     **/
    static void setupFrMid( double i_rhoL,
                            double i_rhoR,
                            double i_lamL,
                            double i_lamR,
                            double i_muL,
                            double i_muR,
                            double o_frMid[5][5] );

    /**
     * Computes the left side's contribution to the three-dimensional Riemann solution in face-aligned coordinates.
     *
     * @param i_rhoL left element's density rho.
     * @param i_rhoR right element's density rho.
     * @param i_lamL left element's Lame parameter lambda.
     * @param i_lamR right element's Lame parameter lambda.
     * @param i_muL left element's Lame parameter mu.
     * @param i_muR right element's Lame parameter mu.
     * @param o_flMid will be set to matrix for the left element's contribution.
     **/
    static void setupFlMid( double i_rhoL,
                            double i_rhoR,
                            double i_lamL,
                            double i_lamR,
                            double i_muL,
                            double i_muR,
                            double o_flMid[9][9] );

    /**
     * Computes the right side's contribution to the three-dimensional Riemann solution in face-aligned coordinates.
     *
     * @param i_rhoL left element's density rho.
     * @param i_rhoR right element's density rho.
     * @param i_lamL left element's Lame parameter lambda.
     * @param i_lamR right element's Lame parameter lambda.
     * @param i_muL left element's Lame parameter mu.
     * @param i_muR right element's Lame parameter mu.
     * @param o_frMid will be set to matrix for the right element's contribution.
     **/
    static void setupFrMid( double i_rhoL,
                            double i_rhoR,
                            double i_lamL,
                            double i_lamR,
                            double i_muL,
                            double i_muR,
                            double o_frMid[9][9] );
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

    /**
     * Sets up the two-dimensional flux solvers for a single face.
     *
     * @param i_rhoL density of the left element.
     * @param i_rhoL density of the right element.
     * @param i_lamL Lame parameter lambda of the left element.
     * @param i_lamR Lame parameter lambda of the right element.
     * @param i_muL Lame parameter mu of the left element.
     * @param i_muR Lame parameter mu of the right element.
     * @param i_t transformation of stresses from face-aligned coordinates to Cartesian.
     * @param i_tm1 transformation of elastic quantities from Cartesian coordinates to face-aligned.
     * @param o_fsEl will be set to elastic flux solver, applied to the left element's DOFs.
     * @param o_fsEr will be set to elastic flux solver, applied to the right element's DOFs.
     * @param i_freeSurface true if free surface boundary conditions are applied to the right element.
     **/
    static void fs( double       i_rhoL,
                    double       i_rhoR,
                    double       i_lamL,
                    double       i_lamR,
                    double       i_muL,
                    double       i_muR,
                    double const i_t[5][5],
                    double const i_tm1[5][5],
                    double       o_fsEl[5*5],
                    double       o_fsEr[5*5],
                    bool         i_freeSurface );

    /**
     * Sets up the three-dimensional flux solvers for a single face.
     *
     * @param i_rhoL density of the left element.
     * @param i_rhoL density of the right element.
     * @param i_lamL Lame parameter lambda of the left element.
     * @param i_lamR Lame parameter lambda of the right element.
     * @param i_muL Lame parameter mu of the left element.
     * @param i_muR Lame parameter mu of the right element.
     * @param i_t transformation of stresses from face-aligned coordinates to Cartesian.
     * @param i_tm1 transformation of elastic quantities from Cartesian coordinates to face-aligned.
     * @param o_fsl will be set to elastic flux solver, applied to the left element's DOFs.
     * @param o_fsr will be set to elastic flux solver, applied to the right element's DOFs.
     * @param i_freeSurface true if free surface boundary conditions are applied to the right element.
     **/
    static void fs( double       i_rhoL,
                    double       i_rhoR,
                    double       i_lamL,
                    double       i_lamR,
                    double       i_muL,
                    double       i_muR,
                    double const i_t[9][9],
                    double const i_tm1[9][9],
                    double       o_fsEl[9*9],
                    double       o_fsEr[9*9],
                    bool         i_freeSurface );
};

#endif