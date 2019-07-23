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
 * Friction laws.
 **/

#ifndef EDGE_SEISMIC_FRICTION_LAWS_HPP
#define EDGE_SEISMIC_FRICTION_LAWS_HPP

#include "FrictionLaws.type"
#include <cmath>

namespace edge {
  namespace seismic {
    namespace solvers {
      template< unsigned short TL_N_DIM, unsigned short TL_N_CRUNS >
      class FrictionLaws;

      template< unsigned short TL_N_CRUNS >
      class FrictionLaws< 2, TL_N_CRUNS>;

      template< unsigned short TL_N_CRUNS >
      class FrictionLaws< 3, TL_N_CRUNS > ;
    }
  }
}

/**
 * Support for 2D friction laws.
 **/
template< unsigned short TL_N_CRUNS >
class edge::seismic::solvers::FrictionLaws< 2, TL_N_CRUNS > {
  private:
    /**
     * Applies the linear slip weakening friction law in two dimensions.
     *
     *         _
     *        | mu_s - ( ( mu_s - mu_d ) / D_c ) * delta    if delta  < D_c
     * mu_f = |
     *        |_mu_d                                        if delta >= D_c
     *
     * Sketch:
     *
     *      *                /|\
     *      *            mu_f |
     * mu_s x                 |---->
     *      **                 delta
     *      * *
     *      *  *
     *      *   *
     * mu_d -----*********************
     *      *    |
     *      *    |                    
     *      *****|********************
     *          D_c
     *
     * Reference: Puente, Ampuero, Kaeser
     *            Dynamic rupture modeling on unstructured meshes using a discontinuous Galerkin method
     *            Journal of Geophysical Research, Vol. 114, 2009
     *
     * @param i_dt time step.
     * @param i_csDmuM shear wave speed divided by Lame parameter mu (minus side).
     * @param i_csDmuP shear wave speed divided by Lame parameter mu (plus side).
     * @param i_mus static friction coefficients.
     * @param i_mud dynamic friction coefficients.
     * @param i_dcInv inverse critical slip distance (1/Dc).
     * @param i_sn0 initial normal stress.
     * @param i_co0 cohesion.
     * @param i_ss0 initial shear stress.
     * @param i_ms middle state.
     * @param io_dd slip, will be updated with slip contribution of this step.
     * @param io_muf friction coefficient, will be updated according to the friction law.
     * @param o_st will be set to stregth of the fault.
     * @param o_sr will be set to slip rate at this point.
     * @param o_tr will be set to traction at this point.
     * @param o_msM will be set to perturbed minus side middle state.
     * @param o_msP will be set to perturbed plus side middle state.
     * @param o_per will be set to true if middle state was perturbed, false otherwise.
     *
     * @paramt TL_T_REAL precision in all computations.
     **/
    template< typename TL_T_REAL >
    static void inline linSlipWeak( TL_T_REAL       i_dt,
                                    TL_T_REAL const i_csDmuM,
                                    TL_T_REAL const i_csDmuP,
                                    TL_T_REAL const i_mus[TL_N_CRUNS],
                                    TL_T_REAL const i_mud[TL_N_CRUNS],
                                    TL_T_REAL const i_dcInv[TL_N_CRUNS],
                                    TL_T_REAL const i_sn0[TL_N_CRUNS],
                                    TL_T_REAL const i_co0[TL_N_CRUNS],
                                    TL_T_REAL const i_ss0[TL_N_CRUNS],
                                    TL_T_REAL const i_ms[5][TL_N_CRUNS],
                                    TL_T_REAL       io_dd[TL_N_CRUNS],
                                    TL_T_REAL       io_muf[TL_N_CRUNS],
                                    TL_T_REAL       o_st[TL_N_CRUNS],
                                    TL_T_REAL       o_sr[TL_N_CRUNS],
                                    TL_T_REAL       o_tr[TL_N_CRUNS],
                                    TL_T_REAL       o_msM[5][TL_N_CRUNS],
                                    TL_T_REAL       o_msP[5][TL_N_CRUNS],
                                    bool            o_per[TL_N_CRUNS] ) {
      // init minus and plus "perturbed" middle states
      for( unsigned short l_qt = 0; l_qt < 5; l_qt++ ) {
        for( unsigned short l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
          o_msM[l_qt][l_ru] = i_ms[l_qt][l_ru];
          o_msP[l_qt][l_ru] = i_ms[l_qt][l_ru];
        }
      }

      // determine if the fault fails
      for( unsigned short l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
        // fault strength
        o_st[l_ru]  = io_muf[l_ru] * ( i_sn0[l_ru] + i_ms[0][l_ru] );
        o_st[l_ru] += i_co0[l_ru];

        // strength is only relevant for negative normal stress (compression) + switch sign
        o_st[l_ru] = (o_st[l_ru] < 0) ? -o_st[l_ru] : 0;

        // total shear stress
        TL_T_REAL l_shear = i_ss0[l_ru] + i_ms[2][l_ru];

        // eval failure criterion
        o_per[l_ru] = ( std::abs(l_shear) > o_st[l_ru] ) ? true : false;
      }

      // perturb middle states and update friction coefficient
      for( unsigned short l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
        // compute traction
        o_tr[l_ru] = o_st[l_ru] - std::abs(i_ss0[l_ru]);
        o_tr[l_ru] = (i_ss0[l_ru] > 0) ? o_tr[l_ru] : -o_tr[l_ru];

        // fall back to middle state if the fault is locked
        o_tr[l_ru] = (o_per[l_ru]) ? o_tr[l_ru] : i_ms[2][l_ru];

        // apply difference in 2nd and 4th wave strength to shear stress
        o_msM[2][l_ru] = o_tr[l_ru];
        o_msP[2][l_ru] = o_tr[l_ru];

        // apply difference in second and fourth wave strength to the remaining non-zero
        // component (fault-tangent particle velocity) of the second and fourth eigenvector
        o_msM[4][l_ru] -= i_csDmuM * (o_tr[l_ru] - i_ms[2][l_ru]);
        o_msP[4][l_ru] += i_csDmuP * (o_tr[l_ru] - i_ms[2][l_ru]);

        // compute slip rate
        o_sr[l_ru] = o_msP[4][l_ru] - o_msM[4][l_ru];

        // compute resulting slip
        io_dd[l_ru] += std::abs( o_sr[l_ru] ) * i_dt;

        // update friction coefficient
        TL_T_REAL  l_arg  = -(i_mus[l_ru] - i_mud[l_ru]) * i_dcInv[l_ru];
                   l_arg *= io_dd[l_ru];
                   l_arg += i_mus[l_ru];
        io_muf[l_ru] = std::max( i_mud[l_ru], l_arg );

        // add background shear stress to traction for output
        o_tr[l_ru] += i_ss0[l_ru];
      }
    }

  public:
    /**
     * Applies the linear slip weakening friction law in two dimensions.
     *
     * @param i_dt time step.
     * @param i_lswGlobal per-simulation global parameters of linear slip weakening, shared among all faces.
     * @param i_lswFace per-face parameters of linear slip weakening, private from face to face but shared on among the fused simulations.
     * @param i_ms middle state.
     * @param io_lswSf private data per-simulation and sub-face of each DG-face.
     * @param o_msL perturbed left-side middle state.
     * @param o_msR perturbed right-side middle state.
     * @param o_per will be set to true if the middle state was perturbed, false otherwise.
     *
     * @paramt TL_T_REAL precision in all computations.
     **/
    template< typename TL_T_REAL >
    static void perturb( TL_T_REAL                                                i_dt,
                         t_LinSlipWeakGlobal<TL_T_REAL, TL_N_CRUNS>        const &i_lswGlobal,
                         t_LinSlipWeakFace<TL_T_REAL>                      const &i_lswFace,
                         TL_T_REAL                                         const  i_ms[5][TL_N_CRUNS],
                         t_LinSlipWeakSubFace<
                           TL_T_REAL,
                           2,
                           TL_N_CRUNS >                                          &io_lswSf,
                         TL_T_REAL                                                o_msL[5][TL_N_CRUNS],
                         TL_T_REAL                                                o_msR[5][TL_N_CRUNS],
                         bool                                                     o_per[TL_N_CRUNS] ) {
      linSlipWeak( i_dt,
                   i_lswFace.csDmuM,
                   i_lswFace.csDmuP,
                   i_lswGlobal.mus,
                   i_lswGlobal.mud,
                   i_lswGlobal.dcInv,
                   io_lswSf.sn0,
                   io_lswSf.co0,
                   io_lswSf.ss0[0],
                   i_ms,
                   io_lswSf.dd[0],
                   io_lswSf.muf,
                   io_lswSf.st,
                   io_lswSf.sr[0],
                   io_lswSf.tr[0],
                   (i_lswFace.lEqM) ? o_msL : o_msR,
                   (i_lswFace.lEqM) ? o_msR : o_msL,
                   o_per );
    }

    /**
     * Applies the given friction law (derived from face data type) by perturbing the middle states.
     *
     * @param i_sf id of the sub-face.
     * @param i_dt "time step" of this pertubation. used for slip computation, not to be confused of the time step of seismic wave propagation.
     * @param i_ms middle states for all fused simulations.
     * @param io_faData data of the friction law at the face.
     * @param o_msL will be set to middle states at the left-side for all fused simulations.
     * @param o_msR will be set to middle states at the right-side for all fused simulations.
     * @param o_per will be set to true if the middle state was perturbed, false otherwise.
     *
     * @paramt TL_T_REAL real type used for arithmetic oprations.
     * @paramt TL_T_FA_DATA face data of the friction law.
     **/
    template< typename TL_T_REAL, typename TL_T_FA_DATA >
    static void perturb( unsigned short      i_sf,
                         TL_T_REAL           i_dt,
                         TL_T_REAL    const  i_ms[5][TL_N_CRUNS],
                         TL_T_FA_DATA       *io_faData,
                         TL_T_REAL           o_msL[5][TL_N_CRUNS],
                         TL_T_REAL           o_msR[5][TL_N_CRUNS],
                         bool                o_per[TL_N_CRUNS] ) {
      perturb( i_dt,
               *(io_faData->gl),
               *(io_faData->fa),
               i_ms,
              (*(io_faData->sf))[i_sf],
               o_msL,
               o_msR,
               o_per );
    }
};

/**
 * Support for 3D friction laws.
 **/
template< unsigned short TL_N_CRUNS >
class edge::seismic::solvers::FrictionLaws< 3, TL_N_CRUNS > {
  private:
    /**
     * Applies the linear slip weakening friction law in three dimensions.
     *         _
     *        | mu_s - ( ( mu_s - mu_d ) / D_c ) * delta    if delta  < D_c
     * mu_f = |
     *        |_mu_d                                        if delta >= D_c
     *
     * Sketch:
     *
     *      *                /|\
     *      *            mu_f |
     * mu_s x                 |---->
     *      **                 delta
     *      * *
     *      *  *
     *      *   *
     * mu_d -----*********************
     *      *    |
     *      *    |                    
     *      *****|********************
     *          D_c
     *
     * Reference: Pelties, Gabriel, Ampuero
     *            Verification of an ADER-DG method for complex dynamic rupture problems.
     *            Geosci. Model Dev., 7, 2014
     *
     * @param i_dt time step.
     * @param i_csDmuM shear wave speed divided by Lame parameter mu (minus side).
     * @param i_csDmuP shear wave speed divided by Lame parameter mu (plus side).
     * @param i_mus static friction coefficients.
     * @param i_mud dynamic friction coefficients.
     * @param i_dcInv inverse critical slip distance (1/Dc).
     * @param i_sn0 initial normal stress.
     * @param i_co0 cohesion.
     * @param i_ss0 initial shear stresses. [0]: along-strike, [1]: along-dip.
     * @param i_ss0A absolute initial shear stress (sqrt( i_ss0[0]*i_ss0[0] + i_ss0[1]*i_ss0[1]).
     * @param i_ms middle state.
     * @param io_dd slip, will be updated with slip contribution of this step. [0]: along-strike slip, [1]: along-dip slip.
     * @param io_muf friction coefficient, will be updated according to the friction law.
     * @param o_st will be set to strength of the fault.
     * @param o_sr will be set to slip rate at this point. [0]: along-strike slip rate, [1]: along-dip slip rate.
     * @param o_tr will be set to traction at this point. [0]: along-strike traction, [1]: along-dip traction.
     * @param o_msM will be set to perturbed minus side middle state.
     * @param o_msP will be set to perturbed plus side middle state.
     * @param o_per will be set to true if the middle state was perturbed, false otherwise.
     **/
    template< typename TL_T_REAL >
    static void inline linSlipWeak( TL_T_REAL       i_dt,
                                    TL_T_REAL const i_csDmuM,
                                    TL_T_REAL const i_csDmuP,
                                    TL_T_REAL const i_mus[TL_N_CRUNS],
                                    TL_T_REAL const i_mud[TL_N_CRUNS],
                                    TL_T_REAL const i_dcInv[TL_N_CRUNS],
                                    TL_T_REAL const i_sn0[TL_N_CRUNS],
                                    TL_T_REAL const i_co0[TL_N_CRUNS],
                                    TL_T_REAL const i_ss0[2][TL_N_CRUNS],
                                    TL_T_REAL const i_ss0A[TL_N_CRUNS],
                                    TL_T_REAL const i_ms[9][TL_N_CRUNS],
                                    TL_T_REAL       io_dd[2][TL_N_CRUNS],
                                    TL_T_REAL       io_muf[TL_N_CRUNS],
                                    TL_T_REAL       o_st[TL_N_CRUNS],
                                    TL_T_REAL       o_sr[2][TL_N_CRUNS],
                                    TL_T_REAL       o_tr[2][TL_N_CRUNS],
                                    TL_T_REAL       o_msM[9][TL_N_CRUNS],
                                    TL_T_REAL       o_msP[9][TL_N_CRUNS],
                                    bool            o_per[TL_N_CRUNS] ) {
      // init minus and plus "perturbed" middle states
      for( unsigned short l_qt = 0; l_qt < 9; l_qt++ ) {
        for( unsigned short l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
          o_msM[l_qt][l_ru] = i_ms[l_qt][l_ru];
          o_msP[l_qt][l_ru] = i_ms[l_qt][l_ru];
        }
      }

      // combined total shear stress
      TL_T_REAL l_shear[TL_N_CRUNS];

      // determine if the fault fails
      for( unsigned short l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
        // fault strength
        o_st[l_ru]  = io_muf[l_ru] * ( i_sn0[l_ru] + i_ms[0][l_ru] );
        o_st[l_ru] += i_co0[l_ru];

        // strength is only relevant for negative normal stress (compression) + switch sign
        o_st[l_ru] = (o_st[l_ru] < 0) ? -o_st[l_ru] : 0;

        // total shear stress
        TL_T_REAL l_shear1 = i_ss0[0][l_ru] + i_ms[3][l_ru];
        TL_T_REAL l_shear2 = i_ss0[1][l_ru] + i_ms[5][l_ru];

        l_shear[l_ru] = std::sqrt( l_shear1*l_shear1 + l_shear2*l_shear2 );

        // eval failure criterion
        o_per[l_ru] = ( l_shear[l_ru] > o_st[l_ru] ) ? true : false;
      }

      // perturb middle states and update friction coefficient
      for( unsigned short l_ru = 0; l_ru < TL_N_CRUNS; l_ru++ ) {
        // compute traction
        TL_T_REAL l_scale = o_st[l_ru] / l_shear[l_ru];

        o_tr[0][l_ru] = (i_ss0[0][l_ru] + i_ms[3][l_ru]) * l_scale - i_ss0[0][l_ru];
        o_tr[1][l_ru] = (i_ss0[1][l_ru] + i_ms[5][l_ru]) * l_scale - i_ss0[1][l_ru];

        // fall back to middle state if the fault is locked
        o_tr[0][l_ru] = (o_per[l_ru]) ? o_tr[0][l_ru] : i_ms[3][l_ru];
        o_tr[1][l_ru] = (o_per[l_ru]) ? o_tr[1][l_ru] : i_ms[5][l_ru];

        // perturb shear stresses
        o_msM[3][l_ru] = o_tr[0][l_ru];
        o_msM[5][l_ru] = o_tr[1][l_ru];
        o_msP[3][l_ru] = o_tr[0][l_ru];
        o_msP[5][l_ru] = o_tr[1][l_ru];

        TL_T_REAL l_diff1 = o_tr[0][l_ru] - i_ms[3][l_ru];
        TL_T_REAL l_diff2 = o_tr[1][l_ru] - i_ms[5][l_ru];

        // apply difference to fault parallel velocities
        o_msM[7][l_ru] -= i_csDmuM * l_diff1;
        o_msM[8][l_ru] -= i_csDmuM * l_diff2;

        o_msP[7][l_ru] += i_csDmuP * l_diff1;
        o_msP[8][l_ru] += i_csDmuP * l_diff2;


        // compute slip rate
        o_sr[0][l_ru] = o_msP[7][l_ru] - o_msM[7][l_ru];
        o_sr[1][l_ru] = o_msP[8][l_ru] - o_msM[8][l_ru];

        // compute resulting slip
        io_dd[0][l_ru] += std::abs(o_sr[0][l_ru]) * i_dt;
        io_dd[1][l_ru] += std::abs(o_sr[1][l_ru]) * i_dt;

        TL_T_REAL l_ddA = std::sqrt( io_dd[0][l_ru] * io_dd[0][l_ru] + io_dd[1][l_ru] * io_dd[1][l_ru] );

        // update friction coefficient
        TL_T_REAL  l_arg  = -(i_mus[l_ru] - i_mud[l_ru]) * i_dcInv[l_ru];
                   l_arg *= l_ddA;
                   l_arg += i_mus[l_ru];
        io_muf[l_ru] = std::max( i_mud[l_ru], l_arg );

        // add background shear stress to traction for output
        o_tr[0][l_ru] += i_ss0[0][l_ru];
        o_tr[1][l_ru] += i_ss0[1][l_ru];
      }
    }

  public:
    /**
     * Applies the linear slip weakening friction law in three dimensions.
     *
     * @param i_dt distance between this quad point and the previous one (or wave prop time step), see linSlipWeak for details.
     * @param i_lswGlobal per-simulation global parameters of linear slip weakening, shared among all faces.
     * @param i_lswFace per-face parameters of linear slip weakening, private from face to face but shared on among the fused simulations.
     * @param i_ms middle state.
     * @param io_lswSf private data per-simulation and sub-face of each DG-face.
     * @param o_msL perturbed left-side middle state.
     * @param o_msR perturbed right-side middle state.
     * @param o_per will be set to true if the middle state was perturbed, false otherwise.
     *
     * @paramt TL_T_REAL precision in all computations.
     **/
    template< typename TL_T_REAL >
    static void perturb( TL_T_REAL                                            i_dt,
                         t_LinSlipWeakGlobal< TL_T_REAL, TL_N_CRUNS > const & i_lswGlobal,
                         t_LinSlipWeakFace< TL_T_REAL >               const & i_lswFace,
                         TL_T_REAL                                    const   i_ms[9][TL_N_CRUNS],
                         t_LinSlipWeakSubFace<
                           TL_T_REAL,
                           3,
                           TL_N_CRUNS >                                     & io_lswSf,
                         TL_T_REAL                                            o_msL[9][TL_N_CRUNS],
                         TL_T_REAL                                            o_msR[9][TL_N_CRUNS],
                         bool                                                 o_per[TL_N_CRUNS] ) {
      linSlipWeak( i_dt,
                   i_lswFace.csDmuM,
                   i_lswFace.csDmuP,
                   i_lswGlobal.mus,
                   i_lswGlobal.mud,
                   i_lswGlobal.dcInv,
                   io_lswSf.sn0,
                   io_lswSf.co0,
                   io_lswSf.ss0,
                   io_lswSf.ss0A,
                   i_ms,
                   io_lswSf.dd,
                   io_lswSf.muf,
                   io_lswSf.st,
                   io_lswSf.sr,
                   io_lswSf.tr,
                   (i_lswFace.lEqM) ? o_msL : o_msR,
                   (i_lswFace.lEqM) ? o_msR : o_msL,
                   o_per );
    }

    /**
     * Applies the given friction law (derived from face data type) by perturbing the middle states.
     *
     * @param i_sf id of the sub-face.
     * @param i_dt "time step" of this pertubation. used for slip computation, not to be confused of the time step of seismic wave propagation.
     * @param i_ms middle states for all fused simulations.
     * @param io_faData data of the friction law at the face.
     * @param o_msL will be set to middle states at the left-side for all fused simulations.
     * @param o_msR will be set to middle states at the right-side for all fused simulations.
     * @param o_per will be set to true if the middle state was perturbed, false otherwise.
     *
     * @paramt TL_T_REAL real type used for arithmetic oprations.
     * @paramt TL_T_FA_DATA face data of the friction law.
     **/
    template< typename TL_T_REAL, typename TL_T_FA_DATA >
    static void perturb( unsigned short         i_sf,
                         TL_T_REAL              i_dt,
                         TL_T_REAL      const   i_ms[9][TL_N_CRUNS],
                         TL_T_FA_DATA         * io_faData,
                         TL_T_REAL              o_msL[9][TL_N_CRUNS],
                         TL_T_REAL              o_msR[9][TL_N_CRUNS],
                         bool                   o_per[TL_N_CRUNS] ) {
      perturb( i_dt,
               *(io_faData->gl),
               *(io_faData->fa),
               i_ms,
              (*(io_faData->sf))[i_sf],
               o_msL,
               o_msR,
               o_per );
    }
};

#endif
