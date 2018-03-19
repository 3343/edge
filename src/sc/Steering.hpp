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
 * Functions for steering of the sub-cell limiter.
 **/
#ifndef EDGE_SC_STEERING_HPP
#define EDGE_SC_STEERING_HPP

namespace edge {
  namespace sc {
    class Steering;
  }
}

/**
 * Steering of the sub-cell limiter.
 **/
class edge::sc::Steering {
  public:
    /**
     * Gets the admissiblity ids based on the number of performed time steps since the last synchronization.
     *   First entry: Admissibility id of the previous solution (corresponds to input time step id).
     *   Second entry: Admissiblity id for the candidate solution.
     *   Third and fourth entry: Admissiblity id for the limited solution.
     * 
     * @param i_nTs number of time steps since last synchronization.
     * @param o_admIds will be set to admissibility ids.
     *
     * @paramt TL_T_TS integral type of the time step id.
     **/
    template< typename TL_T_TS >
    static void getAdmIds( TL_T_TS        i_nTs,
                           unsigned short o_admIds[4] ) {
      if( (i_nTs%2) == 0 ) {
        o_admIds[0] = 0;
        o_admIds[1] = 1;
        o_admIds[2] = 2;
        o_admIds[3] = 3;
      }
      else {
        o_admIds[0] = 2;
        o_admIds[1] = 3;
        o_admIds[2] = 0;
        o_admIds[3] = 1;
      }
    }

    /**
     * Gets the ids of the subcell extreme based on the number of time steps since the last synchronization.
     *   First entry: Id of the previous solution.
     *   Second entry: Id of the candidate/limited solution.
     *
     * @param i_nTs number of time steps since last synchronization.
     * @param o_extIds will be set to ids of extrema.
     *
     * @paramt TL_T_TS integral type of the time step id.
     **/
    template< typename TL_T_TS >
    static void getExtIds( TL_T_TS        i_nTs,
                           unsigned short o_exIds[2] ) {
      o_exIds[0] =  i_nTs   %2;
      o_exIds[1] = (i_nTs+1)%2;
    }

    /**
     * Gets the ids of the subcell dofs based on the number of time steps since the last synchronization.
     *   First entry: Id of the previous solution.
     *   Second entry: Id of the candidate/limited solution.
     *
     * @param i_nTs number of time steps since last synchronization.
     * @param o_dofsIds will be set to ids of extrema.
     *
     * @paramt TL_T_TS integral type of the time step id.
     **/
    template< typename TL_T_TS >
    static void getDofsIds( TL_T_TS       i_nTs,
                           unsigned short o_dofsIds[2] ) {
      o_dofsIds[0] =  i_nTs   %2;
      o_dofsIds[1] = (i_nTs+1)%2;
    }


    /**
     * Resets the position of the limited elements' admissibility information.
     * Afterwards the first pointer holds the admissibility of the previous solution,
     * the second pointer that of the outdated "candidate solution",
     * and the third the that of the outdated "pre-previous solution". 
     *
     * @param i_nTs number of time steps since the last synchronization.
     * @param io_adm admissibility information.
     *
     * @paramt TL_T_TS integral type of the time step id.
     * @paramt TL_N_CRS number of fused simulations.
     **/
    template< typename       TL_T_TS,
              unsigned short TL_N_CRS >
    static void resetAdm( TL_T_TS   i_nTs,
                          bool    (*io_adm[4])[TL_N_CRS] ) {
      // get admissibility ids
      unsigned short l_admIds[4];
      getAdmIds( i_nTs, l_admIds );

      // store current info
      bool (*l_adm[4])[TL_N_CRS];
      for( unsigned short l_am = 0; l_am < 4; l_am++ )
        l_adm[l_am] = io_adm[l_am];

      // reorder the array
      for( unsigned short l_am = 0; l_am < 4; l_am++ )
        io_adm[l_am] = l_adm[ l_admIds[l_am] ];
    }

    /**
     * Resets the position of the limited elements' dofs.
     * Afterwards the first pointer holds the dofs of the previous solution.
     * The second those of the outdated next time step.
     *
     * @param i_nTs number of time steps since the last synchronization.
     * @param io_dofs sub-cell dofs.
     *
     * @paramt TL_T_TS intgral type of the time step ids.
     * @paramt TL_T_REAL floating point arithmetic.
     * @paramt TL_N_FAS number of faces.
     * @paramt TL_N_SFS number sub-faces per face.
     * @paramt TL_N_QTS number of quantities.
     * @paramt TL_N_CRS number of fused simulations.
     **/
    template< typename       TL_T_TS,
              typename       TL_T_REAL,
              unsigned short TL_N_FAS,
              unsigned short TL_N_SFS,
              unsigned short TL_N_QTS,
              unsigned short TL_N_CRS >
    static void resetDofs( TL_T_TS         i_nTs,
                           TL_T_REAL (* (* io_dofs[2]) [TL_N_FAS])[TL_N_QTS][TL_N_SFS][TL_N_CRS] ) {
      // get dof ids
      unsigned short l_dofsIds[2];
      getDofsIds( i_nTs, l_dofsIds );

      // store current info
      TL_T_REAL (* (* l_dofs[2]) [TL_N_FAS])[TL_N_QTS][TL_N_SFS][TL_N_CRS];
      for( unsigned short l_do = 0; l_do < 2; l_do++ ) l_dofs[l_do] = io_dofs[l_do];

      // reorder dofs
      for( unsigned short l_do = 0; l_do < 2; l_do++ ) io_dofs[l_do] = l_dofs[ l_dofsIds[l_do] ];
    }

    /**
     * Resets the position of the limited elements' extrema.
     * Afterwards the first poniter holds the extrema of the previous solution.
     * The second those of the outdated "candidate solution".
     *
     * @param i_nTs number of time steps since the last synchronization.
     * @param io_ext sub-cell extrema.
     *
     * @paramt TL_T_TS intgral type of the time step ids.
     * @paramt TL_T_REAL floating point arithmetic.
     * @paramt TL_N_QTS number of quantities.
     * @paramt TL_N_CRS number of fused simulations.
     **/
    template< typename       TL_T_TS,
              typename       TL_T_REAL,
              unsigned short TL_N_QTS,
              unsigned short TL_N_CRS >
    static void resetExt( TL_T_TS     i_nTs,
                          TL_T_REAL (*io_ext[2])[2][TL_N_QTS][TL_N_CRS] ) {
      // get extrema ids
      unsigned short l_extIds[2];
      getExtIds( i_nTs, l_extIds );

      // store current info
      TL_T_REAL (*l_ext[2])[2][TL_N_QTS][TL_N_CRS];
      for( unsigned short l_ex = 0; l_ex < 2; l_ex++ ) l_ext[l_ex] = io_ext[l_ex];

      // reorder extrema
      for( unsigned short l_ex = 0; l_ex < 2; l_ex++ ) io_ext[l_ex] = l_ext[ l_extIds[l_ex] ];
    }

    /**
     * Resets the counter for the number of times the elements were limited.
     *
     * @param i_lpFirst first limited plus element.
     * @param i_lpSize number of limited plus elements for which the counter is reset.
     * @param i_lpLi limited elements associated with limited plus elements.
     * @param o_limSync counters, defined on limited element, which will be reset.
     **/
    template< typename       TL_T_LID,
              unsigned short TL_N_CRS >
    static void resetLimSync( TL_T_LID             i_lpFirst,
                              TL_T_LID             i_lpSize,
                              TL_T_LID     const  *i_lpLi,
                              unsigned int       (*o_limSync)[TL_N_CRS] ) {
      // iterate over limited plus elements
      for( TL_T_LID l_lp = i_lpFirst; l_lp < i_lpFirst+i_lpSize; l_lp++ ) {
        // determine id of limited element
        TL_T_LID l_li = i_lpLi[l_lp];

        // reset counter for limited element
        if( l_li != std::numeric_limits< TL_T_LID >::max() ) {
          for( unsigned short l_cr = 0; l_cr < TL_N_CRS; l_cr++ )
            o_limSync[l_li][l_cr] = 0;
        }
      }
    }
};

#endif
