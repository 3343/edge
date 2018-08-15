/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2017-2018, Regents of the University of California
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
 * Expressions in EDGE.
 **/
#ifndef EDGE_DATA_EXPRESSION_HPP
#define EDGE_DATA_EXPRESSION_HPP

#include <string>
#include "constants.hpp"
#include "exprtk.wrap.hpp"
#include "io/logging.h"

namespace edge {
  namespace data {
    template< typename TL_T_REAL >
    class Expression;
  }
}

template< typename TL_T_REAL >
class edge::data::Expression {
  private:
    //! symbol table
    exprtk::symbol_table< TL_T_REAL > m_symTab;

    //! expression
    exprtk::expression< TL_T_REAL > m_expr;

    //! parser
    exprtk::parser< TL_T_REAL > m_parser;

  public:
    /**
     * Binds the location of the scalar to the symbol.
     *
     * @param i_sym symbol.
     * @param i_sca scalar.
     **/
    void bind( std::string const & i_sym,
               TL_T_REAL         & i_sca ) {
      bool l_err = m_symTab.add_variable( i_sym, i_sca );
      EDGE_CHECK( l_err ) << "failed binding " << i_sym;
    }

    /**
     * Binds the location of the array to the symbol.
     *
     * @param i_sym symbol.
     * @param i_arr array.
     * @param i_size size of the array.
     **/
    void bind( std::string const & i_sym,
               TL_T_REAL         * i_arr,
               std::size_t         i_size ) {
      bool l_err = m_symTab.add_vector( i_sym, i_arr, i_size );
      EDGE_CHECK( l_err ) << "failed binding " << i_sym;
    }

    /**
     * Shortcut for the binding of a coordinate array.
     *
     * @param i_crds array of the coordinates.
     * @param i_nDims number of dimensions.
     **/
    void bindCrds( TL_T_REAL      * i_crds,
                   unsigned short   i_nDims ) {
      EDGE_CHECK_GE( i_nDims, 1 );

      std::string l_strs[3] = {"x", "y", "z"};

      for( unsigned short l_di = 0; l_di < i_nDims; l_di++ ) {
        bind( l_strs[l_di], i_crds[l_di] );
      }
    }

    /**
     * Compiles the expression.
     *
     * @param i_exprStr expression string which is compiled. if the length of the string is zero, a dummy expression is compiled instead.
     **/
    void compile( std::string const & i_exprStr ) {
      // use empty expression if nothing is given
      std::string l_exprStr;
      if( i_exprStr.length() > 0 ) l_exprStr = i_exprStr;
      else l_exprStr = "0==0;";

      // add constants by default
      bool l_err = m_symTab.add_constants();
      EDGE_CHECK( l_err ) << "failed adding constants";

      // register the symbol table
      m_expr.register_symbol_table( m_symTab );

      // compile
      if( !m_parser.compile( l_exprStr, m_expr ) ) {
        EDGE_LOG_ERROR << "couldn't compile expression: " << l_exprStr;
        EDGE_LOG_ERROR << "here's what went wrong:";

        for( std::size_t l_er = 0; l_er < m_parser.error_count(); l_er++ ) {
          exprtk::parser_error::type l_pErr = m_parser.get_error(l_er);
          EDGE_LOG_ERROR << "  error #" << l_er+1 << ":";
          EDGE_LOG_ERROR << "    position: " << l_pErr.token.position;
          EDGE_LOG_ERROR << "    type:     " << exprtk::parser_error::to_str(l_pErr.mode);
          EDGE_LOG_ERROR << "    message:  " << l_pErr.diagnostic;
          EDGE_LOG_FATAL << "aborting";
        }
      }
    }

    /**
     * Evaluates the expression.
     **/
    void eval() const {
      m_expr.value();
    }
};

// explicit instantiations
extern template class edge::data::Expression<float>;
extern template class edge::data::Expression<double>;

#endif
