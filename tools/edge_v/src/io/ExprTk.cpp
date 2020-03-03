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
 * ExprTk interface.
 **/
#include "ExprTk.h"
#include "logging.h"

void edge_v::io::ExprTk::addVar( std::string const & i_varName,
                                 double            & i_var ) {
  m_symTab.add_variable( i_varName,
                         i_var );
}

void edge_v::io::ExprTk::compile( std::string const & i_expr ) {
  bool l_res = false;

  // add constants by default
  l_res = m_symTab.add_constants();
  EDGE_V_CHECK( l_res );

  // register the symbol table
  m_expr.register_symbol_table( m_symTab );

  // compile
  if( !m_parser.compile( i_expr, m_expr ) ) {
    EDGE_V_LOG_ERROR << "couldn't compile expression: " << i_expr;
    EDGE_V_LOG_ERROR << "here's what went wrong:";

    for( std::size_t l_er = 0; l_er < m_parser.error_count(); l_er++ ) {
      exprtk::parser_error::type l_pErr = m_parser.get_error(l_er);
      EDGE_V_LOG_ERROR << "  error #" << l_er+1 << ":";
      EDGE_V_LOG_ERROR << "    position: " << l_pErr.token.position;
      EDGE_V_LOG_ERROR << "    type:     " << exprtk::parser_error::to_str(l_pErr.mode);
      EDGE_V_LOG_ERROR << "    message:  " << l_pErr.diagnostic;
      EDGE_V_LOG_FATAL << "aborting";
    }
  }
}

void edge_v::io::ExprTk::eval() {
  m_expr.value();
}