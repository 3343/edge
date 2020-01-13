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
#ifndef EDGE_V_IO_EXPRTK_HPP
#define EDGE_V_IO_EXPRTK_HPP

#define exprtk_disable_enhanced_features
#define exprtk_disable_return_statement
#define exprtk_disable_string_capabilities
#define exprtk_disable_rtl_io_file
#include "submodules/exprtk/exprtk.hpp"

namespace edge_v {
  namespace io {
    class ExprTk;
  }
}

/**
 * ExprTk interface.
 **/
class edge_v::io::ExprTk {
  private:
    //! symbols
    exprtk::symbol_table< double > m_symTab;

    //! expressions
    exprtk::expression< double > m_expr;

    //! parser
    exprtk::parser< double > m_parser;

  public:
    /**
     * Adds a variable.
     *
     * @param i_varName name of the variable.
     * @param i_var variable which is added.
     **/
    void addVar( std::string const & i_varName,
                 double            & i_var );
    /**
     * Compiles the expression.
     *
     * @param i_exp exression which is compiled.
     **/
    void compile( std::string const & i_expr );

    /**
     * Evaluates the expression.
     **/
    void eval();
};

#endif