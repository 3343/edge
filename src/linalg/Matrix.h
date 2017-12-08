/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016, Regents of the University of California
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
 * Matrices.
 **/

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iomanip>
#include <sstream>
#include <io/logging.h>
#include <cmath>
#include "constants.hpp"

/*
 * Sparse matrix format: Coordinate.
 *
 * Remark: This high level representation is only used as intermediate format and never as part of FEM computations.
 */
typedef struct {
  // non-zero values
  std::vector< real_base > nz;
  // rows of non-zero entries
  std::vector< unsigned int > ro;
  // cols of non-zero entries
  std::vector< unsigned int > co;
} t_matCrd;

/*
 * Sparse matrix format: Compressed sparse row.
 *
 * Example:
 *    _     _  nnz:    4
 *   | 1 0 2 | rowPtr: 0,2,3,4
 *   | 0 3 0 | colIdx: 0,2,1,0
 *   |_4 0 0_| val:    1,2,3,4
 *
 */
typedef struct {
  // number of non-zero values
  unsigned int nnz;
  // non-zero values
  std::vector< real_base    > val;
  // row pointers
  std::vector< unsigned int > rowPtr;
  // column indices
  std::vector< unsigned int > colIdx;
} t_matCsr;

/*
 * Sparse matrix format: Compressed sparse column.
 *
 * Example:
 *    _     _  nnz:    4
 *   | 1 0 2 | colPtr: 0,2,3,6
 *   | 0 0 3 | rowIdx: 0,2,2,0,1,2
 *   |_4 5 6_| val:    1,4,5,2,3,6
 *
 */
typedef struct {
  // number of non-zero values
  unsigned int nnz;
  // non-zero values
  std::vector< real_base    > val;
  // column pointers
  std::vector< unsigned int > colPtr;
  // row indices
  std::vector< unsigned int > rowIdx;
} t_matCsc;

namespace edge {
  namespace linalg {
    class Matrix;
  }
}

class edge::linalg::Matrix {
  public:
    /**
     * Computes the derminant of 2x2 matrix.
     *
     * @param i_mat matrix which determinant we are computing.
     * @return determinant of the matrix.
     **/
    template <typename T>
    static T det2x2( const T i_mat[2][2] ) {
      return i_mat[0][0]*i_mat[1][1] - i_mat[0][1]*i_mat[1][0];
    }

    /**
     * Computes the derminant of 3x3 matrix.
     *
     * @param i_mat matrix which determinant we are computing.
     * @return determinant of the matrix.
     **/
    template <typename T>
    static T det3x3( const T i_mat[3][3] ) {
      return   i_mat[0][0] * (i_mat[1][1]*i_mat[2][2] - i_mat[1][2]*i_mat[2][1])
             - i_mat[0][1] * (i_mat[1][0]*i_mat[2][2] - i_mat[1][2]*i_mat[2][0])
             + i_mat[0][2] * (i_mat[1][0]*i_mat[2][1] - i_mat[1][1]*i_mat[2][0]);
    }

    /**
     * Computes the determinant of 4x4 matrix by expanding the last column.
     *
     * @param i_mat matrix which determinant we are computing.
     * @return determinant of the matrix.
     **/
    template <typename T>
    static T det4x4( const T i_mat[4][4] ) {
      // determinant of the 4x4
      T l_det4x4 = 0;

      for( unsigned short l_expRo = 0; l_expRo < 4; l_expRo++ ) {
        // minor matrix
        T l_minor[3][3];

        // column in the minor
        unsigned short l_minRo = 0;

        // get minor
        for( unsigned short l_ro = 0; l_ro < 4; l_ro++ ) {
          if( l_ro != l_expRo ) {
            for( unsigned short l_co = 0; l_co < 3; l_co++ ) {
              l_minor[l_minRo][l_co] = i_mat[l_ro][l_co];
            }
            l_minRo++;
          }
        }

        // add effect of minor determinant
        T l_det3x3 = det3x3( l_minor );
        if( l_expRo % 2 == 0 ) l_det4x4 -= i_mat[l_expRo][3] * l_det3x3;
        else                   l_det4x4 += i_mat[l_expRo][3] * l_det3x3;
      }

      return l_det4x4;
    }

    /**
     * Computes the derminant of 2x2 or 3x3 or 4x4 matrix.
     *
     * @parma i_dim number of rows and columns in the matrix.
     * @param i_mat matrix which determinant we are computing.
     * @return determinant of the matrix.
     **/
    template <typename T>
    static T det( const unsigned short  i_dim,
                  const T              *i_mat ) {
      T l_det = 0;

      if(      i_dim == 2 ) l_det = det2x2(  (T(*)[2]) i_mat );
      else if( i_dim == 3 ) l_det = det3x3(  (T(*)[3]) i_mat );
      else if( i_dim == 4 ) l_det = det4x4(  (T(*)[4]) i_mat );
      else EDGE_LOG_FATAL;

      return l_det;
    }

    /**
     * Computes the inverse of a 2x2 matrix
     *
     * @param i_mat matrix which inverse is computed
     * @param o_inv will be set to the inverse of the matrix
     **/
    template <typename T>
    static void inv2x2( const T i_mat[2][2],
                              T o_inv[2][2] ) {
      T l_det = det2x2( i_mat );
      assert( l_det > TOL.LINALG );

      o_inv[0][0] =  i_mat[1][1];
      o_inv[0][1] = -i_mat[0][1];
      o_inv[1][0] = -i_mat[1][0];
      o_inv[1][1] =  i_mat[0][0];

      for( unsigned short l_ro = 0; l_ro < 2; l_ro++ ) {
        for( unsigned short l_co = 0; l_co < 2; l_co++ ) {
          o_inv[l_ro][l_co] /= l_det;
        }
      }
    }

    /**
     * Computes the inverse of a 3x3 matrix
     *
     * @param i_mat matrix which inverse is computed
     * @param o_inv will be set to the inverse of the matrix
     **/
    template <typename T>
    static void inv3x3( const T i_mat[3][3],
                              T o_inv[3][3] ) {
      // get the determinant
      T l_det = det3x3( i_mat );
      assert( std::abs(l_det) > TOL.LINALG );

      // compute the inverse
      o_inv[0][0] = i_mat[1][1]*i_mat[2][2] - i_mat[1][2]*i_mat[2][1];
      o_inv[0][1] = i_mat[0][2]*i_mat[2][1] - i_mat[0][1]*i_mat[2][2];
      o_inv[0][2] = i_mat[0][1]*i_mat[1][2] - i_mat[0][2]*i_mat[1][1];

      o_inv[1][0] = i_mat[1][2]*i_mat[2][0] - i_mat[1][0]*i_mat[2][2];
      o_inv[1][1] = i_mat[0][0]*i_mat[2][2] - i_mat[0][2]*i_mat[2][0];
      o_inv[1][2] = i_mat[0][2]*i_mat[1][0] - i_mat[0][0]*i_mat[1][2];

      o_inv[2][0] = i_mat[1][0]*i_mat[2][1] - i_mat[1][1]*i_mat[2][0];
      o_inv[2][1] = i_mat[0][1]*i_mat[2][0] - i_mat[0][0]*i_mat[2][1];
      o_inv[2][2] = i_mat[0][0]*i_mat[1][1] - i_mat[0][1]*i_mat[1][0];

      for( unsigned short l_ro = 0; l_ro < 3; l_ro++ ) {
        for( unsigned short l_co = 0; l_co < 3; l_co++ ) {
          o_inv[l_ro][l_co] /= l_det;
        }
      }
    }

    /**
     * Performs the operation C = A.B.
     *
     * @param i_m blas identifier M.
     * @param i_n blas identifier N.
     * @param i_k blas identifier K.
     * @param i_a matrix A.
     * @param i_b matrix B.
     * @param o_c matrix C, will bet set to A.B.
     *
     * @paramt TL_T_REAL_A precision of matrix A.
     * @paramt TL_T_REAL_B precision of matrix B.
     * @paramt TL_T_REAL_C precision of matrix C.
     **/
    template <typename TL_T_REAL_A, typename TL_T_REAL_B, typename TL_T_REAL_C>
    static void matMulB0(       unsigned int  i_m,
                                unsigned int  i_n,
                                unsigned int  i_k,
                          const TL_T_REAL_A  *i_a,
                          const TL_T_REAL_B  *i_b,
                                TL_T_REAL_C  *o_c ) {
      // reset result to zero
      for( unsigned int l_m = 0; l_m < i_m; l_m++ ) {
        for( unsigned int l_n = 0; l_n < i_n; l_n++ ) {
            o_c[l_m * i_n + l_n] = 0;
        }
      }

      for( unsigned int l_k = 0; l_k < i_k; l_k++ ) {
        for( unsigned int l_m = 0; l_m < i_m; l_m++ ) {
          for( unsigned int l_n = 0; l_n < i_n; l_n++ ) {
            o_c[l_m*i_n + l_n] += i_a[l_m*i_k +l_k] * i_b[l_k*i_n+l_n];
          }
        }
      }
    }

    /**
     * Performs the operation C[r] = A.B[r], for all 0 =< r =< #matrices.
     * Here, B and C are arrays of matrices, with 0 =< r =< #matrices being the fastest dimensions.
     *
     * @param i_r number of B and C matrices.
     * @param i_m blas identifier M.
     * @param i_n blas identifier N.
     * @param i_k blas identifier K.
     * @param i_a matrix A.
     * @param i_b matrix B.
     * @param o_c matrix C, will bet set to A.B.
     *
     * @paramt TL_T_REAL_A precision of matrix A.
     * @paramt TL_T_REAL_B precision of matrix B.
     * @paramt TL_T_REAL_C precision of matrix C.
     **/
    template <typename TL_T_REAL_A, typename TL_T_REAL_B, typename TL_T_REAL_C>
    static void matMulB0FusedBC(       unsigned short  i_r,
                                       unsigned int    i_m,
                                       unsigned int    i_n,
                                       unsigned int    i_k,
                                 const TL_T_REAL_A    *i_a,
                                 const TL_T_REAL_B    *i_b,
                                       TL_T_REAL_C    *o_c ) {
      // reset result to zero
      for( unsigned int l_m = 0; l_m < i_m; l_m++ ) {
        for( unsigned int l_n = 0; l_n < i_n; l_n++ ) {
          for( unsigned short l_r = 0; l_r < i_r; l_r++ ) {
            o_c[l_m*i_n*i_r + l_n*i_r + l_r] = 0;
          }
        }
      }

      for( unsigned int l_k = 0; l_k < i_k; l_k++ ) {
        for( unsigned int l_m = 0; l_m < i_m; l_m++ ) {
          for( unsigned int l_n = 0; l_n < i_n; l_n++ ) {
            for( unsigned short l_r = 0; l_r < i_r; l_r++ ) {
              o_c[l_m*i_n*i_r + l_n*i_r + l_r] += i_a[l_m*i_k + l_k] * i_b[l_k*i_n*i_r + l_n*i_r + l_r];
            }
          }
        }
      }
    }

    /**
     * Performs the operation C[r] = A[r].B, for all 0 =< r =< #matrices.
     * Here, A and C are arrays of matrices, with 0 =< r =< #matrices being the fastest dimensions.
     *
     * @param i_r number of A and C matrices.
     * @param i_m blas identifier M.
     * @param i_n blas identifier N.
     * @param i_k blas identifier K.
     * @param i_a matrix A.
     * @param i_b matrix B.
     * @param o_c matrix C, will bet set to A.B.
     *
     * @paramt TL_T_REAL_A precision of matrix A.
     * @paramt TL_T_REAL_B precision of matrix B.
     * @paramt TL_T_REAL_C precision of matrix C.
     **/
    template <typename TL_T_REAL_A, typename TL_T_REAL_B, typename TL_T_REAL_C>
    static void matMulB0FusedAC(       unsigned short  i_r,
                                       unsigned int    i_m,
                                       unsigned int    i_n,
                                       unsigned int    i_k,
                                 const TL_T_REAL_A    *i_a,
                                 const TL_T_REAL_B    *i_b,
                                       TL_T_REAL_C    *o_c ) {
      // reset result to zero
      for( unsigned int l_m = 0; l_m < i_m; l_m++ ) {
        for( unsigned int l_n = 0; l_n < i_n; l_n++ ) {
          for( unsigned short l_r = 0; l_r < i_r; l_r++ ) {
            o_c[l_m*i_n*i_r + l_n*i_r + l_r] = 0;
          }
        }
      }

      for( unsigned int l_k = 0; l_k < i_k; l_k++ ) {
        for( unsigned int l_m = 0; l_m < i_m; l_m++ ) {
          for( unsigned int l_n = 0; l_n < i_n; l_n++ ) {
            for( unsigned short l_r = 0; l_r < i_r; l_r++ ) {
              o_c[l_m*i_n*i_r + l_n*i_r + l_r] += i_a[l_m*i_k*i_r + l_k*i_r + l_r] * i_b[l_k*i_n + l_n];
            }
          }
        }
      }
    }

    /**
     * Performs the operation C[r] += A.B[r], for all 0 =< r =< #matrices.
     * Here, B and C are arrays of matrices, with 0 =< r =< #matrices being the fastest dimensions.
     *
     * @param i_r number of B and C matrices.
     * @param i_m blas identifier M.
     * @param i_n blas identifier N.
     * @param i_k blas identifier K.
     * @param i_a matrix A.
     * @param i_b matrix B.
     * @param o_c matrix C, will bet set to A.B.
     *
     * @paramt TL_T_REAL_A precision of matrix A.
     * @paramt TL_T_REAL_B precision of matrix B.
     * @paramt TL_T_REAL_C precision of matrix C.
     **/
    template <typename TL_T_REAL_A, typename TL_T_REAL_B, typename TL_T_REAL_C>
    static void matMulB1FusedBC(       unsigned short  i_r,
                                       unsigned int    i_m,
                                       unsigned int    i_n,
                                       unsigned int    i_k,
                                 const TL_T_REAL_A    *i_a,
                                 const TL_T_REAL_B    *i_b,
                                       TL_T_REAL_C    *o_c ) {
      for( unsigned int l_k = 0; l_k < i_k; l_k++ ) {
        for( unsigned int l_m = 0; l_m < i_m; l_m++ ) {
          for( unsigned int l_n = 0; l_n < i_n; l_n++ ) {
            for( unsigned short l_r = 0; l_r < i_r; l_r++ ) {
              o_c[l_m*i_n*i_r + l_n*i_r + l_r] += i_a[l_m*i_k + l_k] * i_b[l_k*i_n*i_r + l_n*i_r + l_r];
            }
          }
        }
      }
    }

    /**
     * Transposes the given dense matrix.
     *
     * @param i_nModes number of modes; equals number of columns and rows.
     * @param io_matrix matrix which is transposed.
     **/
    static void transposeDense( int_md     i_nModes,
                                real_base *io_matrix ) {
      for( int_md l_ro = 0; l_ro < i_nModes; l_ro++ ) {
        for( int_md l_co = l_ro+1; l_co < i_nModes; l_co++ ) {
          unsigned int l_oId = l_ro * i_nModes + l_co;
          unsigned int l_nId = l_co * i_nModes + l_ro;

          real_base l_tmp  = io_matrix[l_oId];
          io_matrix[l_oId] = io_matrix[l_nId];
          io_matrix[l_nId] = l_tmp;
        }
      }
    }

    /**
     * Transposes the given matrix in coordinate format.
     *   Remark: The resulting matrix is store in row major again.
     *
     * @param i_mat matrix which gets transposed.
     * @param o_matT will be set to resulting, transposed matrix.
     **/
    static void transposeCrd( const t_matCrd &i_mat,
                                    t_matCrd &o_matT ) {
      // derive col and row size, at least of those which are non-zero
      unsigned int l_nColsOld = 0;
      unsigned int l_nRowsOld = 0;
      for( unsigned int l_nzId = 0; l_nzId < i_mat.nz.size(); l_nzId++ ) {
        l_nRowsOld = std::max( l_nRowsOld, i_mat.ro[l_nzId] );
        l_nColsOld = std::max( l_nColsOld, i_mat.co[l_nzId] );
      }

      // clear new matrix
      o_matT.ro.clear();
      o_matT.co.clear();
      o_matT.nz.clear();

      // do the transpose
      for( unsigned int l_rowNew = 0; l_rowNew <= l_nColsOld; l_rowNew++ ) {
        for( unsigned int l_colNew = 0; l_colNew <= l_nRowsOld; l_colNew++ ) {
          for( unsigned int l_nzId = 0; l_nzId < i_mat.nz.size(); l_nzId++ ) {
            if( i_mat.ro[l_nzId] == l_colNew && i_mat.co[l_nzId] == l_rowNew ) {
              o_matT.ro.push_back( l_rowNew );
              o_matT.co.push_back( l_colNew );
              o_matT.nz.push_back( i_mat.nz[l_nzId] );
            }
          }
        }
      }
    }

    /**
     * Prints a dense matrix with the given #rows and #cols.
     *
     * @param i_nRows number of rows.
     * @param i_nCols number of cols.
     * @param i_mat dense matrix.
     **/
    static void printMatrixDense(       unsigned i_nRows,
                                        unsigned i_nCols,
                                  const real_base *i_mat ) {
      // stream for the outout
      std::stringstream l_stream;

      for( unsigned int l_ro = 0; l_ro < i_nRows; l_ro++ ) {
        for( unsigned int l_co = 0; l_co < i_nCols; l_co++ ) {
          l_stream << std::fixed << std::setprecision(2) << i_mat[l_ro*i_nCols + l_co] << "\t";
        }
        EDGE_LOG_INFO << "    " << l_stream.str();
        l_stream.str("");
      }
    }

    /**
     * Prints a given matrix in coordinate format as dense.
     *
     * @param i_nRows number rows.
     * @param i_nCols number of cols.
     * @param i_mat matrix which is printed.
     **/
    static void printMatrixCrd(       unsigned int  i_nRows,
                                      unsigned int  i_nCols,
                                const t_matCrd     &i_mat ) {
      // stream for the outout
      std::stringstream l_stream;

      // id of the current nz entry
      unsigned int l_nzId = 0;

      for( unsigned int l_ro = 0; l_ro < i_nRows; l_ro++ ) {
        for( unsigned int l_co = 0; l_co < i_nCols; l_co++ ) {
          if( i_mat.nz.size() > l_nzId && (i_mat.co[l_nzId] == l_co && i_mat.ro[l_nzId] == l_ro) ) {
            l_stream << std::fixed << std::setprecision(2) << i_mat.nz[l_nzId] << "\t";
            l_nzId++;
          }
          else l_stream << "----\t";
        }
        EDGE_LOG_INFO << "    " << l_stream.str();
        l_stream.str("");
      }
    }

    /**
     * Prints a given matrix in CSR-format.
     *
     * @param i_nRows number of rows in the matrix.
     * @param i_mat matrix which will bet printed.
     **/
    static void printMatrixCsr(        unsigned int  i_nRows,
                                const t_matCsr     &i_mat ) {
      // get number of non-zeros
      unsigned int l_nnz = i_mat.rowPtr[i_nRows];

      // stream for the outout
      std::stringstream l_stream;
      l_stream << std::fixed << std::setprecision(2); 

      l_stream << "vals:   ";
      for( unsigned int l_nzId = 0; l_nzId < l_nnz; l_nzId++ ) {
        l_stream << i_mat.val[l_nzId] << " ";
      }
      EDGE_LOG_INFO << "    " << l_stream.str();
      l_stream.str("");

      l_stream << "cols:   ";
      for( unsigned int l_nzId = 0; l_nzId < l_nnz; l_nzId++ ) {
        l_stream << i_mat.colIdx[l_nzId] << " ";
      }
      EDGE_LOG_INFO << "    " << l_stream.str();
      l_stream.str("");

      l_stream << "rowptr: ";
      for( unsigned int l_ro = 0; l_ro <= i_nRows; l_ro++ ) {
        l_stream << i_mat.rowPtr[l_ro] << " ";
      }
      EDGE_LOG_INFO << "    " << l_stream.str();
    }

    /**
     * Gets the first nonzero-row, last nonzero-row,
     * first nonzero-column and last nonzero-column in the given matrix.
     *
     * Example:
     *
     *  x: non-zero entry
     *    __        __
     *  0| x          |   A call with i_nZBlock = { { 0, 2 }, {-1, -1 } }
     *  1|    x       |   will return o_nZBlock = { { 0, 2 }, { 0,  1 } }
     *  2| x          |
     *  3|    x     x |   A call with i_nZBlock = { {-1, -1}, { 1,  2 } } 
     *  4|__        __|   will return o_nZBlock = { { 1,  3}, { 1,  1 } }
     *     0  1  2  3
     *
     * @param i_mat matrix which nonzero-block is determined.
     * @param i_BlockNz contains the non-zero block considered in the derivation [0][0]: first row, [0][1]: last row, [1][0]: first col, [1][1]: last col; std::numeric_limits<unsigned int>::max() if unbound and the entire matrix dimension is considered.
     * @param o_BlockNz will be set to non-zero block [0][0]: first row, [0][1]: last row, [1][0]: first col, [1][1]: last col; std::numeric_limits<unsigned int>::max() will be set if no entries exist.
     **/
    static void getBlockNz( const t_matCrd     &i_mat,
                            const unsigned int  i_blockNz[2][2],
                                  unsigned int  o_blockNz[2][2] ) {
      // reset output
      o_blockNz[0][0] = o_blockNz[0][1] = o_blockNz[1][0] = o_blockNz[1][1] = std::numeric_limits<unsigned int>::max();

      // iterate over non-zero entries
      for( unsigned int l_nz = 0; l_nz < i_mat.nz.size(); l_nz++ ) {
        // check if the entry is in requested bounds
        bool l_inBo = true;
        l_inBo = l_inBo && (    i_mat.ro[l_nz] >= i_blockNz[0][0]
                             || i_blockNz[0][0] == std::numeric_limits<unsigned int>::max() );
        l_inBo = l_inBo && (    i_mat.ro[l_nz] <= i_blockNz[0][1]
                             || i_blockNz[0][1] == std::numeric_limits<unsigned int>::max() );
        l_inBo = l_inBo && (    i_mat.co[l_nz] >= i_blockNz[1][0]
                             || i_blockNz[1][0] == std::numeric_limits<unsigned int>::max() );
        l_inBo = l_inBo && (    i_mat.co[l_nz] <= i_blockNz[1][1]
                             || i_blockNz[1][1] == std::numeric_limits<unsigned int>::max() );

        // check for new extremum if we are in bounds
        if( l_inBo ) {
          if( o_blockNz[0][0] == std::numeric_limits<unsigned int>::max() || o_blockNz[0][0] > i_mat.ro[l_nz] )
            o_blockNz[0][0] = i_mat.ro[l_nz];

          if( o_blockNz[0][1] == std::numeric_limits<unsigned int>::max() || o_blockNz[0][1] < i_mat.ro[l_nz] )
            o_blockNz[0][1] = i_mat.ro[l_nz];

          if( o_blockNz[1][0] == std::numeric_limits<unsigned int>::max() || o_blockNz[1][0] > i_mat.co[l_nz] )
            o_blockNz[1][0] = i_mat.co[l_nz];

          if( o_blockNz[1][1] == std::numeric_limits<unsigned int>::max() || o_blockNz[1][1] < i_mat.co[l_nz] )
            o_blockNz[1][1] = i_mat.co[l_nz];
        }
      }
    }

    /**
     * Gets the number of non-zeros in the submatrix of size #nrows#ncols.
     *
     * @param i_nRows number of rows in the submatrix.
     * @param i_nCols number of cols in the submatrix.
     * @param i_matrix matrix which number of non-zeros is derived.
     **/
    static unsigned int getNnzCrd(       unsigned int  i_nRows,
                                         unsigned int  i_nCols,
                                   const t_matCrd     &i_matrix ) {
      unsigned int l_nnz = 0;

      for( unsigned int l_nzId = 0; l_nzId < i_matrix.nz.size(); l_nzId++ ) {
          if( i_matrix.ro[l_nzId] < i_nRows &&
              i_matrix.co[l_nzId] < i_nCols ) {
            l_nnz++;
          }
      }

      return l_nnz;
    }

    /**
     * Converts a dense matrix in row-major format to coordinate format.
     *
     * @param i_nRows number of rows.
     * @param i_nCols number of columns.
     * @param i_a dense matrix in row-major storage which gets converted.
     * @param o_crd will be set to result in coordinate format.
     * @param i_tol tolerance/delta which is considered to be zero for the matrix entries.
     *
     **/
    template <typename T>
    static void denseToCrd(     unsigned int  i_nRows,
                                unsigned int  i_nCols,
                          const T            *i_a,
                                t_matCrd     &o_crd,
                                T             i_tol = 0.000001 ) {
      EDGE_CHECK( i_tol > 0 );

      // reset the matrix
      o_crd.nz.clear();
      o_crd.ro.clear();
      o_crd.co.clear();

      for( unsigned int l_ro = 0; l_ro < i_nRows; l_ro++ ) {
        for( unsigned int l_co = 0; l_co < i_nCols; l_co++ ) {
          // derive 1D id
          unsigned int l_id = l_ro * i_nCols + l_co;

          if( std::abs( i_a[l_id] ) > i_tol ) {
            o_crd.nz.push_back( i_a[l_id] );
            o_crd.ro.push_back( l_ro      );
            o_crd.co.push_back( l_co      );
          }
        }
      }
    }

    /**
     * Converts a matrix in coordinate format to CSR-format.
     *
     * @param i_subMatRows number of rows in the sub-matrix extracted from coord format.
     * @param i_subMatCols number of cols in the sub-matrix extracted from coord format.
     * @param i_coord matrix in coordinate format.
     * @param o_csr will be set to sub-matrix in CSR format.
     **/
    static void crdToCsr(       unsigned int  i_subMatRows,
                                unsigned int  i_subMatCols,
                          const t_matCrd     &i_coord,
                                t_matCsr     &o_csr ) {
      // get the number of non-zero entries for the given submatrix
      unsigned int l_nnz = getNnzCrd( i_subMatRows,
                                      i_subMatCols,
                                      i_coord );

      o_csr.val.resize(    l_nnz          );
      o_csr.rowPtr.resize( i_subMatRows+1 );
      o_csr.colIdx.resize( l_nnz          );

      // el gives the number of non-zero elements we have seen, equivalent to id in csr
      unsigned int l_el = 0;

      for( unsigned int l_ro = 0; l_ro < i_subMatRows; l_ro++ ) {
        o_csr.rowPtr[l_ro] = l_el;
        for( unsigned int l_co = 0; l_co < i_subMatCols; l_co++ ) {
          // iterate over non-zeros
          for( unsigned int l_nzId = 0; l_nzId < i_coord.nz.size(); l_nzId++ ) {
            if( l_ro == i_coord.ro[l_nzId] && l_co == i_coord.co[l_nzId] ) {
              o_csr.val[l_el] = i_coord.nz[l_nzId];
              o_csr.colIdx[l_el] = i_coord.co[l_nzId];
              l_el++;
            }
          }
        }
      }

      // set number of non-zeros as last entry
      o_csr.rowPtr[i_subMatRows] = l_nnz;
    }

    /**
     * Converts a matrix in coordinate format to CSC-format.
     *
     * @param i_subMatRows number of rows in the sub-matrix extracted from coord format.
     * @param i_subMatCols number of cols in the sub-matrix extracted from coord format.
     * @param i_coord matrix in coordinate format.
     * @param o_csc will be set to sub-matrix in CSC format.
     **/
    static void crdToCsc(       unsigned int  i_subMatRows,
                                unsigned int  i_subMatCols,
                          const t_matCrd     &i_coord,
                                t_matCsc     &o_csc ) {
      // get the number of non-zero entries for the given submatrix
      unsigned int l_nnz = getNnzCrd( i_subMatRows,
                                      i_subMatCols,
                                      i_coord );

      o_csc.val.resize(    l_nnz          );
      o_csc.colPtr.resize( i_subMatCols+1 );
      o_csc.rowIdx.resize( l_nnz          );

      // el gives the number of non-zero elements we have seen, equivalent to id in csc
      unsigned int l_el = 0;

      for( unsigned int l_co = 0; l_co < i_subMatCols; l_co++ ) {
        o_csc.colPtr[l_co] = l_el;

        for( unsigned int l_ro = 0; l_ro < i_subMatRows; l_ro++ ) {
          // iterate over non-zeros
          for( unsigned int l_nzId = 0; l_nzId < i_coord.nz.size(); l_nzId++ ) {
            if( l_ro == i_coord.ro[l_nzId] && l_co == i_coord.co[l_nzId] ) {
              o_csc.val[l_el] = i_coord.nz[l_nzId];
              o_csc.rowIdx[l_el] = i_coord.ro[l_nzId];
              l_el++;
            }
          }
        }
      }

      // set number of non-zeros as last entry
      o_csc.colPtr[i_subMatCols] = l_nnz;
    }

    /**
     * Converts a dense matrix in row-major format to compressed-sparse-row format.
     *
     * @param i_nRows number of rows.
     * @param i_nCols number of columns.
     * @param i_a dense matrix in row-major storage which gets converted.
     * @param o_csr will be set to result in csr-format.
     * @param i_tol tolerance/delta which is considered to be zero for the matrix entries.
     * @param i_subMatRows numeber of rows in the sub-matrix extracted.
     * @param i_subMatCols numeber of cols in the sub-matrix extracted.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template< typename TL_T_REAL >
    static void denseToCsr( unsigned int        i_nRows,
                            unsigned int        i_nCols,
                            TL_T_REAL    const *i_a,
                            t_matCsr           &o_csr,
                            TL_T_REAL           i_tol = 0.000001,
                            unsigned int        i_subMatRows = std::numeric_limits< unsigned int >::max(),
                            unsigned int        i_subMatCols = std::numeric_limits< unsigned int >::max() ) {
      // temporary coord matrix
      t_matCrd l_tmpCrd;

      // convert to coordinate format
      denseToCrd< TL_T_REAL >( i_nRows, i_nCols, i_a, l_tmpCrd, i_tol );

      // adjust submatrix size if not defined
      unsigned int l_nRows = i_subMatRows;
      unsigned int l_nCols = i_subMatCols;
      if( l_nRows == std::numeric_limits< unsigned int >::max() ) l_nRows = i_nRows;
      if( l_nCols == std::numeric_limits< unsigned int >::max() ) l_nCols = i_nCols;

      // do the conversion to csr
      crdToCsr( l_nRows, l_nCols, l_tmpCrd, o_csr );
    }

    /**
     * Converts a dense matrix in row-major format to compressed-sparse-column format.
     *
     * @param i_nRows number of rows.
     * @param i_nCols number of columns.
     * @param i_a dense matrix in row-major storage which gets converted.
     * @param o_csc will be set to result in csc-format.
     * @param i_tol tolerance/delta which is considered to be zero for the matrix entries.
     * @param i_subMatRows numeber of rows in the sub-matrix extracted.
     * @param i_subMatCols numeber of cols in the sub-matrix extracted.
     *
     * @paramt TL_T_REAL floating point precision.
     **/
    template <typename TL_T_REAL>
    static void denseToCsc( unsigned int        i_nRows,
                            unsigned int        i_nCols,
                            TL_T_REAL    const *i_a,
                            t_matCsc           &o_csc,
                            TL_T_REAL           i_tol = 0.000001,
                            unsigned int        i_subMatRows = std::numeric_limits< unsigned int >::max(),
                            unsigned int        i_subMatCols = std::numeric_limits< unsigned int >::max() ) {
      // temporary coord matrix
      t_matCrd l_tmpCrd;

      // convert to coordinate format
      denseToCrd< TL_T_REAL >( i_nRows, i_nCols, i_a, l_tmpCrd, i_tol );

      // adjust submatrix size if not defined
      unsigned int l_nRows = i_subMatRows;
      unsigned int l_nCols = i_subMatCols;
      if( l_nRows == std::numeric_limits< unsigned int >::max() ) l_nRows = i_nRows;
      if( l_nCols == std::numeric_limits< unsigned int >::max() ) l_nCols = i_nCols;

      // do the conversion to csc
      crdToCsc( l_nRows, l_nCols, l_tmpCrd, o_csc );
    }
};

#endif
