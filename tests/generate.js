#!/usr/bin/env node
////
// @file This file is part of EDGE.
//
// @author Alexander Breuer (breuer AT mytum.de)
//
// @section LICENSE
// Copyright (c) 2020, Alexander Breuer
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// @section DESCRIPTION
// Generates the build configs.
////
'use strict';

var l_nunjucks = require('nunjucks');
var l_fs = require('fs');
const { exit } = require('process');

var l_args = require('minimist')( process.argv.slice(2), {
  default: { isa: 'avx2',
             cfs: [1, 16],
             eqs: ['elastic', 'viscoelastic2', 'viscoelastic3', 'viscoelastic4', 'viscoelastic5'],
             ets: ['tet4', 'tria3'],
             ors: [1, 2, 3, 4, 5],
             prs: [32, 64],
             up:  '' }
});

// convert potential single-entry command line input to arrays
for( var l_key of ['cfs', 'eqs', 'ets', 'ors', 'prs', 'up'] ) {
  l_args[l_key] = Array.isArray( l_args[l_key] ) ? l_args[l_key] : [l_args[l_key]];
}

// assemble builds
var l_builds = []

for( var l_cf of l_args.cfs ) {
  for( var l_eq of l_args.eqs ) {
    for( var l_et of l_args.ets ) {
      for( var l_or of l_args.ors ) {
        for( var l_pr of l_args.prs ) {
          l_builds.push( {
            'cfr':          l_cf,
            'equations':    l_eq,
            'element_type': l_et,
            'order':        l_or,
            'precision':    l_pr,
            'upload':       l_args.up
          } )
        }
      }
    }
  }
}

console.log( l_nunjucks.render( 'build.njk', {
                                  'i_isa': l_args.isa,
                                  'i_builds': l_builds
                                }
                              )
)
