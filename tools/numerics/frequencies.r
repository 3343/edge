##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2016, Regents of the University of California
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
# Analysis of frequency contant from the plane waves convergence benchmark.
##

# name of the error norm and var
error_norm <- 'sigma_yz_linf'

# upper limit of the error to search for
error_limit <- 0.001

# benchmark output
file <- 'all.csv'

# domain size 1d
domain_size <- 100

# number of waves w.r.t to the diagonal
n_waves <- 6

# s-wave velocity
vel <- 1

# diagonal length
diagonal_length <- sqrt(3)*domain_size

# derive wave length and requency
wave_length <- diagonal_length / n_waves
freq <- vel / wave_length

# read benchmarks
benchmarks <- read.csv(file)

# compute ratio of wave length / charateristic length (equals cube length)
benchmarks$cl_ratio <- wave_length / ( domain_size / benchmarks$cubes_per_dim )

# remove single precision
subset <- c(benchmarks$precision == 'd')
benchmarks <- benchmarks[subset,]

# extract indices matching our criterion
subset <- benchmarks[error_norm] < error_limit
indices <- as.numeric( rownames( unique( benchmarks[subset,1:2]) ) )

# print the results
benchmarks[indices,]
