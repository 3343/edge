#!/usr/bin/env Rscript
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2019, Alexander Breuer
# Copyright (c) 2016-2017, Regents of the University of California
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
# Converts the csv-receiver data to plain receivers (1 value per line).
##

# get command line arguments
cArgs <- commandArgs()

if( grepl( '-h', cArgs[6]) ) {
  print('Converts receivers in csv formart to simple text-format. N_FREQ only extract evert N_FREQth value.')
  print('If DIM_TIME is set, the time dimension will be added as first column.')
  print('Expects .csv-files in the inout-directory and writes to same paths but with _QUANTITIY_CRUN.dat file extension.')
  print('Use via: Rscript recvs_csv_to_plain.r  INOUT_DIRECTORY N_FREQ DIM_TIME')
  stop()
}
if( length(cArgs) < 6 ) {
  stop( 'Error: Script requires at least one argument: Rscript recvs_csv_to_plain.r  INOUT_DIRECTORY N_FREQ'  )
}

dir <- cArgs[6]

if( length(cArgs) > 6) {
  freq <- as.integer(cArgs[7])
} else {
  freq <- 1
}

if( length(cArgs) > 7) {
  dimTime <- as.integer(cArgs[8])
} else {
  dimTime <- 0
}

# get the receivers
recvs <- list.files( dir, pattern = ".csv" )

for(  fileIn in recvs  ) {
  print( paste('processing receiver data', fileIn) )
  
  # read receiver data
  fileIn  <- paste( dir, fileIn, sep='/')
  recv <- read.csv(fileIn, comment.char='#')
  
  # process receiver
  qs <- names(recv)
  
  # derive subset of values (int-select)
  subSet <- seq(1,length(recv[, 'time']),freq)
  # derive numeric frequency
  freqNumeric <- recv[,'time'][subSet[2]]-recv[,'time'][1]

  for( q in qs[2:length(qs)]  ) {
    fileOut <- paste0( strsplit( fileIn, split='.csv' )[[1]],'_',q,'.dat' )

    header <- paste( length(subSet), freqNumeric, 1 )
    # write header
    write( header, file=fileOut, sep="\n" )
    
    # extract data
    if( dimTime == 0 ) {
      outData <- recv[,q][subSet]
    } else {
      outData <- recv[c('time',q) ]
      outData <- outData[subSet,]
    }

    # write to disk
    write.table( outData, file=fileOut, row.names=FALSE, col.names=FALSE, append=TRUE )
  }
}
