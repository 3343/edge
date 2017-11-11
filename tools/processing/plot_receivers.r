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
# Receiver plots.
##

##
# Plots the given quantity of the receiver.
#
# @param data data of the receiver.
# @param q quantity which gets plotted.
##
plotRecv <- function( data, q ) {
  layout( matrix( c(1,2), 2, 1, byrow=TRUE) )
  
  min = toString( min(data[,q]) )
  max = toString( max(data[,q]) )
  
  plot( data[, 'time'], data[, q], type='l', xlab='time (s)', ylab=q,
        main = paste( 'min:', min, '\nmax:', max  )
        )

  spec <- spectrum( data[,q], log='no', plot=FALSE )
  spf <- spec$freq * ( length(data[,q]) / ( max(data[,'time']) - min(data[,'time'])) )
  sps <- spec$spec

  plot( spf, sps, type='l', xlim=c(0,20), xlab='frequency (Hz)', ylab='spectrum' )
}

# get command line arguments
cArgs <- commandArgs()

if( grepl( '-h', cArgs[6]) ) {
  print('Plots receivers and the spectrum. Expects .csv-files in the inout-directory and writes to same paths but with .pdf file extension.')
  print('Use via: Rscript plot_receivers.r  INOUT_DIRECTORY')
  stop()
}
if( length(cArgs) != 6 ) {
  stop( 'Error: Script requires exactly one argument: Rscript plot_receivers.r  INOUT_DIRECTORY'  )
}

dir <- cArgs[6]

# get the receivers
recvs <- list.files( dir, pattern = ".csv" )

for(  fileIn in recvs  ) {
  print( paste('processing receiver data', fileIn) )
  
  # read receiver data
  fileIn  <- paste( dir, fileIn, sep='/')
  fileOut <- paste0( strsplit( fileIn, split='.csv' )[[1]], '.pdf' )
  recv <- read.csv(fileIn)
  
  # process receiver
  qs <- names(recv)
  
  pdf( file=fileOut  )
  for( q in qs[2:length(qs)]  ) {
    plotRecv( recv, q )
  }
  dev.off()
}
