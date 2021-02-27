#!/usr/bin/env python3
##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2021, Friedrich Schiller University Jena
# Copyright (c) 2019, Alexander Breuer
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
# Convergence plots.
##
import logging
import argparse
import glob
import os
import re
import xmltodict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
from matplotlib.backends.backend_pdf import PdfPages

# set up logger
logging.basicConfig( level=logging.INFO,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )
# command line arugments
l_parser    = argparse.ArgumentParser( description='Produces convergence plots.' )
l_parser.add_argument( '--xdir',
                       dest     = "xdir",
                       required = True,
                       help     = "Directory where the xml files containing the error norms are located.",
                       metavar  = "XML_DIR" )

l_parser.add_argument( '--xregexp',
                       dest     = "xregexp",
                       required = True,
                       help     = "Regular expression describing the format of the xml file-names.\n\
                                   Tags ORDER_TAG and REFINEMENT_TAG are required to distinct\n\
                                   between the runs.\n\
                                   Example: errors_ORDER_TAG_REFINEMENT_TAG.xml",
                       metavar  = "XML_REGEXP" )

l_parser.add_argument( '--pdf',
                       dest     = "pdf",
                       default  = "",
                       help     = "Produces a pdf and saves it to the given location.",
                       metavar  = "OUT_FILE_PDF" )


l_parser.add_argument( '--legend',
                       dest     = "legend",
                       type     = int,
                       default  = 1,
                       help     = "Shows a legend for the plot." )

l_parser.add_argument( '--qfilter',
                       dest     = "qfilter",
                       type     = str,
                       default  = '',
                       help     = "Show only the given quantities; comma-separated list, e.g. --qfilter=0,4,3." )

l_parser.add_argument( '--cfrfilter',
                       dest     = "cfrfilter",
                       type     = str,
                       default  = '',
                       help     = "Show only the given simulations; comma-separated list, e.g. --cfrfilter=0,4,3." )

l_arguments = vars(l_parser.parse_args())

# parse filters
if( l_arguments['qfilter'] != '' ):
  l_arguments['qfilter'] = list( map( int, l_arguments['qfilter'].split(',') ) )
else:
  l_arguments['qfilter'] = range(0, 100)

if( l_arguments['cfrfilter'] != '' ):
  l_arguments['cfrfilter'] = list( map( int, l_arguments['cfrfilter'].split(',') ) )
else:
  l_arguments['cfrfilter'] = range(0, 100)

logging.info( "reading files" )

# get the files
l_files = glob.glob( l_arguments['xdir']+'/'+l_arguments['xregexp'].replace("ORDER_TAG", "*").replace("REFINEMENT_TAG", "*") )

# set up regular epressions
l_regExp = { 'order':     l_arguments['xregexp'].replace( "ORDER_TAG", "([0-9]+)" ).replace( "REFINEMENT_TAG", "(\d*[.])?\d+"   ),
             'ref':       l_arguments['xregexp'].replace( "ORDER_TAG", "[0-9]+"    ).replace( "REFINEMENT_TAG", "((\d*[.])?\d+)" ) }

# dictionary for all results
l_simRes = {}

# extract results
for l_file in range(len(l_files)):
  l_order = re.search(l_regExp['order'], l_files[l_file]).group(1)
  l_ref   = re.search(l_regExp['ref'], l_files[l_file]).group(1)

  # expand dict if required
  if not l_order in l_simRes:
    l_simRes[l_order] = {}
  if not l_ref in l_simRes[l_order]:
    l_simRes[l_order][l_ref] = {}

  # read the results
  with open(l_files[l_file], 'rb') as l_xmlFile:
    l_simRes[l_order][l_ref] = xmltodict.parse(l_xmlFile)

# get number of quantities and cfrs
l_firstKey = [ list(l_simRes.keys())[0] ]
l_firstKey = l_firstKey + [ list(l_simRes[l_firstKey[0]].keys())[0] ]

l_nQs = len(l_simRes[l_firstKey[0]][l_firstKey[1]]['error_norms']['l1']['q'])

l_tmp = []
if l_nQs > 1:
  l_tmp.extend( [l_simRes[l_firstKey[0]][l_firstKey[1]]['error_norms']['l1']['q'][0]['cfr']] )
  if( type(l_tmp[0]) is list):
    l_nCfr = len( l_tmp[0] )
  else:
    l_nCfr = 1
else:
  l_tmp.extend( [l_simRes[l_firstKey[0]][l_firstKey[1]]['error_norms']['l1']['q']['cfr']] )

if( type(l_tmp[0]) is list):
  l_nCfr = len( l_tmp[0] )
else:
  l_nCfr = 1

logging.info( "using "+str(l_nQs)+" quantities and "+str(l_nCfr)+" runs" )

l_xVals = l_simRes[l_firstKey[0]].keys()

# get the results
l_char = {}
l_yVals = {}

for l_enorm in ['l1', 'l2', 'linf']:
  l_char[l_enorm] = []
  l_yVals[l_enorm] = []
  for l_order in l_simRes.keys():
    for l_q in range(l_nQs):
      if( l_q not in l_arguments['qfilter'] ):
        continue
      for l_cfr in range(l_nCfr):
        if( l_cfr not in l_arguments['cfrfilter'] ):
          continue
        # get y vals
        l_yTmp = []

        for l_x in l_xVals:
          if l_x in l_simRes[l_order]:
            # get quantities
            l_qVals = l_simRes[l_order][l_x]['error_norms'][l_enorm]['q']

            if( l_nQs > 1 ):
              if( l_nCfr > 1 ):
                l_yTmp = l_yTmp + [ l_qVals[l_q]['cfr'][l_cfr] ]
              else:
                l_yTmp = l_yTmp + [ l_qVals[l_q]['cfr'] ]
            else:
              if( l_nCfr > 1 ):
                l_yTmp = l_yTmp + [ l_qVals['cfr'][l_cfr] ]
              else:
                l_yTmp = l_yTmp + [ l_qVals['cfr'] ]
          else:
            l_yTmp += ['0']

        l_yVals[l_enorm] = l_yVals[l_enorm] + [ l_yTmp ]
        l_char[l_enorm]  = l_char[l_enorm] + [ {'order': int(l_order), 'q': l_q+1, 'cfr': l_cfr+1 } ]

# sort the results
for l_enorm in ['l1', 'l2', 'linf']:
  for l_config in range(len(l_yVals[l_enorm])):
    l_xTmp = map(float, l_xVals    )
    l_yTmp = map(float, l_yVals[l_enorm][l_config] )

    l_tmp = sorted( zip( l_xTmp, l_yTmp ) )
    l_tmp, l_yVals[l_enorm][l_config] = zip(*l_tmp)

l_xVals = sorted( map(float, l_xVals) )

# plot the results
l_style = { 'markers': {  1: '.',
                          2: '^',
                          3: 'o',
                          4: '*',
                          5: 'p', 
                          6: 's',
                          7: 'v',
                          8: '+',
                          9: '>',
                         10: '<',
                         11: ',',
                         12: '8',
                         13: 'h',
                         14: 'D',
                         15: 'H',
                         16: '|' },
           'line':     {  1: '-',
                          2: '--',
                          3: '-.',
                          4: ':' },
           'color':    {  1: 'blue',
                          2: '#0065BD',
                          3: '#989898',
                          4: 'red',
                          5: 'orange',
                          6: '#808080',
                          7: 'c',
                          8: 'black',
                          9: '#00FFFF',
                         10: '#A9A9A9',
                         11: '#BEBEBE',
                         12: '#D0D0D0',
                         13: '#DCDCDC',
                         14: '#E8E8E8',
                         15: '#F5F5F5',
                         16: '#F8F8F8' }
          }

if l_arguments['pdf'] != "":
  logging.info( "  file: "+l_arguments['pdf'] )
  l_pdfFile = PdfPages( l_arguments['pdf'] )

for l_enorm in ['l1', 'l2', 'linf']:
  logging.info( "plotting " + l_enorm + " norm" )

  l_figure = matplotlib.pyplot.figure( figsize=(10, 7) )

  l_leg = []

  # derive legend
  for l_config in range(len(l_yVals[l_enorm])):
    l_leg = l_leg + [ "O" +str(l_char[l_enorm][l_config]['order'])
                      + " Q" + str(l_char[l_enorm][l_config]['q'])
                      + " C" + str(l_char[l_enorm][l_config]['cfr']) ]

  # sort legend and plots
  l_plotRange = range(len(l_yVals[l_enorm]))
  l_plotRange = [l for (r,l) in sorted(zip(l_leg,l_plotRange))]
  l_leg.sort()

  for l_config in l_plotRange:
    matplotlib.pyplot.loglog( [ l_xVals[i] for i in range(len(l_yVals[l_enorm][l_config])) ],
                              l_yVals[l_enorm][l_config],
                              marker=l_style['markers'][l_char[l_enorm][l_config]['order']],
                              color=l_style['color'][l_char[l_enorm][l_config]['cfr']]  )

  # add infors
  matplotlib.pyplot.title( l_enorm )
  l_nCols = col=(len(l_yVals[l_enorm])/60+1)

  if l_arguments['legend']:
    matplotlib.pyplot.legend( l_leg, loc='lower right', fontsize=4, ncol=int(l_nCols)  )
  matplotlib.pyplot.xlabel('refinement')
  matplotlib.pyplot.ylabel('error')
  matplotlib.pyplot.grid( which='both' )
  if l_arguments['pdf'] == "":
    matplotlib.pyplot.show()
  else:
    l_pdfFile.savefig( l_figure )

  # close figure
  matplotlib.pyplot.close( l_figure )

# close file
if l_arguments['pdf'] != "":
  l_pdfFile.close()
