##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2017, Regents of the University of California
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
# Generates sampled slip-rates for kinematic source benchmarks.
##
import logging
import argparse
import math
import netCDF4
import subprocess
import tempfile
import xml.etree.ElementTree

##
# Parses an XML-config containing a finite fault description.
#
# @param i_inFile XML-file which is parsed
# @return 1) Parsed config (dict) 2) number of assumed dimensions.
##
def parseConfig( i_inFile ):
  l_tree = xml.etree.ElementTree.parse( i_inFile )
  l_root = l_tree.getroot()

  l_dict = { 'subfaults': [] }

  # determine number of dimensions
  l_nDims = 3
  l_dimsN = ['x', 'y', 'z']
  l_dirsN = ['n', 's', 't']
  for l_subfault in l_root.iter('subfault'):
    if( l_subfault.find('center').find( 'z') == None or
        l_subfault.find(  'dirs').find('sz') == None or
        l_subfault.find(  'dirs').find('tx') == None or
        l_subfault.find(  'dirs').find('ty') == None or
        l_subfault.find(  'dirs').find('tz') == None ):
      l_nDims = 2
      l_dimsN = ['x', 'y']
      l_dirsN = ['n', 's']
      logging.info( "couldn't find all 3D attributes, assuming a 2D kinematic source." )
    else:
      logging.info( "found all 3D attributes, assuming a 3D kinematic source." )

  # get common options
  l_dict['cdl'] = l_root.find('path_to_cdl').text
  l_dict['out']  = l_root.find('path_to_out').text

  # iterate over subfaults
  for l_subfault in l_root.iter('subfault'):
    l_dict['subfaults'] = l_dict['subfaults'] + [{}]

    for l_par in ['on', 'dt', 'dur', 'mu', 'a', 'm0', 'T']:
      l_find = l_subfault.find(l_par)
      if( l_find != None ):
        l_dict['subfaults'][-1][l_par] = float( l_find.text )

    # read type of moment-rate time history
    if( l_subfault.find('mrth') != None ):
      l_mrth = l_subfault.find('mrth').text
    else:
      if( l_nDims == 2 ): l_mrth = 'ricker'
      else:               l_mrth = 'wp'
    # store
    l_dict['subfaults'][-1]['mrth'] = l_mrth

    l_dict['subfaults'][-1]['center'] = []
    l_center = l_subfault.find('center')
    for l_di in l_dimsN:
      l_dict['subfaults'][-1]['center'] = l_dict['subfaults'][-1]['center'] +\
                                          [ float( l_center.find(l_di).text ) ]

    l_dict['subfaults'][-1]['dirs'] = {}
    l_dirs = l_subfault.find('dirs')
    for l_dir in l_dirsN:
      l_dict['subfaults'][-1]['dirs'][l_dir] = []
      for l_di in l_dimsN:
        l_str = l_dir+l_di
        l_dict['subfaults'][-1]['dirs'][l_dir] = l_dict['subfaults'][-1]['dirs'][l_dir] +\
                                                 [ float( l_dirs.find(l_str).text ) ]

  return l_nDims, l_dict

##
# Creates a netCDF-file based on the given cdl-source template.
#
# @param i_inFile cdl-template, with the two template paramters TEMPLATE_CENTERS and TEMPLATE_SUBFAULTS.
# @param i_outFile netCDF-output file which will be created.
# @param i_nSrcStr string replacing TEMPLATE_N_SRC;
# @param i_centerStr string replacing TEMPLATE_CENTERS.
# @param i_subStr string replacing TEMPLATE_SUBFAULTS.
##
def createNcFile( i_inFile, i_outFile, i_nSrcStr, i_centerStr, i_subStr ):
  # operate on a temporary file
  with tempfile.NamedTemporaryFile() as l_tmpFile:
    # open template
    with open(i_inFile) as l_tmpl:
      # copy contents and replace templates
      for l_line in l_tmpl:
        l_outStr = l_line

        l_outStr = l_outStr.replace( 'TEMPLATE_N_SRC',
                                     i_nSrcStr )

        l_outStr = l_outStr.replace( 'TEMPLATE_CENTERS',
                                     i_centerStr )
        l_outStr = l_outStr.replace( 'TEMPLATE_SUBFAULTS',
                                     i_subStr )
        l_tmpFile.write( str.encode( l_outStr ) )

      # jump to beginning
      l_tmpFile.seek(0)

      # translate to netCDF
      l_proc = subprocess.check_call( ['ncgen -o '+i_outFile+' '+l_tmpFile.name], shell=True )

##
# Writes the given slip-rate samples of a single point source the corresponding dimensions and metainfo.
# In the case of three-dimensional source, the slip-rate samples are written to var sliprates1.
# In the case of a two-dimnesional source, the slip-rate samples are written to vars sliprates1 and sliprates2 (explosive source).
#
# @param i_nDims number of slip-directions.
# @param i_outFile netCDF-file where the slip-rate samples are written to.
# @param i_srs slip-rate samples which are written. list of lists [*][]: subfault, [][*]: sample
##
def writeSrs( i_nDims, i_outFile, i_srs ):
  # warn user about upcoming warning
  logging.info( 'warning \'unsupported Compound type, skipping...\' can be ignored, we are not touching this part of the netCDF-file' )

  # open netCDF file
  l_rootGroup = netCDF4.Dataset( i_outFile, "r+", format="NETCDF4" )

  # create offset dimension
  l_rootGroup.createDimension( 'sroffset',  len(i_srs)+1 )
  l_rootGroup.createDimension( 'direction', i_nDims )

  # create offset variable
  l_rootGroup.createVariable( varname='sroffsets',
                              datatype='u4',
                              dimensions=('sroffset', 'direction') )

  # determine total size
  l_nSrs = 0
  for l_sr in i_srs: l_nSrs = l_nSrs + len(l_sr)

  # create slip-rate dimensions
  l_rootGroup.createDimension( 'sample1', l_nSrs )
  l_rootGroup.createDimension( 'sample2', None       )
  if( i_nDims > 2 ):
    l_rootGroup.createDimension( 'sample3', None       )

  # create slip-rate variables
  l_rootGroup.createVariable( varname='sliprates1',
                              datatype='f8',
                              dimensions=('sample1') )
  l_rootGroup.createVariable( varname='sliprates2',
                              datatype='f8',
                              dimensions=('sample2') )
  if( i_nDims > 2 ):
    l_rootGroup.createVariable( varname='sliprates3',
                                datatype='f8',
                                dimensions=('sample3') )

  # set units
  l_rootGroup['sliprates1'].units = "m/s"
  l_rootGroup['sliprates2'].units = "m/s"
  if( i_nDims > 2 ):
    l_rootGroup['sliprates3'].units = "m/s"

  # write slip-rates
  l_first = 0
  l_flush = 0;
  for l_su in i_srs:
    l_rootGroup['sliprates1'][l_first:l_first+len(l_su)] = l_su
    if( i_nDims > 2 ): l_rootGroup['sliprates2'][l_first:l_first+len(l_su)] = l_su
    l_first = l_first + len(l_su)

    # flush periodically
    l_flush = l_flush + 1
    if l_flush % 10 == 0:
      l_rootGroup.sync

  # write offsets
  l_rootGroup['sroffsets'][0, :] = ( [ 0, 0, 0 ] if( i_nDims > 2 ) else [0, 0] )
  l_off = 0
  for l_so in range( 0, len(i_srs) ):
    l_off = l_off + len(i_srs[l_so])
    l_rootGroup['sroffsets'][l_so+1, :] = ( [ l_off, 0, 0 ] if( i_nDims > 2 ) else [ l_off, l_off ] )

  l_rootGroup.close()

##
# Generates slip rate samples as given by the Ricker wavelet.
#
# Source: Geophysical Journal International 166.2 (2006): 855-877.
#         An arbitrary high-order discontinuous Galerkin method for
#         elastic waves on unstructured meshes I. The two-dimensional isotropic case with external source terms
#         Eq. (64)
#
# Moment rate time history:  a*(0.5+b(t-T)^2)*exp(b*(t-T)^2) )
#
# @param i_t1 start time of the sampling.
# @param i_t2 end time of the sampling.
# @param i_dt time step of the sampling.
# @param i_a width parameter of the wavelet.
# @param i_M0 scaling used for the moment-rate time history
# @return sampled slip rates.
##
def ricker( i_t1,
            i_t2,
            i_dt,
            i_f=1.0,
            i_M0=1.0):
  # set parameters as used in the equation above
  l_a = i_M0
  l_b = -(math.pi * i_f)**2

  # slip rate samples
  l_sr = []

  # local time
  l_t = i_t1

  while( l_t < i_t2 ):
    l_sr.append( l_a * ( 0.5 + l_b * l_t**2 ) * math.exp( l_b * l_t**2 ) )
    l_t = l_t + i_dt

  return l_sr

##
# Generates slip-rate samples for the SISMOWINE wave propagation benchmarks (see http://www.sismowine.org/model.html )
#
# @param i_t1 start time of the sampling.
# @param i_t2 end time of the sampling.
# @param i_dt time step of the sampling.
# @param i_T parameter T in the analytical slip-rate function.
# @param i_M0 scaling used for the moment-rate time history
# @return samples slip rates.
##
def wp( i_t1,
        i_t2,
        i_dt,
        i_T=0.1,
        i_M0=-1E18 ):
  l_sr = []

  # local time
  l_t = i_t1

  l_tI  = 1.0 / i_T
  l_tsI = l_tI * l_tI

  while( l_t < i_t2 ):
    l_sr.append( l_t * l_tsI * math.exp( -l_t * l_tI ) * i_M0 )
    l_t = l_t+i_dt
  return l_sr

##
# Generates slip-rate samples for boxcar.
#
# @param i_t1 start time of the sampling.
# @param i_t2 end time of the sampling (i_t2-i_t1 is the rise-time).
# @param i_dt time step of the sampling.
# @param i_m0 scaling used for the moment-rate time history
# @return samples slip rates.
##
def boxcar( i_t1,
            i_t2,
            i_dt,
            i_m0 ):
  l_sr = []

  # local time
  l_t = i_t1
  # intergral of boxcar is 1 + scale by moment
  l_m0DivDur = i_m0 / (i_t2-i_t1)

  while( l_t < i_t2 ):
    l_sr.append( l_m0DivDur )
    l_t = l_t+i_dt
  return l_sr

##
# Generates slip-rate samples for the Can1 benchmark of the SISMOWINE benchmark suite.
#
# Moment rate time history (WolframAlpha):
# in: derivative of exp( - ( w * ( t - T ) /y )^2 ) * cos( w * ( t - T) + v )
# out: -((w (2 (t - T) w Cos[v + (t - T) w] + y^2 Sin[v + (t - T) w]))/(E^(((t - T)^2 w^2)/y^2) y^2))
#
# @param i_t1 start time of the sampling.
# @param i_t2 end time of the sampling.
# @param i_dt time step of the sampling.
# @param i_M0 scaling used for the moment-rate time history.
# @param i_w parameter omega: 2 * PI * f_Pin the analytical slip-rate function.
# @param i_y parameter gamma in the analytical slip-rate function.
# @param i_v parameter v in the analytical slip-rate function.
# @param i_ts parameter t_S in the analytical slip-rate function.
# @return samples slip rates.
##
def can1( i_t1,
          i_t2,
          i_dt,
          i_M0=-1E18,
          i_w = 2.0 * math.pi * 1.5,
          i_y = 2.0,
          i_v = math.pi * 0.5,
          i_ts = 1.0 ):
  # slip rate samples
  l_sr = []

  # local time
  l_t = i_t1

  # iterate over samples
  while( l_t < i_t2 ):
    # compute sample
    l_sa = i_y**2 * math.sin( i_w * ( l_t - i_ts ) + i_v )
    l_sa = l_sa +   2.0 * i_w * ( l_t - i_ts ) \
                  * math.cos( i_w * ( l_t - i_ts ) + i_v )
    l_sa = l_sa * -i_w * math.exp( - ( i_w**2 * ( l_t - i_ts )**2 ) / i_y**2 )
    l_sa = l_sa / i_y**2

    # scale
    l_sa = l_sa * i_M0

    # attach sample
    l_sr.append( l_sa )
    l_t = l_t + i_dt

  # return slip rate samples
  return l_sr

# set up logger
logging.basicConfig( level=logging.DEBUG,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# command line arguments
l_parser = argparse.ArgumentParser( description='Generates a NRF file from a finite fault config and a CDL-template. Writes sampled slip-rates for kinematic source benchmarks.',
                                    formatter_class=argparse.RawTextHelpFormatter, )

l_parser.add_argument( '--xml',
                       dest     = 'xml',
                       required = True,
                       type     = str,
                       help     =
'\
XML-configuration for the kinematic source description.\n\
The expected syntax is as follows:\n\
  <?xml version="1.0"?>\n\
  <!-- finite fault description consisting out of one or more subfaults -->\n\
  <finite_fault>\n\
    <!-- subfault struct, arbitrary number of subfaults are allowed -->\n\
    <subfault>\n\
      <!-- onset time, in 2D setups this is used as time shift for the Ricker wavelet -->\n\
      <on>0.0</on>\n\
      <!-- time step of the slip-rate sampling -->\n\
      <dt>0.001</dt>\n\
      <!-- duration for which the subfault is active -->\n\
      <dur>4.0</dur>\n\
      <!-- Lame parameter mu. If set to 0, the background velocity model is used -->\n\
      <mu>1.0</mu>\n\
      <!-- area of the subfault -->\n\
      <a>1.0</a>\n\
\n\
      <!-- scaling, used for the computation of the moment-rate time history,\n\
           typically called \'M0\' -->\n\
      <m0>-1E18</m0>\n\
      <!-- smoothing time T, used for the generation of the moment-rate time series -->\n\
      <T>0.1</T>\n\
\n\
      <!-- moment-rate time history, 2D: ricker (default); 3D: wp (default) or boxcar or can1 -->\n\
      <mrth>wp</mrth>\n\
\n\
      <!-- center of the subfault (x, y, z coordinates) -->\n\
      <center>\n\
        <x>   0</x>\n\
        <y>   0</y>\n\
        <z>2000</z>\n\
      </center>\n\
\n\
      <dirs>\n\
        <!-- direction of fault-normal (x, y, z) -->\n\
        <nx>0</nx>\n\
        <ny>1</ny>\n\
        <nz>0</nz>\n\
\n\
        <!-- direction of first tangent (x, y, z) -->\n\
        <sx>1</sx>\n\
        <sy>0</sy>\n\
        <sz>0</sz>\n\
\n\
        <!-- direction of second tangent (x, y, z) -->\n\
        <tx>0</tx>\n\
        <ty>0</ty>\n\
        <tz>1</tz>\n\
      </dirs>\n\
    </subfault>\n\
\n\
    <!-- path to cdl-file, which is used to generate the basic netCDF-header\n\
         and which sliprates1-array is populated -->\n\
    <path_to_cdl>kinematic_bench.cdl</path_to_cdl>\n\
    <!-- path to output file (netCDF) to which the slip-rates the output is written -->\n\
    <path_to_out>example_src.nc</path_to_out>\n\
  </finite_fault>\n\
')

logging.info( "running kinematic benchmark script.." )

l_args = vars(l_parser.parse_args())

l_nDims, l_cfg = parseConfig( l_args['xml'] )
assert( len( l_cfg['subfaults']) > 0 )

# assemble center str
l_ctrStr = ''
l_first = True;
for l_sf in l_cfg['subfaults']:
  if not l_first: l_ctrStr = l_ctrStr + ', '
  else: l_first = False

  l_ctrStr = l_ctrStr + '{'
  l_ctrStr = l_ctrStr +         str( l_sf['center'][0] )                       +\
                         ', ' + str( l_sf['center'][1] )                       +\
                        (', ' + str( l_sf['center'][2] ) if l_nDims>2 else '')
  l_ctrStr = l_ctrStr + '}'

# assemble subfault string
l_subStr = ''
l_first = True;
for l_sf in l_cfg['subfaults']:
  if not l_first: l_subStr = l_subStr + ', '
  else: l_first = False

  l_subStr = l_subStr + '{'
                          # onset time shifts ricker wavelet in 2D setups
  l_subStr = l_subStr +   str( l_sf['on'] if l_nDims>2 else 0 ) + ', ' +\
                          str( l_sf['dt'] ) + ', ' +\
                          str( l_sf['mu'] ) + ', ' +\
                          str( l_sf[ 'a'] ) + ', ' +\
                          '{' +\
                                    str( l_sf['dirs'][ 's'][0] )                       +\
                             ', ' + str( l_sf['dirs'][ 's'][1] )                       +\
                            (', ' + str( l_sf['dirs'][ 's'][2] ) if l_nDims>2 else '') +\
                          '}' + ', ' +\
                          ('{' if l_nDims>2 else '') +\
                            (       str( l_sf['dirs'][ 't'][0] ) if l_nDims>2 else '') +\
                            (', ' + str( l_sf['dirs'][ 't'][1] ) if l_nDims>2 else '') +\
                            (', ' + str( l_sf['dirs'][ 't'][2] ) if l_nDims>2 else '') +\
                          ('},' if l_nDims>2 else '') +\
                          '{' +\
                             str( l_sf['dirs'][ 'n'][0] )                              +\
                             ', ' + str( l_sf['dirs'][ 'n'][1] )                       +\
                            (', ' + str( l_sf['dirs'][ 'n'][2] ) if l_nDims>2 else '') +\
                          '}' +\
                        '}'

createNcFile( l_cfg['cdl'],
              l_cfg['out'],
              str( len(l_cfg['subfaults']) ),
              l_ctrStr,
              l_subStr )

logging.info( "created netCDF-file "+l_cfg['out'] );

# slip rate samples
l_srs = []

# create slip rates
l_nSrs = 0
for l_sf in l_cfg['subfaults']:
  if( l_sf['mrth'] == 'wp' ):
    l_sr = wp( 0.0,
               l_sf['dur'],
               l_sf['dt'],
               l_sf['T'],
               l_sf['m0'] )

  if( l_sf['mrth'] == 'boxcar' ):
    l_sr = boxcar( l_sf['on'],
                   l_sf['dur']+l_sf['on'],
                   l_sf['dt'],
                   l_sf['m0'] )

  elif( l_sf['mrth'] == 'can1'):
    l_sr = can1( l_sf['on'],
                 l_sf['dur']-l_sf['on'],
                 l_sf['dt'],
                 l_sf['m0'] )

  elif( l_sf['mrth'] == 'ricker' ):
    l_sr = ricker( l_sf['on'],
                   l_sf['dur']+l_sf['on'],
                   l_sf['dt'],
                   l_sf['T'],
                   l_sf['m0'] )
  l_srs  = l_srs + [l_sr]
  l_nSrs = l_nSrs + len( l_sr )

logging.info( "generated a total of "+str(l_nSrs)+" slip-rate samples for " + str(len(l_cfg['subfaults'])) + ' subfaults' )
writeSrs( l_nDims, l_cfg['out'], l_srs )

logging.info( "finished, enjoy your benchmark source" )
