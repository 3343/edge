##
# @file This file is part of EDGE.
#
# @author Alexander Breuer (anbreuer AT ucsd.edu)
#         Alexander Heinecke (alexander.heinecke AT intel.com)
#
# @section LICENSE
# Copyright (c) 2015-2019, Regents of the University of California
# Copyright (c) 2016, Intel Corporation
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
# Build file.
##

import os
import subprocess
import xml.etree.ElementTree
import warnings
import SCons

##
# Gets the EDGE_version
##
def getEdgeVersion():
  # check if dirty
  l_dirty = subprocess.check_output(['git', 'diff', '--', '.', '\':(exclude)submodules\''])

  # assemble version
  l_version = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).strip()
  if( l_dirty != '' ):
    l_version = l_version + '_dirty'

  return l_version

##
# Adjust the given variable by turning relative paths to absolute paths
#
# i_var variable which is adjusted
##
def adjustPath( i_var ):
  l_var = i_var

  # only adjust if not boolean
  if( i_var != True and i_var != False ):
    # relative path is input
    if( i_var[0] != '/' ):
      l_var = os.path.join( Dir( '#'+i_var ).abspath )

  return l_var

def simpleWarning(message, category, filename, lineno, file=None, line=None):
    return '%s\n' % (message)
warnings.formatwarning = simpleWarning

# checkLibWithHeader with link flags
def CheckLibWithHeaderFlags( io_context, i_lib, i_header='', i_lang='CXX', i_flagsBefore=[''], i_flagsAfter=[''], i_dynamic=False ):

  l_lang, l_suffix, l_msg = SCons.Conftest._lang2suffix(i_lang)

  l_msg = 'Checking for '+l_lang
  if i_dynamic: l_msg=l_msg+' dynamic'
  else:         l_msg=l_msg+' static'
  l_msg = l_msg+' library '+i_lib+'..'
  io_context.Message( l_msg )

  # assemble source
  l_srcFile = "int main(int i_argc, char **i_argv) { return 0; }"
  if i_header != '':
    l_srcFile = "#include <"+ i_header+">\n"+l_srcFile

  # store old values
  if 'LIBS' in io_context.env:
    l_oldLibs = io_context.env['LIBS']
  else:
    l_oldLibs = []
  if 'LINKFLAGS' in io_context.env:
    l_oldBefore = io_context.env['LINKFLAGS']
  else:
    l_oldBefore = []
  if '_LIBFLAGS' in io_context.env:
    l_oldAfter = io_context.env['_LIBFLAGS']
  else:
    l_oldAfter = []

  # add link flags
  if i_dynamic:
    io_context.env.AppendUnique( _LIBFLAGS = ['-l'+i_lib] )
  else:
    io_context.env.PrependUnique( LIBS = [i_lib] )

  io_context.env.Prepend(  LINKFLAGS = i_flagsBefore     )
  io_context.env.Append(  _LIBFLAGS  = i_flagsAfter      )

  # test if it exists
  l_result = io_context.TryLink( l_srcFile, l_suffix )
  io_context.Result(l_result)

  # fall back to previous settings
  io_context.env['LIBS']      = l_oldLibs
  io_context.env['LINKFLAGS'] = l_oldBefore
  io_context.env['_LIBFLAGS'] = l_oldAfter

  # set library
  if l_result == 1:
    if( i_dynamic == False ):
      io_context.env.PrependUnique( LIBS = [i_lib] )
    else:
      # this is dirty: full support of static vs. dynamic linking in scons would be appreciated..
      io_context.env.AppendUnique( _LIBFLAGS = ['-l'+i_lib] )

  return l_result

# configuration
vars = Variables()

# add possible xml config
vars.AddVariables(
  PathVariable( 'xml',
                'xml configuration of the build, command line arguments have priority over xml.',
                 None,
                 PathVariable.PathIsFile )
)

# create environment
env = Environment( variables = vars )

# parse xml
xmlArgs = {}
if 'xml' in env:
  xmlTree = xml.etree.ElementTree.parse( env['xml'] )
  build = xmlTree.getroot().find('build')

  # get the xml args
  for option in build:
    xmlArgs[option.tag] = [option.text]

# add command line arguments
vars.AddVariables(
  EnumVariable( 'cfr',
                'concurrent forward runs',
                '16',
                 allowed_values=( '1', '2', '4', '8', '16' )
              ),
  EnumVariable( 'equations',
                'equations solved',
                'elastic',
                 allowed_values=( 'advection', 'elastic', 'swe' )
              ),
  EnumVariable( 'element_type',
                'element type used',
                'tet4',
                 allowed_values=('line', 'quad4r', 'tria3', 'hex8r', 'tet4' )
              ),
  EnumVariable( 'order',
                'order of convergence',
                '4',
                 allowed_values=('1', '2', '3', '4', '5', '6', '7', '8', '9')
              ),
  EnumVariable( 'mode',
                'compile modes, option \'san\' enables address and undefined behavior sanitizers',
                'release',
                 allowed_values=('release', 'debug', 'release+san', 'debug+san' )
              ),
  EnumVariable( 'arch',
                'architecture to compile for',
                'snb',
                 allowed_values=('snb', 'hsw', 'knl', 'skx', 'avx512')
              ),
  EnumVariable( 'precision',
                'floating point precision (bit)',
                '32',
                 allowed_values=('32', '64')
              ),
  EnumVariable( 'parallel',
                'used parallelization',
                'omp',
                 allowed_values=('none', 'omp', 'mpi', 'mpi+omp')
              ),
  BoolVariable( 'cov',
                'enable code coverage',
                 False ),
  BoolVariable( 'tests',
                'enable unit tests.',
                 False ),
  PathVariable( 'build_dir',
                'location where the code is build',
                'build',
                PathVariable.PathIsDirCreate ),
  PackageVariable( 'xsmm',
                   'enable libxsmm',
                   'yes' ),
  PackageVariable( 'zlib',
                   'enable zlib',
                   'no' ),
  PackageVariable( 'hdf5',
                   'enable HDF5',
                   'no' ),
  PackageVariable( 'moab',
                   'Enables the use of MOAB (The Mesh-Oriented datABase) and thus support for unstructured meshes if set. Otherwise regular meshes are used.',
                   'no' ),
  BoolVariable( 'easylogging',
                'Enables the use of Easylogging.',
                True ),
  BoolVariable( 'inst',
                'enable instrumentation',
                False )
)

# command args have priority
cArgs=vars.args
for cArg in cArgs:
  xmlArgs[cArg] = cArgs[cArg]

# forward combined args to scons
vars.args=xmlArgs

# include environment
env = Environment( variables = vars )

# exit in the case of unknown variables
if vars.UnknownVariables():
  print( "build configuration corrupted, don't know what to do with: " + str(vars.UnknownVariables().keys()) )
  exit(1)

# generate help message
Help( vars.GenerateHelpText(env) )

print( "##########################################################################" )
print( "##############   ##############            ###############  ##############" )
print( "##############   ###############         ################   ##############" )
print( "#####            #####       #####      ######                       #####" )
print( "#####            #####        #####    #####                         #####" )
print( "#############    #####         #####  #####                  #############" )
print( "#############    #####         #####  #####      #########   #############" )
print( "#####            #####         #####  #####      #########           #####" )
print( "#####            #####        #####    #####        ######           #####" )
print( "#####            #####       #####      #####       #####            #####" )
print( "###############  ###############         ###############   ###############" )
print( "###############  ##############           #############    ###############" )
print( "##########################################################################" )

# print welcome message
print( 'Running build script of EDGE.' )

# configuration
conf = Configure(env, custom_tests = {'CheckLibWithHeaderFlags': CheckLibWithHeaderFlags})

# include environment
env['ENV'] = os.environ

# enable default tool (otherwise SCons breaks on some systems); for some reason a repeat is required on some systems
env.Tool('default')
env.Tool('default')

# adjust path variables
for l_va in [ 'xsmm', 'zlib', 'hdf5', 'moab' ]:
  env[l_va] = adjustPath( env[l_va] )

# forward compiler
if 'CC' in env['ENV'].keys():
  env['CC'] = env['ENV']['CC']
if 'CXX' in env['ENV'].keys():
  env['CXX'] = env['ENV']['CXX']

# forward flags
if 'CFLAGS' in env['ENV'].keys():
  env['CFLAGS'] = env['ENV']['CFLAGS']
if 'CXXFLAGS' in env['ENV'].keys():
  env['CXXFLAGS'] = env['ENV']['CXXFLAGS']
if 'LIBS' in env['ENV'].keys():
  env['LIBS'] = env['ENV']['LIBS']
if 'LINKFLAGS' in env['ENV'].keys():
  env['LINKFLAGS'] = env['ENV']['LINKFLAGS']

# forward paths
if 'CPLUS_INCLUDE_PATH' in env['ENV'].keys():
  for l_incP in env['ENV']['CPLUS_INCLUDE_PATH'].split(':'):
    if l_incP != '':
      l_incP = adjustPath( l_incP )
      env.AppendUnique( CPPPATH = [l_incP] )
if 'LIBRARY_PATH' in env['ENV'].keys():
  for l_libP in env['ENV']['LIBRARY_PATH'].split(':'):
    if l_libP != '':
      l_libP = adjustPath( l_libP )
      env.AppendUnique( LIBPATH = [l_libP] )
      env.AppendUnique( RPATH   = [l_libP] )


# forward EDGE version
env.Append( CPPDEFINES='PP_EDGE_VERSION=\\"'+getEdgeVersion()+'\\"' )

# detect compilers
try:
  compVer = subprocess.check_output( [env['CXX'], "--version"] )
except:
  compVer = 'unknown'
if( 'g++' in compVer ):
  compilers='gnu'
elif( 'icpc' in compVer ):
  compilers='intel'
elif( 'clang' in compVer):
  compilers='clang'
else:
  compilers='gnu'
  print( '  no compiler detected' )
env['compilers']=compilers
print( '  using ' + compilers + ' as compiler suite' )

if env['PLATFORM'] == 'darwin':
  env['op_sys'] = 'macos'
else:
  env['op_sys'] = 'linux'
print( '  using ' + env['op_sys'] + ' as operating system' )

# use static linking for EDGE's direct dependencies (if possible) and dynamic for the rest
if env['op_sys'] != 'macos': # not supported by mac's ld
  env.PrependUnique( LINKFLAGS = ['-Wl,-Bstatic'] )
  env.AppendUnique( _LIBFLAGS = ['-Wl,-Bdynamic'] )

# disable libxsmm if not build elastic
if( env['xsmm'] ):
  if( 'elastic' not in env['equations'] ):
    warnings.warn('  Warning: LIBXSMM is not supported for equations other than elastic, continuing without' )
    env['xsmm'] = False
  if( env['order'] == '1' ):
    warnings.warn( '  Warning: LIBXSMM is not supported for first order runs, continuing without' )
    env['xsmm'] = False
  if( int(env['cfr']) > 1 and ( env['arch'] != 'snb' and env['arch'] != 'hsw' and env['arch'] != 'knl' and env['arch'] != 'skx' and env['arch'] != 'avx512')  ):
    warnings.warn( '  Warning: LIBXSMM not supported for fused simulations and arch != (snb, hsw, knl, skx or avx512), continuing without' )
    env['xsmm'] = False

  # disable libxsmm, if the number of fused simulations does not match the target-architecture
  if( env['cfr'] != '1' ):
    # Sandy Bridge
    if( env['arch'] == 'snb' ):
      if( env['precision'] == '32' and env['cfr'] != '8' ):
        warnings.warn( '  Warning: LIBXSMM disabled. Use 8 fused simulations for 32-bit precision and Sandy Bridge (snb)' )
        env['xsmm'] = False
      elif( env['precision'] == '64' and env['cfr'] != '4' ):
        warnings.warn( '  Warning: LIBXSMM disabled. Use 4 fused simulations for 64-bit precision and Sandy Bridge (snb)' )
        env['xsmm'] = False
    # Haswell
    elif( env['arch'] == 'hsw' ):
      if( env['precision'] == '32' and env['cfr'] != '8' ):
        warnings.warn( '  Warning: LIBXSMM disabled. Use 8 fused simulations for 32-bit precision and Haswell (hsw)' )
        env['xsmm'] = False
      elif( env['precision'] == '64' and env['cfr'] != '4' ):
        warnings.warn( '  Warning: LIBXSMM disabled. Use 4 fused simulations for 64-bit precision and Haswell (hsw)' )
        env['xsmm'] = False
    # Knights Landing
    elif( env['arch'] == 'knl' ):
      if( env['precision'] == '32' and env['cfr'] != '16' ):
        warnings.warn( '  Warning: LIBXSMM disabled. Use 16 fused simulations for 32-bit precision and Knights Landing (knl)' )
        env['xsmm'] = False
      elif( env['precision'] == '64' and env['cfr'] != '8' ):
        warnings.warn( '  Warning: LIBXSMM disabled. Use 8 fused simulations for 64-bit precision and Knights Landing (knl)' )
        env['xsmm'] = False
    # Skylake
    elif( env['arch'] == 'skx' ):
      if( env['precision'] == '32' and env['cfr'] != '16' ):
        warnings.warn( '  Warning: LIBXSMM disabled. Use 16 fused simulations for 32-bit precision and Skylake (skx)' )
        env['xsmm'] = False
      elif( env['precision'] == '64' and env['cfr'] != '8' ):
        warnings.warn( '  Warning: LIBXSMM disabled. Use 8 fused simulations for 64-bit precision and Skylake (skx)' )
        env['xsmm'] = False
    # AVX-512
    elif( env['arch'] == 'avx512' ):
      if( env['precision'] == '32' and env['cfr'] != '16' ):
        warnings.warn( '  Warning: LIBXSMM disabled. Use 16 fused simulations for 32-bit precision and AVX-512 (avx512)' )
        env['xsmm'] = False
      elif( env['precision'] == '64' and env['cfr'] != '8' ):
        warnings.warn( '  Warning: LIBXSMM disabled. Use 8 fused simulations for 64-bit precision and AVX-512 (avx512)' )
        env['xsmm'] = False

# forward number of forward runs to compiler
env.Append( CPPDEFINES='PP_N_CRUNS='+env['cfr'] )

# enable libhugetlbfs if available (dynamic because static misses functions, e.g., gethugepagesize())
if( conf.CheckLibWithHeaderFlags('hugetlbfs', '', 'CXX', [], [], True) ):
  env.AppendUnique( CPPDEFINES =['PP_USE_HUGETLBFS'] )

# enable libnuma if available
if conf.CheckLibWithHeaderFlags('numa', 'numa.h', 'CXX', [], [], True) or\
   conf.CheckLibWithHeaderFlags('numa', 'numa.h', 'CXX', [], [], False):
  env.AppendUnique( CPPDEFINES =['PP_USE_NUMA'] )

# set architecture and alignment
if env['arch'] == 'snb':
  env.Append( CPPFLAGS = ['-mavx'] )
elif env['arch'] == 'hsw':
  env.Append( CPPFLAGS = ['-march=core-avx2'] )
elif env['arch'] == 'avx512':
  if compilers=='gnu' or compilers=='clang':
    env.Append( CPPFLAGS = ['-mavx512f', '-mavx512cd'] )
  elif compilers=='intel':
    env.Append( CPPFLAGS = ['-xCOMMON-AVX512'] )
elif env['arch'] == 'skx':
  if compilers=='gnu' or compilers=='clang':
    env.Append( CPPFLAGS = ['-mavx512f', '-mavx512cd', '-mavx512bw', '-mavx512dq', '-mavx512vl'] )
  elif compilers=='intel':
    env.Append( CPPFLAGS = ['-xCORE-AVX512'] )
elif env['arch'] == 'knl' or env['arch'] == 'knm':
  if compilers=='gnu' or compilers=='clang':
    env.Append( CPPFLAGS = ['-mavx512f', '-mavx512cd', '-mavx512er', '-mavx512pf'] )
  elif compilers=='intel':
    env.Append( CPPFLAGS = ['-xMIC-AVX512'] )
  # add memkind if available (static is preferred)
  if conf.CheckLibWithHeaderFlags('memkind', 'hbwmalloc.h', 'CXX') or \
     conf.CheckLibWithHeaderFlags('memkind', 'hbwmalloc.h', 'CXX', [], [], True):
    env.AppendUnique( LINKFLAGS=['-pthread'] )
    env.Append( CPPDEFINES =['PP_USE_MEMKIND'] )

# forward equations to build
if env['equations'] == 'advection':
  env.Append( CPPDEFINES=['PP_T_EQUATIONS_ADVECTION']  )
elif 'elastic' in env['equations']:
  env.Append( CPPDEFINES=['PP_T_EQUATIONS_ELASTIC'] )
elif env['equations'] == 'swe':
  env.Append( CPPDEFINES=['PP_T_EQUATIONS_SWE'] )

# forward element type to build
if env['element_type'] == 'line':
  env.Append( CPPDEFINES=['PP_T_ELEMENTS_LINE'] )
elif env['element_type'] == 'quad4r':
  env.Append( CPPDEFINES=['PP_T_ELEMENTS_QUAD4R'] )
elif env['element_type'] == 'tria3':
  env.Append( CPPDEFINES=['PP_T_ELEMENTS_TRIA3'] )
elif env['element_type'] == 'hex8r':
  env.Append( CPPDEFINES=['PP_T_ELEMENTS_HEX8R'] )
elif env['element_type'] == 'tet4':
  env.Append( CPPDEFINES=['PP_T_ELEMENTS_TET4'] )
else:
  assert( False );

# forward order
env.Append( CPPDEFINES=['PP_ORDER='+env['order']] )

# forward precision
env.Append( CPPDEFINES=['PP_PRECISION='+env['precision']] )

# enable omp
if 'omp' in env['parallel']:
  env.Append( CPPDEFINES = ['PP_USE_OMP'] )

  if compilers == 'intel':
    env.Append( CPPFLAGS = ['-qopenmp'] )
    env.Append( LINKFLAGS = ['-qopenmp'] )
    env.Append( LINKFLAGS = ['-qopenmp-link=static'] )
  elif compilers == 'gnu' or compilers == 'clang':
    env.Append( CPPFLAGS = ['-fopenmp'] )
    env.Append( LINKFLAGS = ['-fopenmp'] )

# fix clang's confusion with old libstdc++ at runtime (no proper automatic rpath)
if compilers == 'clang':
  try:
    rpath = subprocess.check_output( [env['CXX'], "--print-file-name=libstdc++.so"] )
    rpath = os.path.dirname( rpath )
    env.Append( LINKFLAGS = [ '-Wl,-rpath='+rpath ] )
  except:
    print( '  failed getting dir of libstdc++.so, not setting rpath' )

# enable mpi
if 'mpi' in env['parallel']:
  env.Append( CPPDEFINES = ['PP_USE_MPI'] )

# add current path to seach path
env.Append( CPPPATH = ['#', '#/src'] )

# add default flags
env.Append( CXXFLAGS = ["-std=c++11", "-Wall", "-Wextra", "-Wno-unknown-pragmas", "-Wno-unused-parameter", "-Werror"] )

if env['inst'] == False:
  env.Append( CXXFLAGS = ["-pedantic", "-Wshadow"] ) # some strict flags break compilation with opari..
if compilers != 'intel':
  env.Append( CXXFLAGS = ["-Wundef"] ) # intel compiler gets this flag back if we can define system headers as in GCC..
else:
  # silence annoying warnings
  env.Append( CXXFLAGS = ["-diag-disable", "186,11074,11076"] )

# set optimization mode
if 'debug' in env['mode']:
  env.Append( CXXFLAGS = ['-g','-O0'] )
else:
  env.Append( CXXFLAGS = ['-O2'] )
  if compilers=='gnu':
    env.Append( CXXFLAGS = '-ftree-vectorize' )
# add sanitizers
if 'san' in  env['mode']:
  env.Append( CXXFLAGS =  ['-g', '-fsanitize=address', '-fsanitize=undefined', '-fno-omit-frame-pointer'] )
  env.Append( LINKFLAGS = ['-g', '-fsanitize=address', '-fsanitize=undefined'] )

# enable code coverage, if requested
if env['cov'] == True:
  env.Append( CXXFLAGS = ['-coverage', '-fno-inline', '-fno-inline-small-functions', '-fno-default-inline'] )
  env.Append( LINKFLAGS = ['-coverage'] )

# add Catch
if env['tests']:
  env.Append( CXXFLAGS = ['-Isubmodules/Catch/include'] )

# get source files
VariantDir( env['build_dir']+'/src', 'src')
VariantDir( env['build_dir']+'/submodules', 'submodules')

env.sources = []
env.tests = []

Export('env')
Export('conf')
SConscript( env['build_dir']+'/submodules/SConscript' )
Import('conf')
Import('env')

Export('env')
SConscript( env['build_dir']+'/src/SConscript' )
Import('env')

# prepend scorep in case of instrumentation and add flag
# WARNING: leave scorep at the very bottom (outherwise it will be used in CheckLib-tests
if env['inst']:
  # disable most of scorep
  scorep = "scorep --thread=omp --nocompiler --user "

  env['CC']  = scorep+env['CC']
  env['CXX'] = scorep+env['CXX']
  env.Append( CPPDEFINES = ['PP_USE_INSTR'] )

# add a new line
print( '' )
env.Program( env['build_dir']+'/edge', source = env.sources )

if env['tests']:
  env.Program( env['build_dir']+'/tests', source = env.tests )
