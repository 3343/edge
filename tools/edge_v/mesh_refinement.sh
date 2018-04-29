#!/bin/bash
##
# @file This file is part of EDGE.
#
# @author Rajdeep Konwar (rkonwar AT ucsd.edu)
#
# @section LICENSE
# Copyright (c) 2018, Regents of the University of California
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
# Script to generate a velocity-based mesh refinement given its geo file
# Tested on: Ubuntu 18.04 LTS
##

# Usage info
function show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-m MODEL -c CONFDIR -o MSHDIR] [-p MSHFILES] [-n ITERATIONS] [-r REMOTEMSH -u USR -d DOMAIN -g REMOTEGMSH -t REMOTEMSHDIR]
Generates velocity-model (UCVM) based refined mesh iteratively.
    -h This help message.
    -m MODEL model name (name of geo file).
    -c CONFDIR config directory.
    -o MSHDIR mesh directory.
    -p MSHFILES handling intermediate mesh files, i.e. msh, pos & logs (optional)
      1: Generate and zip intermediate files (by default).
      2: Generate but don't zip intermediate files.
      3: Do not generate intermediate files (i.e. only generate refined mesh).
    -n ITERATIONS number of iterations (optional, by default 10).
    -r REMOTEMSH remote meshing option (optional, by default 0).
    -u USR remote username (optional).
    -d DOMAIN remote domain (optional).
    -g REMOTEGMSH remote location of Gmsh executable (optional).
    -t REMOTEMSHDIR remote mesh directory (optional).
EOF
}

# Function to mesh on remote client
function mesh() {
  declare bgmPos=''

  if [ $REMOTEMSH -eq 1 ]
  then
    # Remote location of geo file
    printf -v l_geoFile '%s%s.geo'           "$REMOTEMSHDIR" "$MODEL"

    # Remote location of mesh file
    if [ $iter -eq 0 ]
    then
      printf -v l_mshFile '%s%s.msh'         "$REMOTEMSHDIR" "$MODEL"
      printf -v l_posFile '%s%s.pos'         "$REMOTEMSHDIR" "$MODEL"
    elif [ $iter -eq $ITERATIONS ]
    then
      printf -v l_mshFile '%s%s_refined.msh' "$REMOTEMSHDIR" "$MODEL"
      printf -v l_posFile '%s%s_refined.pos' "$REMOTEMSHDIR" "$MODEL"
    else
      printf -v l_mshFile '%s%s_%d.msh'      "$REMOTEMSHDIR" "$MODEL" "$iter"
      printf -v l_posFile '%s%s_%d.pos'      "$REMOTEMSHDIR" "$MODEL" "$iter"
    fi

    # Send pos file from local to remote
    if [ $iter -ne 0 ]
    then
      printf "\nSending $posFile from local m/c to remote client...\n"

      SECONDS=0
      scp -v $posFile $l_remoteLoc
      duration=$SECONDS

      printf "Time taken: $duration s\n\n"
      printf -v bgmPos '%s %s' "-bgm" "$l_posFile"
    fi

    # Remote machine mesh
    printf -v l_login '%s@%s' "$USR" "$DOMAIN"
    ssh $l_login nohup $REMOTEGMSH $l_geoFile $bgmPos -3 -o $l_mshFile -optimize_threshold $1 2>&1 | tee $mshLog

    # Copy mesh file back to local machine
    printf "\nSending $l_mshFile from remote client to local m/c...\n"
    printf -v l_retMshFile '%s:%s' "$l_login" "$l_mshFile"

    SECONDS=0
    scp -v $l_retMshFile $MSHDIR
    duration=$SECONDS
    printf "Time taken: $duration s\n\n"
  else
    # Mesh on local machine
    if [ $iter -eq 0 ]
    then
      gmsh $geoFile -3 -o $mshFile -optimize_threshold $1 2>&1 | tee $mshLog
    else
      gmsh $geoFile -bgm $posFile -3 -o $mshFile -optimize_threshold $1 2>&1 | tee $mshLog
    fi
  fi
}


#######################
##### Entry Point #####
#######################

# Global variables
declare -g iter=0           # Iteration counter
declare -g ITERATIONS=10    # Number of iterations
declare -g geoFile          # geo file
declare -g mshFile          # msh file
declare -g mshLog           # msh log file
declare -g posFile          # pos file
declare -g MSHFILES=1       # Handling intermediate files, see show_help() above for more info
declare -g shPath=$(pwd)    # Current directory (from where this script was run)

# ssh variable (use REMOTEMSH=1 to mesh on a remote client; by default 0)
# IMPORTANT: use only if you have a working SSH public key with the remote client!!
declare -g REMOTEMSH=0

# Parse command line arguments
OPTIND=1

while getopts "hm:c:o:p:n:r:u:d:g:t:" opt; do
  case "$opt" in
    h)
      show_help                         # Show help
      exit 0
      ;;
    m)
      declare -g MODEL=$OPTARG          # Model name
      ;;
    c)
      CONFDIR=$OPTARG                   # Config directory (local)
      ;;
    o)
      MSHDIR=$OPTARG                    # msh directory (local)
      ;;
    p)
      MSHFILES=$OPTARG                  # Handling intermediate files
      ;;
    n)
      ITERATIONS=$OPTARG                # Number on iterations
      ;;
    r)
      REMOTEMSH=$OPTARG                 # Remote mesh functionality
      ;;
    u)
      declare -g USR=$OPTARG            # ssh login user name
      ;;
    d)
      declare -g DOMAIN=$OPTARG         # Remote client domain name
      ;;
    g)
      declare -g REMOTEGMSH=$OPTARG     # gmsh path on remote client
      ;;
    t)
      declare -g REMOTEMSHDIR=$OPTARG   # msh directory on remote client
      ;;
    \?)
      show_help >&2
      exit 1
      ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

# -m, -c, -o are minimum options required to run script
# If -r 1 is used, we need rest of 4 parameters i.e. -u, -d, -g, -t
if [[ ( $OPTIND -lt 7 ) || ( $REMOTEMSH -ne 0  && ( -z "$USR"         || \
                                                    -z "$DOMAIN"      || \
                                                    -z "$REMOTEGMSH"  || \
                                                    -z "$REMOTEMSHDIR" ) ) ]]
then
  show_help >&2
  exit 1
fi

# Export write_pos directory (./bin) to PATH if not set
if ! command -v "write_pos" > /dev/null
then
  printf "Trying to export write_pos to PATH automatically... "
  if [ -d "bin" ] && [ -f "bin/write_pos" ]
  then
    printf -v posPath "export PATH=$shPath/bin:$PATH"
    $posPath
    printf "Done!\n"
  else
    printf "unable to locate!\nPlease export write_pos directory (typically ./bin) to PATH and retry.\n"
    exit 1
  fi
fi

# print info on the run
echo "$(date) Running mesh refinement..."

# Affected files
printf -v confFile '%s%s.conf' "$CONFDIR" "$MODEL"
printf -v geoFile  '%s%s.geo'  "$MSHDIR"  "$MODEL"
printf -v mshFile  '%s%s.msh'  "$MSHDIR"  "$MODEL"

# Send the geo file to remote server
if [ $REMOTEMSH -eq 1 ]
then
  # Create remote msh dir (if doesn't exist)
  printf -v l_login '%s@%s' "$USR" "$DOMAIN"
  ssh $l_login mkdir -p $REMOTEMSHDIR

  # Copy geo file to remote client
  printf "\nSending $geoFile from local m/c to remote client...\n"
  printf -v l_remoteLoc '%s:%s' "$l_login" "$REMOTEMSHDIR"

  SECONDS=0
  scp -v $geoFile $l_remoteLoc
  duration=$SECONDS
  printf "Time taken: $duration s\n\n"
fi

# Intial coarse mesh
printf "Generating initial (coarse) mesh...\n"
printf -v mshLog '%s.log' "$mshFile"

# Optimize tetrahedral elements that have a quality less than a threshold (0.3)
mesh 0.3

for ((iter=1; iter<$ITERATIONS; iter++))
do
  # Create a config file for each iteration
  printf -v confFileIt '%s%s_%d.conf' "$CONFDIR" "$MODEL" "$iter"
  cp $confFile $confFileIt

  printf "\n###########\nIteration $iter\n###########\n\n"

  # Pos file name for current iteration
  printf -v posFile '%s%s_%d.pos' "$MSHDIR"  "$MODEL" "$iter"

  # Replace input mesh_file and output pos_file in config file
  printf -v msh '/mesh_file=/c\mesh_file=%s' "$mshFile"
  sed -i -e "$msh" $confFileIt
  
  printf -v pos '/pos_file=/c\pos_file=%s'   "$posFile"
  sed -i -e "$pos" $confFileIt

  # Mesh file (to be generated) which will be used in next iteration
  printf -v mshFile '%s%s_%d.msh' "$MSHDIR"  "$MODEL" "$iter"
  printf -v mshLog  '%s.log'      "$mshFile"

  # Generate intermediate pos file to be used by gmsh in next step to produce a 
  # more refined mesh than previous iteration (using pos as background mesh)
  printf -v posLog '%s.log' "$posFile"
  printf "Generating pos...\n"
  write_pos -f $confFileIt 2>&1 | tee $posLog

  # Generate intermediate meshes
  printf "\nGenerating intermediate ($iter) mesh...\n"
  mesh 0.3

  # Handling intermediate files
  if [ $MSHFILES -eq 1 ]
  then
    printf -v mshFileZ '%s_%d.msh'    "$MODEL" "$iter"
    printf -v mshLogZ  '%s.log'       "$mshFileZ"
    printf -v posFileZ '%s_%d.pos'    "$MODEL" "$iter"
    printf -v posLogZ  '%s.log'       "$posFileZ"
    printf -v zipFile  '%s_%d.tar.gz' "$MODEL" "$iter"
    printf "\nCompressing intermediate files to %s...\n" "$zipFile"
    cd $MSHDIR
    tar -zcvf $zipFile $mshFileZ $mshLogZ $posFileZ $posLogZ
    cd $shPath
  fi

  if [ $MSHFILES -eq 1 ] || [ $MSHFILES -eq 3 ]
  then
    # Remove previous iteration intermediate files
    if [ $iter -ne 1 ]
    then
      printf "\nRemoving iteration %d intermediate files..." $(($iter-1))
      printf -v mshFileRm '%s%s_%d.msh' "$MSHDIR"  "$MODEL"  $(($iter-1))
      printf -v mshLogRm  '%s.log'      "$mshFileRm"
      printf -v posFileRm '%s%s_%d.pos' "$MSHDIR"  "$MODEL"  $(($iter-1))
      printf -v posLogRm  '%s.log'      "$posFileRm"
      rm -f $mshFileRm $mshLogRm $posFileRm $posLogRm
      printf " Done!\n"
    fi
  fi
done

# Create a config file for the final refined mesh
printf -v confFileRe '%s%s_refined.conf' "$CONFDIR" "$MODEL"
cp $confFile $confFileRe

# Modifying file paths in config file
printf -v msh '/mesh_file=/c\mesh_file=%s' "$mshFile"
sed -i -e "$msh" $confFileRe

printf -v posFile '%s%s_refined.pos' "$MSHDIR" "$MODEL"
printf -v pos '/pos_file=/c\pos_file=%s' "$posFile"
sed -i -e "$pos" $confFileRe

printf -v mshFile '%s%s_refined.msh' "$MSHDIR" "$MODEL"
printf -v mshLog '%s_refined.log' "$mshFile"

# Final refined pos file
printf -v posLog '%s_refined.log' "$posFile"
printf "\nGenerating pos...\n"
write_pos -f $confFileRe 2>&1 | tee $posLog

# Final mesh-refinement
printf '\nGenerating final (refined) mesh...\n'
mesh 0.3

# Remove last iteration's intermediate files
if [ $ITERATIONS -gt 1 ]
then
  if [ $MSHFILES -eq 1 ] || [ $MSHFILES -eq 3 ]
  then
    printf "Removing iteration %d intermediate files..."  $(($ITERATIONS-1))
    printf -v mshFileRm '%s%s_%d.msh' "$MSHDIR"  "$MODEL" $(($ITERATIONS-1))
    printf -v mshLogRm  '%s.log'      "$mshFileRm"
    printf -v posFileRm '%s%s_%d.pos' "$MSHDIR"  "$MODEL" $(($ITERATIONS-1))
    printf -v posLogRm  '%s.log'      "$posFileRm"
    rm -f $mshFileRm $mshLogRm $posFileRm $posLogRm
    printf " Done!\n"
  fi
fi

# Print info on the run
echo "$(date) Finished mesh refinement."