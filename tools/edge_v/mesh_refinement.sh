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
# Tested on: Ubuntu 17.10
##


# Global variables
declare -g iter=0
declare -g iterations=10
declare -g model='PLACEHOLDER'
declare -g geoFile
declare -g mshFile
declare -g mshLog
declare -g posFile

# Affected files and directories
confDir='PLACEHOLDER'
mshDir='PLACEHOLDER'
writePos='PLACEHOLDER'

# ssh variables (use remoteMsh=1 to mesh on a remote client; by default 0)
# IMPORTANT: use only if you have a working SSH public key with the remote client!!
declare -g remoteMsh=1
declare -g remoteUsr='PLACEHOLDER'
declare -g remoteClient='PLACEHOLDER'
declare -g remoteGmshDir='PLACEHOLDER'
declare -g remoteMshDir='PLACEHOLDER'


# Function to send geo file from local m/c to remote client
function sendGeo() {
  printf -v l_login '%s@%s' "$remoteUsr" "$remoteClient"

  # Copy geo file to remote client
  printf "\nSending $geoFile from local m/c to remote client ...\n"
  printf -v l_remoteLoc '%s:%s' "$l_login" "$remoteMshDir"

  SECONDS=0
  scp -v $geoFile $l_remoteLoc
  duration=$SECONDS
  printf "Time taken: $duration s\n\n"
}


# Function to mesh on remote client
function mesh() {
  declare bgmPos=''

  if [ $remoteMsh -eq 1 ]
  then
    # Create remote msh dir (if doesn't exist)
    printf -v l_login '%s@%s' "$remoteUsr" "$remoteClient"
    ssh $l_login mkdir -p $remoteMshDir

    # Remote location of geo file
    printf -v l_geoFile '%s%s.geo'           "$remoteMshDir" "$model"

    # Remote location of mesh file
    if [ $iter -eq 0 ]
    then
      printf -v l_mshFile '%s%s.msh'         "$remoteMshDir" "$model"
      printf -v l_posFile '%s%s.pos'         "$remoteMshDir" "$model"
    elif [ $iter -eq $iterations ]
    then
      printf -v l_mshFile '%s%s_refined.msh' "$remoteMshDir" "$model"
      printf -v l_posFile '%s%s_refined.pos' "$remoteMshDir" "$model"
    else
      printf -v l_mshFile '%s%s_%d.msh'      "$remoteMshDir" "$model" "$iter"
      printf -v l_posFile '%s%s_%d.pos'      "$remoteMshDir" "$model" "$iter"
    fi

    # Send pos file from local to remote
    if [ $iter -ne 0 ]
    then
      printf "\nSending $posFile from local m/c to remote client ...\n"

      SECONDS=0
      scp -v $posFile $l_remoteLoc
      duration=$SECONDS

      printf "Time taken: $duration s\n\n"
      printf -v bgmPos '%s %s' "-bgm" "$l_posFile"
    fi

    # Remote machine mesh
    ssh $l_login nohup $remoteGmshDir $l_geoFile $bgmPos -3 -o $l_mshFile 2>&1 | tee $mshLog

    # Copy mesh file back to local machine
    printf "\nSending $l_mshFile from remote client to local m/c ...\n"
    printf -v l_retMshFile '%s:%s' "$l_login" "$l_mshFile"

    SECONDS=0
    scp -v $l_retMshFile $mshDir
    duration=$SECONDS
    printf "Time taken: $duration s\n\n"
  else
    # Mesh on local machine
    if [ $iter eq 0 ]
    then
      gmsh $geoFile -3 -o $mshFile 2>&1 | tee $mshLog
    else
      gmsh $geoFile -bgm $posFile -3 -o $mshFile 2>&1 | tee $mshLog
    fi
  fi
}


#######################
##### Entry Point #####
#######################

# Affected files
printf -v confFile '%s%s.conf' "$confDir" "$model"
printf -v geoFile  '%s%s.geo'  "$mshDir"  "$model"
printf -v mshFile  '%s%s.msh'  "$mshDir"  "$model"

# Turn off optimization for initial (coarse) and intermediate meshes
sed -i -e '/Mesh.Optimize  /c\Mesh.Optimize       = 0;'     $geoFile
sed -i -e '/Mesh.OptimizeNetgen/c\Mesh.OptimizeNetgen = 0;' $geoFile

sendGeo

# Intial coarse mesh
printf "Generating initial (coarse) mesh ...\n"
printf -v mshLog '%s.log' "$mshFile"
mesh

for ((iter=1; iter < $iterations; iter++))
do
  printf "###########\nIteration $iter\n###########\n\n"

  # Pos file name for current iteration
  printf -v posFile '%s%s_%d.pos' "$mshDir"  "$model" "$iter"

  # Replace input mesh_file and output pos_file in config file
  printf -v msh '/mesh_file=/c\mesh_file=%s' "$mshFile"
  sed -i -e "$msh" $confFile
  
  printf -v pos '/pos_file=/c\pos_file=%s'   "$posFile"
  sed -i -e "$pos" $confFile

  # Mesh file (to be generated) which will be used in next iteration
  printf -v mshFile '%s%s_%d.msh' "$mshDir"  "$model" "$iter"
  printf -v mshLog  '%s.log'      "$mshFile"

  # Generate intermediate pos file to be used by gmsh in next step to produce a 
  # more refined mesh than previous iteration (using pos as background mesh)
  printf -v posLog '%s.log' "$posFile"
  printf "Generating pos ...\n"
  $writePos -f $confFile 2>&1 | tee $posLog

  # Generate intermediate meshes
  printf "\nGenerating intermediate ($iter) mesh ...\n"
  mesh
done

# Turn on optimization for final (refined) mesh
sed -i -e '/Mesh.Optimize  /c\Mesh.Optimize       = 1;'     $geoFile
sed -i -e '/Mesh.OptimizeNetgen/c\Mesh.OptimizeNetgen = 1;' $geoFile

sendGeo

printf -v posFile '%s%s_refined.pos' "$mshDir" "$model"

printf -v msh '/mesh_file=/c\mesh_file=%s' "$mshFile"
sed -i -e "$msh" $confFile

printf -v pos '/pos_file=/c\pos_file=%s'   "$posFile"
sed -i -e "$pos" $confFile

printf -v mshFile '%s%s_refined.msh' "$mshDir" "$model"
printf -v mshLog '%s_refined.log' "$mshFile"

# Final refined pos file
printf -v posLog '%s_refined.log' "$posFile"
printf "Generating pos ...\n"
$writePos -f $confFile 2>&1 | tee $posLog

# Final mesh-refinement
printf '\n\nGenerating final (refined) mesh ...\n'
mesh