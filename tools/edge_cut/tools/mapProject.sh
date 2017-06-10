#!/bin/bash
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
# Extracts the given topo and performs a map projection.
##

poi=data/poi.csv
tools_dir=../../../tools/meshing/
region_proj=-120/-116/33/35
region_cut=330000/500000/3700000/3850000
tmp_dir=tmp
#grid=grids/gebco_2014/GEBCO_2014_2D.nc
grid=grids/southern_calif_crm_v1.nc

mkdir $tmp_dir

# cut region of interest
gmt grdcut -R$region_proj $grid -G$tmp_dir/region.nc

# generate a color tagble
gmt makecpt -Cglobe -T-10000/10000/200 > $tmp_dir/colors.cpt

# plot topography 2D
gmt grdimage -C$tmp_dir/colors.cpt -JM13c $tmp_dir/region.nc -Ba -K -C$tmp_dir/colors.cpt > $tmp_dir/region.ps
#gmt grdcontour -C1000 -JM13c $tmp_dir/region.nc -Ba -O >> $tmp_dir/region.ps

# add points of interests
gmt psxy $poi -R"$region_proj" -JM13c -Sd.3 -Gred -O -K -V >> $tmp_dir/region.ps
gmt pstext $poi -R"$region_proj" -JM13c -F+f7+jL -D0.2/0c -O -V >> $tmp_dir/region.ps

# write poi's to ascii
gmt mapproject $poi -R$region_proj -J$JM13c -F -V > tmp/poi.txt

# plot topo 3D
gmt grdview $tmp_dir/region.nc -JM15c -p155/35 -Qi250 -C$tmp_dir/colors.cpt -Ba -JZ2c > $tmp_dir/region_pers.ps

# convert to PDF
mkdir output
ps2pdf $tmp_dir/region.ps output/region.pdf
ps2pdf $tmp_dir/region_pers.ps output/region_pers.pdf
mv $tmp_dir/poi.txt output/poi.txt

# do the map projection
gmt grdproject $tmp_dir/region.nc -Jutm/+11/1:1 -C -G$tmp_dir/map_proj.nc -F -V2
gmt grdsample $tmp_dir/map_proj.nc -I550 -R${region_cut} -G$tmp_dir/sampled.nc -V2
gmt grd2xyz $tmp_dir/sampled.nc > output/sampled.xyz
