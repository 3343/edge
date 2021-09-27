#!/bin/bash
##
# @file This file is part of EDGE.
#
# @author Sarah Bachinger (sarah.bachinger AT uni-jena.de)
#
# @section LICENSE
# Copyright (c) 2020, Friedrich Schiller University Jena
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
# A script to generate convergence plots from the error norms of an EDGE computation.
# Identifies error norms in xml files and extracts them for plotting according to the error norm. 
# The extracted files and the plots will be placed inside a "plots" folder in the same directory as the error directory.
# MAKE SURE there is no other directory named "plots" in there!
# The goal is to use display one plot for every error norm and quantity and display the values for gts, lts, parallel and not parallel.
#
# Usage ./error_norms_convergence_plots.sh /path/to/error/folder NumberOfQuantities 
# 
# @section DEPENDENCIES
# This file depends on the following projects: 
# xmlstarlet: http://xmlstar.sourceforge.net/
# gnuplot: http://www.gnuplot.info/index.html


echo "Hi, I am starting now!"

# remove and make new files
# $1 is /path/to/error/folder
# $2 is NumberOfQuantities

#rm -r -f $1/plots
#mkdir $1/plots

# find relevant data and extract into different files  


echo 'start extracting ... '

#for timestepping in 'gts' 'lts'
: '
do
	for cl in 25 20 15 10 9 8 7 6 5 4 3 2 1
	do
		for parallel in 1 13
		do 
			#echo $1/errors/${timestepping}_cl_${cl}_pa_${parallel}.xml
			for i in $(seq 1 $2)
			do
				echo ${cl}, $(xmlstarlet sel -t -v "error_norms/l1/q[${i}]" -nl $1/errors/${timestepping}_cl_${cl}_pa_${parallel}.xml) >> $1/plots/l1_${i}_${timestepping}_pa_${parallel}.csv
				echo ${cl}, $(xmlstarlet sel -t -v "error_norms/l2/q[${i}]" -nl $1/errors/${timestepping}_cl_${cl}_pa_${parallel}.xml) >> $1/plots/l2_${i}_${timestepping}_pa_${parallel}.csv
				echo ${cl}, $(xmlstarlet sel -t -v "error_norms/linf/q[${i}]" -nl $1/errors/${timestepping}_cl_${cl}_pa_${parallel}.xml) >> $1/plots/linf_${i}_${timestepping}_pa_${parallel}.csv
			done
		done
	done 
done 
'

echo 'All data extracted... start plotting'


# plot data to pdf with gnuplot
gnuplot -e "set terminal pdf ;
set output '$1/plots/output.pdf'; 
set style line 1 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.5 ;
set style line 2 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 2 \
    pointtype 5 pointsize 1.5;
set style line 3 \
    linecolor rgb '#51ff2e' \
    linetype 1 linewidth 2 \
    pointtype 2 pointsize 1.5;
set style line 4 \
    linecolor rgb '#132f38' \
    linetype 1 linewidth 2 \
    pointtype 3 pointsize 1.5;
do for [l in \"1 2 inf\" ] {
	do for [ q = 1:$2]{
		set multiplot title 'L'.l.' Q'.q; 
		set datafile separator comma;
		set grid;
		set logscale xy;
		set size ratio 0.5;
		set ylabel 'error';
		set xlabel 'mesh width';
		plot '$1/plots/l'.l.'_'.q.'_gts_pa_1.csv' with linespoints linestyle 2 title 'L'.l.' Q'.q.' gts non parallel', \
		'$1/plots/l'.l.'_'.q.'_gts_pa_13.csv' with linespoints linestyle 4 title 'L'.l.' Q'.q.' gts parallel', \
		'$1/plots/l'.l.'_'.q.'_lts_pa_1.csv' with linespoints linestyle 3 title 'L'.l.' Q'.q.' lts non parallel', \
		'$1/plots/l'.l.'_'.q.'_lts_pa_13.csv' with linespoints linestyle 1 title 'L'.l.' Q'.q.' lts parallel'; 
		unset logscale;
		unset multiplot;
	}
}
"
echo 'finished plotting'