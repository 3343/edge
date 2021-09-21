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
# The extracted files and the plots will be placed inside a "plots" folder in the same directory as the error directory
# Usage ./error_norms_convergence_plots.sh /path/to/error/folder NumberOfQuantities NumberOf
# 
# @section DEPENDENCIES
# This file depends on the following projects: 
# xmlstarlet: http://xmlstar.sourceforge.net/
# gnuplot: http://www.gnuplot.info/index.html


echo "Hi, I am starting now!"

# remove and make new files

rm -r norms
mkdir norms

for norm_val in '1' '2' 'inf'
do
	for quantity_val in $(seq 1 5)
	do 
		echo "#L${norm_val} Data Quantity ${quantity_val}" > norms/l${norm_val}_${quantity_val}.csv
	done
done

#echo '#L1 Data Quantity 1' > l1_1.txt
#echo '#L1 Data Quantity 1' > l1_2.txt
#echo '#L1 Data Quantity 1' > l1_3.txt
#echo '#L1 Data Quantity 1' > l1_4.txt
#echo '#L1 Data Quantity 1' > l1_5.txt
#echo '#L2 Data' > l2.txt<
#echo '#Linf Data' > linf.txt

# find relevant data and extract into different files 

for timestepping in gts #lts
do
	for cl in 25 20 15 10 9 8 7 6 5 4 3 2 1
	do
		for parallel in 1 #13
		do 
			echo errors/${timestepping}_cl_${cl}_pa_${parallel}.xml
			for i in $(seq 1 5)
			do
				echo ${cl}, $(xmlstarlet sel -t -v "error_norms/l1/q[${i}]" -nl errors/${timestepping}_cl_${cl}_pa_${parallel}.xml) >> norms/l1_${i}.csv
				echo ${cl}, $(xmlstarlet sel -t -v "error_norms/l2/q[${i}]" -nl errors/${timestepping}_cl_${cl}_pa_${parallel}.xml) >> norms/l2_${i}.csv
				echo ${cl}, $(xmlstarlet sel -t -v "error_norms/linf/q[${i}]" -nl errors/${timestepping}_cl_${cl}_pa_${parallel}.xml) >> norms/linf_${i}.csv
			done
		done
	done 
done 

# plot data to pdf with gnuplot
gnuplot -e "set terminal pdf ;
set output 'norms/output.pdf'; 
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
set style line 5 \
    linecolor rgb '#ff30ff' \
    linetype 1 linewidth 2 \
    pointtype 4 pointsize 1.5;
set multiplot; 
set yrange [0:3700];
set datafile separator comma;
plot [25:0] 'norms/l1_1.csv' with linespoints linestyle 1; 
plot [25:0] 'norms/l1_2.csv' with linespoints linestyle 2; 
plot [25:0] 'norms/l1_3.csv' with linespoints linestyle 3; 
plot [25:0] 'norms/l1_4.csv' with linespoints linestyle 4; 
plot [25:0] 'norms/l1_5.csv' with linespoints linestyle 5; 
unset multiplot;
set multiplot; 
set yrange [0:45];
set datafile separator comma;
plot [25:0] 'norms/l2_1.csv' with linespoints linestyle 1; 
plot [25:0] 'norms/l2_2.csv' with linespoints linestyle 2; 
plot [25:0] 'norms/l2_3.csv' with linespoints linestyle 3; 
plot [25:0] 'norms/l2_4.csv' with linespoints linestyle 4; 
plot [25:0] 'norms/l2_5.csv' with linespoints linestyle 5; 
unset multiplot;
set multiplot; 
set yrange [0:0.8];
set datafile separator comma;
plot [25:0] 'norms/linf_1.csv' with linespoints linestyle 1; 
plot [25:0] 'norms/linf_2.csv' with linespoints linestyle 2; 
plot [25:0] 'norms/linf_3.csv' with linespoints linestyle 3; 
plot [25:0] 'norms/linf_4.csv' with linespoints linestyle 4; 
plot [25:0] 'norms/linf_5.csv' with linespoints linestyle 5; 
unset multiplot ;
"
