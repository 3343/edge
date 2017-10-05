#!/usr/bin/env python
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
# This is the entry point of EDGEflow.
##
import edge_flow.jobs.Schedule
import logging
import argparse
import edge_flow.io.Config

# set up logger
logging.basicConfig( level=logging.INFO,
                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' )

# welcome our users
logging.info( '##########################################################################' )
logging.info( '##############   ##############            ###############  ##############' )
logging.info( '##############   ###############         ################   ##############' )
logging.info( '#####            #####       #####      ######                       #####' )
logging.info( '#####            #####        #####    #####                         #####' )
logging.info( '#############    #####         #####  #####                  #############' )
logging.info( '#############    #####         #####  #####      #########   #############' )
logging.info( '#####            #####         #####  #####      #########           #####' )
logging.info( '#####            #####        #####    #####        ######           #####' )
logging.info( '#####            #####       #####      #####       #####            #####' )
logging.info( '###############  ###############         ###############   ###############' )
logging.info( '###############  ##############           #############    ###############' )
logging.info( '######################################################################flow' )
logging.info( 'you reached the workflow engine of EDGE' )

# parse command line options
l_parser = argparse.ArgumentParser( description='Execution of workflows in EDGE.' )

l_parser.add_argument( '-x', '--xml',
                       dest     = 'xml',
                       required = True,
                       type     = str,
                       help     = 'XML configuration of EDGEflow')
l_args = vars( l_parser.parse_args() )

# parse XML-config
l_conf = edge_flow.io.Config.Config( l_args['xml'] )

# create job schedulers
l_scheds = {}
for l_k in l_conf.m_machines.keys():
  l_ma = l_conf.m_machines[l_k]
  l_scheds[l_k] = edge_flow.jobs.Schedule.Schedule( l_ma['architecture'],
                                                    l_ma['project'],
                                                    l_ma['user'],
                                                    l_ma['service_url'],
                                                    l_ma['environment'],
                                                    l_ma['max_jobs'] )

# iterate over stages
for l_st in l_conf.m_stages:
  logging.info( 'running stage: ' + l_st['tag'] )

  # iterate over jobs and add them to the scheduler
  for l_jo in l_st['jobs']:
    # assemble base name
    l_tag = l_st['tag'] + '_' + l_jo['tag'] + '_' + str(l_jo['id'])

    # add job to scheduler
    l_scheds[ l_st['machine'] ].addJob( l_jo['queue'],
                                        l_jo['nodes'],
                                        l_jo['time_limit'],
                                        [], # TODO: add job-specific env
                                        l_jo['working_dir'],
                                        l_jo['executable'],
                                        l_jo['arguments'],
                                        l_tag )
  # run job of stage's scheduler
  l_scheds[ l_st['machine'] ].runWaitAll()

logging.info( 'EDGEflow finished, see you next time!' )
