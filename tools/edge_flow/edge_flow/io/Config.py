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
# XML configuration.
##
import xml.etree.ElementTree
import logging
import itertools
import copy

class Config:
  ##
  # Creates template instantiations by replacing template parameters in the job template.
  #
  # @param i_firstId id of the first job
  # @param i_tmplPars template parameters.
  # @param i_jobTmpl job template
  # @return instantiated jobs. 
  ##
  def instJobTmpl( self,
                   i_firstId,
                   i_tmplPars,
                   i_jobTmpl ):
    # job id
    l_jobId = i_firstId

    # assemble list of template tags which get replaced
    l_tags = [ x['tag'] for x in i_tmplPars]

    # assemble list of replacement values
    l_vals = [ x['values'] for x in i_tmplPars ]

    # compute cartersian product of values
    l_vals = list( itertools.product( *l_vals ) )

    # generate jobs by replaces values
    l_jobs = []

    for l_va in l_vals:
      l_job = copy.deepcopy( i_jobTmpl )

      #  iterate over all keys
      for l_ke in l_job.keys():
        # replace values
        for l_ta in range( len(l_tags) ):
          if l_job[l_ke] != None:
            # try key as value
            try:
              l_job[l_ke] = l_job[l_ke].replace( l_tags[l_ta], l_va[l_ta] )
            # try key as list
            except AttributeError:
              for l_it in range( len(l_job[l_ke]) ):
                l_job[l_ke][l_it] = l_job[l_ke][l_it].replace( l_tags[l_ta], l_va[l_ta] )

      # convert number of nodes to int
      if l_job['nodes'] != None:
        l_job['nodes'] = int( l_job['nodes'] )

      # set job-id and increase counter
      l_job['id'] = l_jobId
      l_jobId = l_jobId + 1

      # add job
      l_jobs = l_jobs + [l_job]

    return l_jobId, l_jobs

  ##
  # Prints the configuration.
  ##
  def log( self ):
    logging.info( 'EDGEflow\'s config:' )
    logging.info( '  machines:'        )
    for l_ma in self.m_machines.keys():
      logging.info( '    ' + l_ma + ':' )
      for l_en in ['architecture', 'project', 'user', 'service_url']:
        if self.m_machines[l_ma][l_en] != None:
          logging.info( '      ' + l_en + ': ' + self.m_machines[l_ma][l_en] )
      logging.info( '      max_jobs: ' + str( self.m_machines[l_ma]['max_jobs'] ) )
      logging.info( '      environment: ' )
      for l_var in self.m_machines[l_ma]['environment']:
        logging.info( '        ' + l_var[0] + '=' + l_var[1] )

    logging.info(   '  stages:' )
    for l_st in self.m_stages:
      logging.info( '    ' + l_st['tag'] + ':' )
      for l_en in ['machine']:
        if l_st[l_en] != None:
          logging.info( '      ' + l_en + ': ' + l_st[l_en] )
      logging.info( '      #jobs: ' + str(len(l_st['jobs'])) )

  ##
  # Constructor, which parses the XML-config.
  #
  # @param i_pathToXml xml config of the workflow.
  ##
  def __init__(self, i_pathToXml):
    # parse the xml
    l_xml = xml.etree.ElementTree.parse( i_pathToXml )

    # get root
    l_conf = l_xml.getroot()

    #
    #  read machine info
    #
    self.m_machines = {}
    l_machs = list( l_conf.find( 'machines' ) )

    # iterate over machines
    for l_ma in l_machs:
      # create dict
      self.m_machines[l_ma.tag] = {}

      # popular dict
      for l_en in ['architecture', 'project', 'user', 'service_url']:
        if l_ma.find(l_en) != None:
          self.m_machines[l_ma.tag][l_en] = l_ma.find(l_en).text
        else:
          self.m_machines[l_ma.tag][l_en] = None

      self.m_machines[l_ma.tag]['max_jobs'] = int( l_ma.find('max_jobs').text )

      l_env = l_ma.find( 'environment' )
      self.m_machines[l_ma.tag]['environment'] = []
      if l_env != None:
        for l_var in l_env.findall('variable'):
          self.m_machines[l_ma.tag]['environment'] = self.m_machines[l_ma.tag]['environment'] + \
            [ ( l_var.find('name').text, l_var.find('value').text ) ]

    #
    # read stages
    #
    l_stages = l_conf.find( 'stages' )
    # create stages
    self.m_stages = []
    
    # unique job-id
    l_jobId = 0

    # iterate over stages
    for l_st in list(l_stages):
      # create temporary directory for the stage
      l_stage = { 'tag': l_st.tag }

      # stage propoerties
      for l_en in ['machine']:
        if l_st.find(l_en) != None:
          l_stage[l_en] = l_st.find(l_en).text
        else: l_stage[l_en] = None

      # template parameters
      l_tmplPars = []
      for l_tp in l_st.findall('template_parameter'):
        # assemble single template parameter
        l_tmplPar = { 'tag': l_tp.find('tag').text, 'values': [] }
        # read values
        for l_va in l_tp.findall('value'):
          l_tmplPar['values'] = l_tmplPar['values'] + [ l_va.text ]
        # save template parameter
        l_tmplPars = l_tmplPars + [ l_tmplPar ]

      # create jobs
      l_stage['jobs'] = []
      # iterate over job-templats
      for l_jo in list( l_st.find('jobs') ):
        l_jobTmpl = { 'id': None,
                      'tag': l_jo.tag,
                      'arguments': [] }

        # job properties
        for l_en in ['time_limit', 'executable', 'config', 'queue', 'nodes', 'working_dir']:
          if l_jo.find(l_en) != None:
            l_jobTmpl[l_en] = l_jo.find(l_en).text
          else: l_jobTmpl[l_en] = None

        # arguments
        for l_ar in l_jo.findall('argument'):
          l_jobTmpl['arguments'] = l_jobTmpl['arguments'] + [ l_ar.text ]

        # create non-template job
        if len(l_tmplPars) == 0:
          # set job-id
          l_jobTmpl['id'] = l_jobId

          # #nodes are int
          if l_jobTmpl['nodes'] != None:
            l_jobTmpl['nodes'] = int(l_jobTmpl['nodes'])

          l_stage['jobs'] = l_stage['jobs'] + [ l_jobTmpl ]

          # increase job-id
          l_jobId = l_jobId + 1

        # create templated jobs
        else:
          l_jobId, l_jobs = self.instJobTmpl( l_jobId, l_tmplPars, l_jobTmpl )
          l_stage['jobs'] = l_stage['jobs'] + l_jobs

      # append stage
      self.m_stages = self.m_stages + [ l_stage ]

    # log config
    self.log()