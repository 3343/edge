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
# Job scheduling.
##
import saga
import logging
import time
import sys

class Schedule:
  # jobs
  m_jobs = { 'new':      [],
             'queued':  [],
             'done':     [],
             'failed':   [] }

  ##
  # Error handling of SAGA-exceptions.
  #
  # @param i_ex exception.
  ##
  def sagaEx( self,
              i_ex ):
    logging.fatal( 'Saga exception: ' + i_ex.type + ' ' + str(i_ex ) )
    logging.fatal( 'Jobs: ' + str(self.m_jobs) )
    logging.fatal( 'Backtrace: ' + str(i_ex.traceback) )

  ##
  # Constructor, which initializes the scheduler and corresponding machine config.
  #
  # @param i_arch architecture of the machine.
  # @param i_project project under which the jobs are run.
  # @param i_user user id.
  # @param i_serviceUrl service url.
  # @param i_env environment variables set for jobs.
  # @param i_maxJobs maximum number of jobs submitted in parallel.
  ##
  def __init__( self,
                i_arch,
                i_project,
                i_user,
                i_serviceUrl,
                i_env,
                i_maxJobs ):
    try:
      # set architecture
      self.m_architecture = i_arch

      # set project
      self.m_project = i_project

      # set environment
      self.m_environment = i_env

      # set maximum number of jobs
      self.m_maxJobs = i_maxJobs

      # create context
      self.m_context = saga.Context('ssh')
      self.m_context.user_id = i_user

      # create saga session
      self.m_session = saga.Session()
      self.m_session.add_context( self.m_context )

      # create service
      self.m_service = saga.job.Service( i_serviceUrl , session = self.m_session )
    
    # catch saga exceptions
    except saga.SagaException, i_ex:
      self.sagaEx( i_ex )

  ##
  # Runs the longest queued job.
  ##
  def run( self ):
    if len( self.m_jobs['new'] ) > 0:
      try:
        # add to service
        self.m_jobs['queued'] = self.m_jobs['queued'] + \
          [ { 'tag':  self.m_jobs['new'][0]['tag'],
              'saga': self.m_service.create_job( self.m_jobs['new'][0]['saga'] )
            } ]

        # run job
        self.m_jobs['queued'][-1]['saga'].run()
        logging.debug( '  queued: ' + \
                      self.m_jobs['queued'][-1]['tag'] + ', ' +\
                      self.m_jobs['queued'][-1]['saga'].id )

        # remove job from new
        self.m_jobs['new'] = self.m_jobs['new'][1:]
      except saga.SagaException, i_ex:
        self.sagaEx( i_ex )
  
  ##
  # Clears the queue by moving finished jobs to finished.
  ##
  def clear( self ):
    for l_ru in self.m_jobs['queued']:
      l_state = l_ru['saga'].get_state()

      if l_state in [ saga.job.UNKNOWN,
                      saga.job.NEW,
                      saga.job.PENDING,
                      saga.job.RUNNING ]:
        pass # white list, nothing to do
      elif l_state == saga.job.DONE:
        logging.debug( '  done: ' + l_ru['tag'] + ', ' + l_ru['saga'].id )
        self.m_jobs['done'] = self.m_jobs['done'] + [ l_ru ]
        self.m_jobs['queued'].remove( l_ru )
      
      elif l_state in [ saga.job.FAILED,
                        saga.job.CANCELED,
                        saga.job.SUSPENDED ]:
        logging.warning( '  failed: ' + l_ru['tag'] + ', ' + l_ru['saga'].id )
        self.m_jobs['failed'] = self.m_jobs['failed'] + [ l_ru ]
        self.m_jobs['queued'].remove( l_ru )
      else:
        logging.fatal( '  dont know nothing about job the state: ' + l_state )

  ##
  # Generates a report.
  ##
  def report( self ):
    logging.info( '  status report coming through.. #new: '     + str( len(self.m_jobs['new']) ) +\
                                                 ', #queued: '  + str( len(self.m_jobs['queued']) ) +\
                                                 ', #done: '    + str( len(self.m_jobs['done']) ) +\
                                                 ', #failed: '  + str( len(self.m_jobs['failed']) ) )

  ##
  # Runs and waits until all jobs are finished.
  #
  # @param i_timeOut time out times between job-status checks.
  # @param i_report generate status report after the given number of timeouts.
  ##
  def runWaitAll( self,
                  i_timeOut=30,
                  i_report=2 ):
    try:
      logging.info( '  waiting for all jobs to finish')

      # report counter
      l_re = 0

      while len( self.m_jobs['queued'] ) > 0 or len( self.m_jobs['new'] ) > 0:
        # clear queued jobs
        self.clear()

        # add jobs if possible
        if len( self.m_jobs['queued'] ) < self.m_maxJobs:
          for l_re in range( self.m_maxJobs - len(self.m_jobs['queued']) ):
            self.run()

        if l_re % i_report == 0:
          self.report()
        l_re = l_re+1

        # sleep
        time.sleep( i_timeOut )

      # one final report
      self.report()
    
    # catch interrupts and try to clean up the jobs
    except KeyboardInterrupt:
      # get number of queued jobs
      l_nr = len( self.m_jobs['queued'] )

      # try to abort jobs
      logging.info( 'caught interrupt, trying to cancel ' + str( l_nr ) + ' jobs'  )
      for l_ru in self.m_jobs['queued']:
        l_ru['saga'].cancel()

      # exit
      logging.info( 'finished cancelling jobs, exiting now' )
      sys.exit(-1)

  ##
  # Adds a job to the scheduler.
  #
  # @param i_queue queue used for execution of the job.
  # @param i_nNodes number of nodes requested.
  # @param i_timeLim time limit of the job.
  # @param i_env environemnt used during job execution.
  # @param i_wd working directory.
  # @param i_exe executable of the job.
  # @param i_args arguments of the job.
  # @param i_tag tag of the job, used for job-logs.
  ##
  def addJob( self,
              i_queue,
              i_nNodes,
              i_timeLim,
              i_env,
              i_wd,
              i_exe,
              i_args,
              i_tag ):
    # create job description
    l_jd = saga.job.Description()

    # scheduler specific variables
    l_jd.cpu_architecture = self.m_architecture
    l_jd.project = self.m_project

    # set job-specific variables
    if len(i_env + self.m_environment ) > 0: l_jd.environment = dict( i_env + self.m_environment )

    l_jd.executable = i_exe
    l_jd.working_directory = i_wd
    l_jd.arguments = i_args

    l_jd.number_of_processes = i_nNodes
    l_jd.processes_per_host = 1 # we request nodes as smallest entity

    l_jd.queue = i_queue
    l_jd.wall_time_limit = i_timeLim
    l_jd.name = i_tag
    l_jd.output = i_tag + '.out'
    l_jd.error = i_tag + '.err'

    l_job = { 'tag': i_tag,
              'saga': l_jd  }

    # enqueue job
    self.m_jobs['new'] = self.m_jobs['new'] + [ l_job ]
