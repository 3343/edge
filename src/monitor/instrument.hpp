/**
 * @file This file is part of EDGE.
 *
 * @author Alexander Breuer (anbreuer AT ucsd.edu)
 *
 * @section LICENSE
 * Copyright (c) 2016-2018, Regents of the University of California
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Instrumentation macros.
 **/

#ifndef EDGE_MONITOR_INSTRUMENT_H_
#define EDGE_MONITOR_INSTRUMENT_H_

#ifdef PP_USE_INSTR

#include <scorep/SCOREP_User.h>

#ifdef __GNUC__
// gcc throws unused variable errors and the ignored diagnostic doesn't work (5.4.1) -> touch manually
class gccDisableUnusedScoreP {
  void dummy() { SCOREP_User_LastFileHandle+=0; SCOREP_User_LastFileName+=0; }
};
#endif

// forward scorep-macros for instrumentation
#define PP_INSTR_FUN(str)              SCOREP_USER_REGION(str,SCOREP_USER_REGION_TYPE_FUNCTION)
#define PP_INSTR_REG_DEF(str)          SCOREP_USER_REGION_DEFINE(str)
#define PP_INSTR_REG_BEG(str1,str2)    SCOREP_USER_REGION_BEGIN(str1,str2,SCOREP_USER_REGION_TYPE_COMMON)
#define PP_INSTR_PAR_UINT64(str1,str2) SCOREP_USER_PARAMETER_UINT64(str1,str2)
#define PP_INSTR_REG_END(str)          SCOREP_USER_REGION_END(str)
#define PP_INSTR_REG_NAME_BEG(str)     SCOREP_USER_REGION_BY_NAME_BEGIN(str,SCOREP_USER_REGION_TYPE_COMMON)
#define PP_INSTR_REG_NAME_END(str)     SCOREP_USER_REGION_BY_NAME_END(str)

#else

// set empty macros in the case of disabled instrumentation
#define PP_INSTR_FUN(str)
#define PP_INSTR_REG_DEF(str)
#define PP_INSTR_REG_BEG(str1,str2)
#define PP_INSTR_PAR_UINT64(str1,str2)
#define PP_INSTR_REG_END(str)
#define PP_INSTR_REG_NAME_BEG(str)
#define PP_INSTR_REG_NAME_END(str)

#endif

#endif
