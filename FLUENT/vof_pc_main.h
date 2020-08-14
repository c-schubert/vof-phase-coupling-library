/* 
vof_pc_main.h
Main header file for the VoF phase coupling library for coupling between 
ANSYS Fluent and ANSYS Mechanical APDL.

You have to modifiy vof_pc_case.h and vof_pc_nn_coupling.c to fit your case!

License (MIT):

Copyright (c) 2016-2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#ifndef VOF_PC_H

#define VOF_PC_H
#include "udf.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "malloc.h"
#include "udf_helpers.h"

/* Enable more I/O messages for debbuging purposes */
#define VOF_PC_DEBUG 0

#define VOF_MAX_REL_CHANGE 0.25 /* Max value that any VOF of a Cell in Fluent can cange until recoupling with ANSYS if loose coupling is choosen */

#define MAX_COUPLING_TRIALS 10000 
#define COUPLING_SLEEP_TIME_IN_S 1

enum couplingStates {COUPLING_INIT=0, ANSYS_READY, FLUENT_READY, STOP_SIM, SYNC_ERROR};
enum coupledProperties {NONE=0, JOULE_HEAT, JOULE_HEAT_PLUS_LORENTZ, VOF};
enum udmis {UDM_JH=0, UDM_LFx, UDM_LFy, UDM_LFz, UDM_VOF_old};

#define LINUX 0

#if LINUX
#include "unistd.h"
#endif

#endif
