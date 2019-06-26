/* 
Main header file for the VoF phase coupling library for coupling between 
ANSYS Fluent and ANSYS Mechanical APDL.


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

#define FLUID_ID 11 /* Adopt!*/

#define MAX_COUPLING_TRIALS 1200 
#define COUPLING_SLEEP_TIME_IN_S 1

#define _XC_FOLDER_PATH_ "Z:\\ESU\\Cases\\Coupling\\Neue_Kopplung\\XC\\" /* Adopt! */

#define _SYNC_DAT_ _XC_FOLDER_PATH_ "SYNC.TXT"
#define _FLUENT_TO_ANSYS_VOFOUT_DAT_ _XC_FOLDER_PATH_ "FLUENT_TO_ANSYS_VOF_OUT.DAT"
#define _FLUENT_ALLOUT_DAT_ _XC_FOLDER_PATH_ "FLUENT_DEBUG_ALL_OUT.DAT"
#define _ANSYS_TO_FLUENT_OUT_DAT_ _XC_FOLDER_PATH_ "ANSYS_TO_FLUENT_OUT.DAT"
#define _ANSYS_ALLOUT_DAT_ _XC_FOLDER_PATH_ "ANSYS_TO_FLUENT_ALLOUT.DAT"
#define _ANSYS_TO_FLUENT_MAPPING_DAT_ _XC_FOLDER_PATH_ "ANSYS_TO_FLUENT_MAP.DAT"
#define _FLUENT_TO_ANSYS_MAPPING_DAT_ _XC_FOLDER_PATH_ "FLUENT_TO_ANSYS_MAP.DAT"


enum couplingStates {COUPLING_INIT=0, ANSYS_READY, FLUENT_READY, STOP_SIM, SYNC_ERROR};
enum udmis {Jouleheat=0, LFx, LFy, LFz};

#endif
