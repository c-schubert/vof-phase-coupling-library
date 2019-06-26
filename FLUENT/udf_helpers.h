/*
 Collection of divers and generally usefull functions
 for makink writing udf functions a little bit more comfortable

License (MIT):

Copyright (c) 2016-2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef MY_UDS_HELPERS_H
#define MY_UDS_HELPERS_H
#include "udf.h"

/* Defined error states */

#define _STATE_OK 0
#define _STATE_WARNING 1
#define _STATE_ERROR -1

/* Settings */

/* Tolerance values for numerical comparison of floating point numbers*/
#define MAC_ZERO 1e-10 
#define MAC_SMALL 1e-4 

#define MAC_IS_ZERO(A) ( (A < MAC_ZERO) && (A > -MAC_ZERO) )
#define MAC_IS_SMALL(A) ( (A < MAC_SMALL) && (A > -MAC_SMALL) )
#define MAC_IS_EQUAL(A, B) ( ((A-B) < MAC_ZERO) && ((A-B) > -MAC_ZERO) )
#define MAC_IS_SIMILAR(A, B) ( ((A-B) < MAC_SMALL) && ((A-B) > -MAC_SMALL) )
#define MAC_IS_IN_BETWEEN(A, B, C) ( ((A-B) > -MAC_ZERO) && ((A-C) < MAC_ZERO) )

/* Enable boolean values */
typedef int bool;
enum { false, true };

#define IS_SINGLE_PRECISION (sizeof(real) == 4)
#define UDMI_ERROR_VALUE -100000


/*  helpful functions */

int countLinesOfFile(char filename[], int *line_count);

int reorderRealArr(
                     real *distribute_arr, 
                     int size_distribute_arr, 
                     int *receive_from_distribute_arr_idx_arr,
                     real *receive_arr, 
                     int size_receive_arr
                  );

int reorderRealND_ND_Arr(
                           real (*distribute_ND_ND_arr)[ND_ND], 
                           int size_distribute_arr, 
                           int *receive_from_distribute_arr_idx_arr,
                           real (*receive_ND_ND_arr)[ND_ND], 
                           int size_receive_arr
                        );

int writeRealArrToFile( 
                        char filename[],
                        real *arr,
                        int size_arr
                      );

int writeIntegerArrToFile(
                            char filename[],
                            int *arr,
                            int size_arr
                         );

#endif
