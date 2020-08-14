/*
Functions to get fields and node distribution information from nodes to host

License (MIT):

Copyright (c) 2016-2019 Christian Schubert 

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef VOF_PC_FLUENT_GET_FIELDS_H
#include "vof_pc_main.h"
#include "vof_pc_case.h"

int hostGetCellCountPerNodeInCellZone( 
                                      int **cells_per_node,
                                      int *cell_sum_over_nodes,
                                      int cellZoneID
                                      );

real get_c_vof(cell_t c, Thread *t);
real get_coord(cell_t c,Thread *t,int dim);
real get_x_coord(cell_t c,Thread *t);
real get_y_coord(cell_t c,Thread *t);
real get_z_coord(cell_t c,Thread *t);

void hostGetCellCoordsFromNodesInCellZone(
                                          real (**coord_arr_full)[ND_ND],
                                          int length_arrs_full,
                                          int cell_zone
                                          );

void hostGetOrderingArraysFromNodesInCellZone(
                                                int **cid_arr_full, 
                                                int **myid_arr_full,
                                                int *length_arrs_full,
                                                int cell_zone
                                              );
                                              
void hostGetOrderedFieldValueArrayFromNodesInCellZone(
                                                      real **val_arr_full,
                                                      real (*C_VAL_WRAPPER_FUN)(cell_t, Thread*),
                                                      int *length_arrs_full, 
                                                      int cell_zone
                                                    );

#define VOF_PC_FLUENT_GET_FIELDS_H

#endif
