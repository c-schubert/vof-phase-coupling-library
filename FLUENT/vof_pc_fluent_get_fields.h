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

int hostGetCellCountPerNodeInCellZone( 
                                      int **cells_per_node,
                                      int *cell_sum_over_nodes,
                                      int cellZoneID
                                      );

void hostGetCouplingFieldsFromNodesinCellZone(
                                                real (**coord_arr_full)[ND_ND],
                                                real **vof_arr_full, 
                                                int **cell_id_arr_full, 
                                                int **compute_node_id_arr_full,
                                                int *arr_full_size,
                                                int cell_zone
                                             );

#define VOF_PC_FLUENT_GET_FIELDS_H

#endif
