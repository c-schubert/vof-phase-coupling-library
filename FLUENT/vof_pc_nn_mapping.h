/*
C-functions containing the implementation of NN mapping algorithm and exchange 
from mapped array back to ANSYS Fluent nodes and to UDMI values.

License (MIT):

Copyright (c) 2016-2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#ifndef VOF_PC_NN_MAPPING_H
#include "vof_pc_main.h"
#include "vof_pc_fluent_get_fields.h" 
#define VOF_PC_NN_MAPPING_H


int debugWriteMappings(
                          int *mapping_arr1_to_arr2, 
                          int size_arr_1,
                          int *mapping_arr2_to_arr1,
                          int size_arr_2, 
                          char filename_arr1_to_arr2_mapping[],
                          char filename_arr2_to_arr1_mapping[]
                      );

int getAnsysElementCountFromFileLength(
                                        char filename[]
                                      );


void nearestNeighborMatching(
                                real (*coord_arr_1)[ND_ND],
                                int size_arr_1,
                                real (*coord_arr_2)[ND_ND],
                                int size_arr_2, 
                                int **mappings_arr1_to_arr2, 
                                real **weights_arr1_to_arr2, 
                                int **mappings_arr2_to_arr1,
                                real **weights_arr2_to_arr1
                            );

int safeDistributeMappedArrayToNodesInCellZone(
                                                real *mapped_arr,
                                                int size_mapped_arr,
                                                int noUDMI,
                                                int *fluent_compute_node_id_arr_full,
                                                int *fluent_cell_id_arr_full,
                                                int *cells_per_node,
                                                int cell_zone_id
                                                );

int safeNodeRealArrayToCUDMI(
                                Thread *t, 
                                int noUDMI,
                                real *node_n_mapped_arr,
                                int *node_n_compute_node_id_arr,
                                int *node_n_cell_id_arr,
                                int cells_in_node_n
                            );

int writePropertyArrayToFileMapped(
                                    char filename[],
                                    real *property_arr_full,
                                    int size_property_arr_full,
                                    int *mapping_arr,
                                    int size_mapping_arr
                                );

#endif
