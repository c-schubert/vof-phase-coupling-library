/*
C-functions containing the implementation of NN mapping algorithm and exchange 
from mapped array back to ANSYS Fluent nodes and to UDMI values.

License (MIT):

Copyright (c) 2016-2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "vof_pc_nn_mapping.h"


int safeDistributeMappedArrayToNodesInCellZone(
                                                real *mapped_arr,
                                                int size_mapped_arr,
                                                int noUDMI,
                                                int *fluent_compute_node_id_arr_full,
                                                int *fluent_cell_id_arr_full,
                                                int *cells_per_node,
                                                int cell_zone_id
                                                )
{
/* Checks if distributing of mapped_arr see distributeMappedArrayFromHostToUDMI()
   would be in accordance with the mapping 
 */
    int pe;
    int state = _STATE_OK;
    int i = 0;
    real *node_n_mapped_arr = NULL;
    int *node_n_compute_node_id_arr = NULL;
    int *node_n_cell_id_arr = NULL;
    int cells_in_node_n = 0;

    #if !RP_HOST
    Thread *t;
    Domain *domain = Get_Domain(1); // mixture domain
    t = Lookup_Thread(domain, cell_zone_id);
    #endif

    #if RP_HOST
    int ii;
    int i_start = 0;

    compute_node_loop (pe) 
    { 
        if (pe < compute_node_count)
        {
            cells_in_node_n = cells_per_node[pe];
            #if VOF_PC_DEBUG
            Message("Info safeDistributeMappedArrayToNodes(): Compute "
                    "node count %i is %i.\n", pe, cells_in_node_n);
            #endif
        }
        else
        {
            cells_in_node_n = 0;
        }

        if(
            cells_in_node_n > 0 && 
            pe < compute_node_count && 
            state != _STATE_ERROR
            )
        {
            node_n_mapped_arr = (real *) calloc(cells_in_node_n, sizeof(real));
            node_n_compute_node_id_arr = (int *) calloc(cells_in_node_n, sizeof(int));
            node_n_cell_id_arr = (int *) calloc(cells_in_node_n, sizeof(int));

            if(
                node_n_compute_node_id_arr == NULL ||
                node_n_cell_id_arr == NULL ||
                node_n_mapped_arr == NULL
                )
            {
                Message("Error safeDistributeMappedArrayToNodes(): Memory "
                "allocation failed?!\n");
                state = _STATE_ERROR;
            }
            else
            {
                /* Bound check i - outside*/
                if(i_start >= 0 && (i_start + cells_in_node_n) <= size_mapped_arr)
                {
                    for(i = i_start; i < (i_start + cells_in_node_n); ++i)
                    {
                        ii = (i-i_start);
                        if(ii < cells_in_node_n && ii >= 0)
                        {
                            node_n_mapped_arr[ii] = mapped_arr[i];
                            node_n_compute_node_id_arr[ii] = fluent_compute_node_id_arr_full[i];
                            node_n_cell_id_arr[ii] = fluent_cell_id_arr_full[i];
                        }
                        else
                        {
                            Message("Error safeDistributeMappedArrayToNodes(): Index "
                            "out of bounds 1!\n");
                            state = _STATE_ERROR;
                        }
                    }
                    
                    i_start += cells_in_node_n;
                }
                else
                {
                    Message("Error safeDistributeMappedArrayToNodes(): Index "
                            "out of bounds 2!\n");
                    state = _STATE_ERROR;
                }
            }
        }
        else
        {
            if(pe > compute_node_count)
            {
                Message("Error safeDistributeMappedArrayToNodes(): "
                        "Wrong compute node count...\n");
                state = _STATE_ERROR;
            }
            if(cells_in_node_n <= 0)
            {
                Message("Warning safeDistributeMappedArrayToNodes(): "
                        "No cells in node %i correct?\n", pe);
                state = _STATE_WARNING;
            }
        }

        PRF_CSEND_INT(node_zero, &state, 1, node_host);

        if(state != _STATE_ERROR)
        {
            #if VOF_PC_DEBUG
            Message("host sent to node zero for %i!\n", pe);
            #endif

            PRF_CSEND_INT(node_zero, &cells_in_node_n, 1, node_host);
            PRF_CSEND_REAL(node_zero, node_n_mapped_arr, cells_in_node_n, node_host);
            PRF_CSEND_INT(node_zero, node_n_compute_node_id_arr, cells_in_node_n, node_host);
            PRF_CSEND_INT(node_zero, node_n_cell_id_arr, cells_in_node_n, node_host);
        }

        if(node_n_mapped_arr != NULL)
        {
            free(node_n_mapped_arr);
            node_n_mapped_arr = NULL;
        }
        
        if(node_n_compute_node_id_arr != NULL)
        {
            free(node_n_compute_node_id_arr);
            node_n_compute_node_id_arr = NULL;
        }
        
        if(node_n_cell_id_arr != NULL)
        {
            free(node_n_cell_id_arr);
            node_n_cell_id_arr = NULL;
        }
    }

    if(state == _STATE_ERROR)
    {
        Message("Error safeDistributeMappedArrayToNodes(): "
                "In mapping from host!\n");
    }
    
    #endif


    #if !RP_HOST
    if (I_AM_NODE_ZERO_P)
    {
        compute_node_loop (pe) 
        { 
            PRF_CRECV_INT(node_host, &state, 1, node_host);

            if(state != _STATE_ERROR)
            {   
                #if VOF_PC_DEBUG
                Message("node zero reveice for %i!\n", pe);
                #endif

                PRF_CRECV_INT(node_host, &cells_in_node_n, 1, node_host);

                node_n_mapped_arr = (real *) calloc(cells_in_node_n, sizeof(real));
                node_n_compute_node_id_arr = (int *) calloc(cells_in_node_n, sizeof(int));
                node_n_cell_id_arr = (int *) calloc(cells_in_node_n, sizeof(int));

                PRF_CRECV_REAL(node_host, node_n_mapped_arr, cells_in_node_n, node_host);
                PRF_CRECV_INT(node_host, node_n_compute_node_id_arr, cells_in_node_n, node_host);
                PRF_CRECV_INT(node_host, node_n_cell_id_arr, cells_in_node_n, node_host);

                if (pe == myid)
                {
                    state = safeNodeRealArrayToCUDMI(
                                                        t,
                                                        noUDMI,
                                                        node_n_mapped_arr,
                                                        node_n_compute_node_id_arr,
                                                        node_n_cell_id_arr,
                                                        cells_in_node_n
                                                    );
                }

                if(pe != myid)
                {
                    #if VOF_PC_DEBUG
                    Message("node zero sent to %i!\n", pe);
                    #endif
                    PRF_CSEND_INT(pe, &state, 1, myid);

                    PRF_CSEND_INT(pe, &cells_in_node_n, 1, myid);
                    PRF_CSEND_REAL(pe, node_n_mapped_arr, cells_in_node_n, myid);
                    PRF_CSEND_INT(pe, node_n_compute_node_id_arr, cells_in_node_n, myid);
                    PRF_CSEND_INT(pe, node_n_cell_id_arr, cells_in_node_n, myid);
                }

                if(node_n_mapped_arr != NULL)
                {
                    free(node_n_mapped_arr);
                    node_n_mapped_arr = NULL;
                }
            
                if(node_n_compute_node_id_arr != NULL)
                {
                    free(node_n_compute_node_id_arr);
                    node_n_compute_node_id_arr = NULL;
                }
                
                if(node_n_cell_id_arr != NULL)
                {
                    free(node_n_cell_id_arr);
                    node_n_cell_id_arr = NULL;
                }
            }
            else
            {
                
                if(pe != myid)
                {
                    PRF_CSEND_INT(pe, &state, 1, myid);
                }
            }
            
        }
    }
    else
    {
        PRF_CRECV_INT(node_zero, &state, 1, node_zero);

        if(state != _STATE_ERROR)
        {
            #if VOF_PC_DEBUG
            Message("node %i reveice from node zero!\n", pe);
            #endif

            PRF_CRECV_INT(node_zero, &cells_in_node_n, 1, node_zero);

            node_n_mapped_arr = (real *) calloc(cells_in_node_n, sizeof(real));
            node_n_compute_node_id_arr = (int *) calloc(cells_in_node_n, sizeof(int));
            node_n_cell_id_arr = (int *) calloc(cells_in_node_n, sizeof(int));

            PRF_CRECV_REAL(node_zero, node_n_mapped_arr, cells_in_node_n, node_zero);
            PRF_CRECV_INT(node_zero, node_n_compute_node_id_arr, cells_in_node_n, node_zero);
            PRF_CRECV_INT(node_zero, node_n_cell_id_arr, cells_in_node_n, node_zero);

            state = safeNodeRealArrayToCUDMI(
                                                t, 
                                                noUDMI,
                                                node_n_mapped_arr,
                                                node_n_compute_node_id_arr,
                                                node_n_cell_id_arr,
                                                cells_in_node_n
                                            );

            if(node_n_mapped_arr != NULL)
            {
                free(node_n_mapped_arr);
                node_n_mapped_arr = NULL;
            }
        
            if(node_n_compute_node_id_arr != NULL)
            {
                free(node_n_compute_node_id_arr);
                node_n_compute_node_id_arr = NULL;
            }
            
            if(node_n_cell_id_arr != NULL)
            {
                free(node_n_cell_id_arr);
                node_n_cell_id_arr = NULL;
            }
        }
    }
    #endif

    return state;
}


int safeNodeRealArrayToCUDMI(
                                Thread *t, 
                                int noUDMI,
                                real *node_n_mapped_arr,
                                int *node_n_compute_node_id_arr,
                                int *node_n_cell_id_arr,
                                int cells_in_node_n
                            )
{
int state = _STATE_OK;

#if RP_NODE
    cell_t c; 
    int cells_in_thread = 0;
    int i = 0;
    cells_in_thread = THREAD_N_ELEMENTS_INT(t);

if(
    N_UDM > noUDMI &&
    cells_in_thread == cells_in_node_n &&
    node_n_mapped_arr != NULL &&
    node_n_compute_node_id_arr != NULL &&
    node_n_cell_id_arr != NULL
)
{
    i = 0;
    begin_c_loop_int(c, t)
    {
        if(node_n_cell_id_arr[i] != (int) c)
        {
            Message("Error safeNodeArrayToCUDMI(): "
                    "Wrong cell mapping during "
                    "redistribution to nodes!\n");
            state = _STATE_ERROR;
        }

        if(
            node_n_compute_node_id_arr[i] != myid 
            && state != _STATE_ERROR
            )
        {
            Message("Error safeNodeArrayToCUDMI(): "
                    "Wrong comute node mapping during "
                    "redistribution to nodes!\n");
            state = _STATE_ERROR;
        }

        if(state != _STATE_ERROR)
        {
            if (i < cells_in_thread)
                C_UDMI(c,t,noUDMI) = node_n_mapped_arr[i];
        }
        else
        {
            C_UDMI(c,t,noUDMI) = UDMI_ERROR_VALUE;
        }

        ++i;
    }
    end_c_loop_int(c, t)
}
else
{   

    if(cells_in_thread != cells_in_node_n)
    {
            Message("Error safeNodeArrayToCUDMI(): "
            "Missmatch of cell count in "
            "distributeMappedArrayFromHostToUDMI()!\n");
            state = _STATE_ERROR;
    }

    if(
        node_n_mapped_arr == NULL ||
        node_n_compute_node_id_arr == NULL ||
        node_n_cell_id_arr == NULL
      )
    {
        Message("Error safeNodeArrayToCUDMI(): "
        "One or multiple input arrays pointing to NULL!\n");
        state = _STATE_ERROR;

    }

    if( N_UDM <= noUDMI)
    {
        Message("Error safeNodeArrayToCUDMI(): "
            "Not enough UDMI's for UDMI index!\n", noUDMI);
            state = _STATE_ERROR;
    }
}
#endif

return state;
}



void nearestNeighborMatching(
                                real (*coord_arr_1)[ND_ND],
                                int size_arr_1,
                                real (*coord_arr_2)[ND_ND],
                                int size_arr_2, 
                                int **mapping_arr1_to_arr2, 
                                int **mapping_arr2_to_arr1
                            )
{
    int i = 0;
    int j = 0;

    real x[ND_ND];
    real squared_dist = 0;
    real min_arr1_to_arr2 = 0;
    real *min_arr2_to_arr1= NULL;

    *mapping_arr1_to_arr2 = NULL;
    *mapping_arr2_to_arr1 = NULL;

    min_arr2_to_arr1 = (real *) calloc(size_arr_2, sizeof(real));

    *mapping_arr1_to_arr2 = (int *) calloc(size_arr_1, sizeof(int));
    *mapping_arr2_to_arr1 = (int *) calloc(size_arr_2, sizeof(int));
    
    if( min_arr2_to_arr1 == NULL || 
        *mapping_arr1_to_arr2 == NULL ||
        *mapping_arr2_to_arr1 == NULL
     )
    {
        Message("Error at allocation in nearestNeighborMatching"
                " , not enough Memory?\n");
    }
    else
    {
        Message("Start NN Mapping ...\n");		
        /*Init distances for min_arr2_to_arr1*/
        for (j = 0; j < size_arr_2; ++j)
        {
            NV_VV(x, =, coord_arr_1[0], - , coord_arr_2[j]);
            min_arr2_to_arr1[j] = NV_MAG2(x);

            (*mapping_arr2_to_arr1)[j] = 0;
        }


        for (i = 0; i < size_arr_1; ++i)
        {
            /* j = 0 */
            j = 0;
            {
                NV_VV(x, =, coord_arr_1[i], - , coord_arr_2[j]);
                squared_dist =  NV_MAG2(x);
                min_arr1_to_arr2 = squared_dist;
                (*mapping_arr1_to_arr2)[i] = j;

                if (squared_dist < min_arr2_to_arr1[j])
                {
                    (*mapping_arr2_to_arr1)[j] = i;
                    min_arr2_to_arr1[j] = squared_dist;
                }
            }
        
            for (j = 1; j < size_arr_2; ++j)
            {
                NV_VV(x, =, coord_arr_1[i], - , coord_arr_2[j]);
                squared_dist = NV_MAG2(x);

                if(squared_dist < min_arr1_to_arr2)
                {
                    (*mapping_arr1_to_arr2)[i] = j;
                    min_arr1_to_arr2 = squared_dist;
                }

                if (squared_dist < min_arr2_to_arr1[j])
                {
                    (*mapping_arr2_to_arr1)[j] = i;
                    min_arr2_to_arr1[j] = squared_dist;
                }
            }
            
            if( (i+1)%((int) size_arr_1/5) == 0 )
            {
                Message("NN Mapping running...\n");
            }

        }
        Message("NN Mapping finished!\n");
    }
}


int getAnsysElementCountFromFileLength(char filename[])
{
    int lines = 0;
    int state = _STATE_OK;

    state = countLinesOfFile(filename, &lines);

    if(state != _STATE_OK)
    {
        Message("Error (getAnsysElementCountFromFileLength()): Posible error in "
               "countLinesOfFile()!\n");
    }
    else if(lines < 1)
    {
        Message("Warning (getAnsysElementCountFromFileLength()): Lines smaller 1 "
               "should not happen!\n");
        state = _STATE_WARNING;
    }

    return state;
}


int debugWriteMappings(
                            int *mapping_arr1_to_arr2, 
                            int size_arr_1,
                            int *mapping_arr2_to_arr1,
                            int size_arr_2, 
                            char filename_arr1_to_arr2_mapping[],
                            char filename_arr2_to_arr1_mapping[]
                        )
{

    int state = 0;

    state = writeIntegerArrToFile(
                                    filename_arr1_to_arr2_mapping,
                                    mapping_arr1_to_arr2,
                                    size_arr_1
                                 );

    if(state != _STATE_OK)
    {
        Message("Warning (debugWriteMappings()): Possible error in "
                "writeIntegerArrToFile()-1\n");
    }
    state  = writeIntegerArrToFile( 
                                    filename_arr2_to_arr1_mapping,
                                    mapping_arr2_to_arr1,
                                    size_arr_2
                                  );
    if(state != _STATE_OK)
    {
        Message("Warning (debugWriteMappings()): Possible error in "
                "writeIntegerArrToFile()-2\n");
    }
    else
    {
        Message("Info (debugWriteMappings()): Debug NN Mappings written to %s "
                "and %s.\n"
                , filename_arr1_to_arr2_mapping, filename_arr2_to_arr1_mapping);
    }

    return state;
}


int writePropertyArrayToFileMapped(
                                    char filename[],
                                    real *property_arr_full,
                                    int size_property_arr_full,
                                    int *mapping_arr,
                                    int size_mapping_arr
                                )
{
    real *mapped_property_arr = NULL;
    int state = _STATE_OK;

    mapped_property_arr = (real *) calloc(size_mapping_arr, sizeof(real));

    if (mapped_property_arr == NULL)
    {
        Message("Error (writePropertyArrayToFileMapped()): Memory allocation in "
                "writeVofToAnsys()!\n");
        state = _STATE_ERROR;
    }
    else
    {
        state = reorderRealArr(
                    property_arr_full, 
                    size_property_arr_full, 
                    mapping_arr,
                    mapped_property_arr, 
                    size_mapping_arr
                );

        if( state != _STATE_ERROR)
        {
            state = writeRealArrToFile(
                                            filename,
                                            mapped_property_arr,
                                            size_mapping_arr
                                      );
            if(state == _STATE_ERROR)
            {
                Message("Error (writePropertyArrayToFileMapped()): Error in "
                        "writeRealArrToFile()!\n");
            }
        }
        else
        {
            Message("Error (writePropertyArrayToFileMapped()): Error in "
                    "reorderArr()!\n");
        }
    }
    
    return state;
}
