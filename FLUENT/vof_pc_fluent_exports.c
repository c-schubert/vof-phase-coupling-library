/*
Functions to write fluent fields to files for debug purposes

License (MIT):

Copyright (c) 2016-2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include "vof_pc_fluent_exports.h" 

DEFINE_ON_DEMAND(Debug_Export_Fluent_Fields)
{
    int i = 0;
    char filename[] = _FLUENT_ALLOUT_DAT_;
    real (*coord_arr_full)[ND_ND];
    real *vof_arr_full;
    int *cell_id_arr_full;
    int *compute_node_id_arr_full;
    int arr_full_size = -1;

    hostGetCouplingFieldsFromNodesinCellZone(
                                                &coord_arr_full,
                                                &vof_arr_full, 
                                                &cell_id_arr_full,
                                                &compute_node_id_arr_full,
                                                &arr_full_size,
                                                FLUID_ID
                                              );
    
    hostWriteDebugField (
                            filename,
                            coord_arr_full,
                            vof_arr_full, 
                            cell_id_arr_full,
                            compute_node_id_arr_full,
                            arr_full_size
                        );
    #if RP_HOST
        Message("Export Done!\n");
    #endif
}


int hostWriteDebugField(
                            char filename[],
                            real (*coord_arr_full)[ND_ND],
                            real *vof_arr_full, 
                            int *cell_id_arr_full, 
                            int *compute_node_id_arr_full,
                            int arr_full_size
                        )
{
/*
    Write a file with the follwing columns:

    cell_id compute_node_id volfrac x-coord y-coord (z-coord only if 3D)
*/
int state = _STATE_OK;

#if RP_HOST
FILE *fp = NULL;
int i = 0;

if ((fp = fopen(filename, "w")) == NULL)
{
    Message("Error (hostWriteDebugField()): Unable to open %s for writing!\n"
                , filename);
    state = _STATE_ERROR;
}
else
{
    //fprintf(fp, "cell_ids node_ids volfrac x y z\n");

    if(arr_full_size > 0)
    {
        if ( vof_arr_full == NULL || cell_id_arr_full == NULL 
            || coord_arr_full == NULL || compute_node_id_arr_full == NULL)
        {
            Message("Error (hostWriteDebugField()): Can not export empty arrays!\n"
                        , filename);
            state = _STATE_ERROR;
        }
        else
        {
            #if RP_3D
            for (i = 0; i < arr_full_size; ++i)
            {
                if(IS_SINGLE_PRECISION)
                    fprintf(fp, "%i %i %f %f %f %f\n", cell_id_arr_full[i],
                             compute_node_id_arr_full[i], vof_arr_full[i], 
                             coord_arr_full[i][0], coord_arr_full[i][1], 
                             coord_arr_full[i][2]);
                else
                    fprintf(fp, "%i %i %lf %lf %lf %lf\n", cell_id_arr_full[i], 
                            compute_node_id_arr_full[i], vof_arr_full[i], 
                            coord_arr_full[i][0], coord_arr_full[i][1], 
                            coord_arr_full[i][2]);

            }
            #else
            for (i = 0; i < arr_full_size; ++i)
            {
                if(IS_SINGLE_PRECISION)
                    fprintf(fp, "%i %i %f %f %f\n", cell_id_arr_full[i], 
                             compute_node_id_arr_full[i], vof_arr_full[i], 
                             coord_arr_full[i][0], coord_arr_full[i][1]);
                else
                    fprintf(fp, "%i %i %lf %lf %lf\n", cell_id_arr_full[i], 
                             compute_node_id_arr_full[i],vof_arr_full[i], 
                             coord_arr_full[i][0], coord_arr_full[i][1]);
            }
            #endif /* RP_3D */

            free(coord_arr_full);
            free(vof_arr_full);
            free(cell_id_arr_full);
            free(compute_node_id_arr_full);
        }
    }
    else
    {
        Message("Warning (hostWriteDebugField()): Array sizes are zero!\n");
        state = _STATE_ERROR;
    }
    fclose(fp);
}
#endif

return state;
}
