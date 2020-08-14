/*
C-functions for reading the ANSYS APDL script coupling files.

License (MIT):

Copyright (c) 2016-2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include "vof_pc_read_ansys.h"

int readElemValueVecFromAnsysOut(
                                    char filename[],
                                    real (**e_vec_prop)[ND_ND],
                                    int no_e
                                )
{
int state = _STATE_OK;

#if RP_HOST
FILE *fp;
int i = 0;
*e_vec_prop = NULL;
*e_vec_prop = (real (*)[ND_ND]) calloc(ND_ND * no_e, sizeof(real));

if(*e_vec_prop == NULL)
{
    Message("Error (readElemValueVecFromAnsysOut()):  Memory "
            "allocation error!Not enough Memory? \n");
    state = _STATE_ERROR;
}
else
{
    if ((fp = fopen(filename, "r")) == NULL)
    {
        Message("Error (readElemValueVecFromAnsysOut()): Unable to "
                "open %s, create Ansys output first!\n", filename);
        state = _STATE_ERROR;
    }
    else
    {
        Message("Reading %i vector properties from file %s\n", no_e, filename);

        #if RP_3D
        if (IS_SINGLE_PRECISION) 
        {
            while(fscanf(fp, " %e %e %e",    &(*e_vec_prop)[i][0], 
                                            &(*e_vec_prop)[i][1],
                                            &(*e_vec_prop)[i][2]
                        ) == 3)
            {
                ++i;

                if(i > no_e)
                {
                    Message("Warning (readElemValueVecFromAnsysOut()): "
                            "Reading %s, file has more lines than there is space "
                            "in coupling arrays", filename);
                    state = _STATE_ERROR;
                    break;
                }
            }
        } 
        else 
        {
            while(fscanf(fp, " %lE %lE %lE",&(*e_vec_prop)[i][0], 
                                            &(*e_vec_prop)[i][1],
                                            &(*e_vec_prop)[i][2]
                        ) == 3)
            {
                ++i;

                if(i > no_e)
                {
                    Message("Warning (readElemValueVecFromAnsysOut()): "
                            "Reading %s, file has more lines than there is space "
                            "in coupling arrays", filename);
                    break;
                }
            }
        }
        #else
        if (IS_SINGLE_PRECISION) 
        {
            while(fscanf(fp, " %e %e",   &(*e_vec_prop)[i][0], 
                                        &(*e_vec_prop)[i][1]
                ) == 2)
            {
                ++i;

                if(i > no_e)
                {
                    Message("Warning (readElemValueVecFromAnsysOut()): "
                            "Reading %s, file has more lines than there is space "
                            "in coupling arrays", filename);
                    break;
                }
            }
        } 
        else 
        {
            while(fscanf(fp, " %lE %lE", &(*e_vec_prop)[i][0], 
                                        &(*e_vec_prop)[i][1]
                        ) == 2)
            {
                ++i;

                if(i > no_e)
                {
                    Message("Warning (readElemValueVecFromAnsysOut()): "
                            "Reading %s, file has more lines than there is space "
                            "in coupling arrays", filename);
                    break;
                }
            }
        }
        #endif
        fclose(fp);
    }
}
#endif

return state;
}


int readElemValueAndVolumeFromAnsysOut(
                                            char filename[],
                                            real **e_prop,
                                            real **e_vol,
                                            int no_e
                                        )
{
int state = _STATE_OK;
#if RP_HOST
FILE *fp;
int i = 0;
*e_prop = NULL;
(*e_prop) = (real *) calloc(no_e, sizeof(real));
*e_vol = NULL;
(*e_vol) = (real *) calloc(no_e, sizeof(real));

if( *e_prop == NULL || *e_vol == NULL)
{
    Message("Error (readElemValueAndVolumeFromAnsysOut()):  Memory "
            "allocation error!Not enough Memory? \n");
    state = _STATE_ERROR;
}
else
{
    if ((fp = fopen(filename, "r")) == NULL)
    {
        Message("Error (readElemValueAndVolumeFromAnsysOut()): Unable to "
                "open %s, create Ansys output first!\n", filename);
        state = _STATE_ERROR;
    }
    else
    {
        Message("Reading %i volume properties from file %s\n",no_e, filename);

        if (IS_SINGLE_PRECISION) 
        {
            while(fscanf(fp, " %e %e", &(*e_prop)[i],&(*e_vol)[i]) == 2)
            {
                ++i;

                if(i > no_e)
                {
                    Message("Warning (readElemValueAndVolumeFromAnsysOut()): "
                            "Reading %s, file has more lines than there is space "
                            "in coupling arrays", filename);
                    state = _STATE_ERROR;
                    break;
                }
            }
        } 
        else 
        {
            while(fscanf(fp, " %lE %lE", &(*e_prop)[i],&(*e_vol)[i]) == 2)
            {
                ++i;

                if(i > no_e)
                {
                    Message("Warning (readElemValueAndVolumeFromAnsysOut()): "
                            "Reading %s, file has more lines than there is space "
                            "in coupling arrays", filename);
                    break;
                }
            }
        }
        fclose(fp);
    }
}
#endif
return state;
}



int readCoordinatesFromAnsysOut(
                                char filename[],
                                real (**coord_arr_ansys)[ND_ND],
                                int *size_coord_arr_ansys
                            )
{
    /* Where is the error? */
    FILE *fp = NULL;
    int state = _STATE_OK;
    real dummy;
    int i = 0;
    int ansys_elements = 0;

    state = countLinesOfFile(filename, &ansys_elements);
    (*coord_arr_ansys) = NULL;

    if (state != _STATE_ERROR)
    {
        (*coord_arr_ansys) = (real (*)[ND_ND]) calloc(ND_ND * ansys_elements, sizeof(real));
    }
    else
    {
        Message("Error (readCoordinatesAnsysAllOut()):  Error in "
                "getAnsysElementCountFromFileLength().\n");
    }

    if((*coord_arr_ansys) == NULL)
    {
        Message("Error (readCoordinatesAnsysAllOut()): Memory "
                "allocation error! Not enough Memory? \n");
        state = _STATE_ERROR;
    }
    else if (state != _STATE_ERROR)
    {
        if ((fp = fopen(filename, "r")) == NULL)
        {
            Message("Error (readCoordinatesAnsysAllOut()): Unable to "
                "open %s, create Ansys output first!\n", filename);
            state = _STATE_ERROR;
        }
        else
        {
            Message("Info (readCoordinatesAnsysAllOut()):Reading element "
                    "coordinates from ANSYS table %s...\n", filename);

            #if RP_3D
            if (IS_SINGLE_PRECISION)
            {
                while(fscanf(fp, " %e %e %e", 
                &(*coord_arr_ansys)[i][0], &(*coord_arr_ansys)[i][1], &(*coord_arr_ansys)[i][2]) == 3)
                {
                    ++i;

                    if(i > ansys_elements)
                    {
                        Message("Warning (readCoordinatesAnsysAllOut()): "
                        "Problem while reading file %s," 
                        "count of read lines not equal to assumed element "
                        "count please check file or function "
                        "getAnsysElementCountFromFileLength()\n", filename);
                        state = _STATE_ERROR;
                        break;
                    }
                }
                
            } 
            else 
            {
                while(fscanf(fp, " %lE %lE %lE",
                &(*coord_arr_ansys)[i][0], &(*coord_arr_ansys)[i][1], &(*coord_arr_ansys)[i][2]) == 3)
                {
                    ++i;

                    if(i > ansys_elements)
                    {
                        Message("Warning (readCoordinatesAnsysAllOut()): "
                        "Problem while reading file %s," 
                        "count of read lines not equal to assumed element "
                        "count please check file or function "
                        "getAnsysElementCountFromFileLength()\n", filename);
                        state = _STATE_ERROR;
                        break;
                    }
                }
            }
            #else 
            if (IS_SINGLE_PRECISION) 
            {
                while(fscanf(fp, " %e %e", 
                    &(*coord_arr_ansys)[i][0], &(*coord_arr_ansys)[i][1]) == 2)
                {
                    ++i;

                    if(i > ansys_elements)
                    {
                        Message("Warning (readCoordinatesAnsysAllOut()): "
                        "Problem while reading file %s," 
                        "count of read lines not equal to assumed element "
                        "count please check file or function "
                        "getAnsysElementCountFromFileLength()\n", filename);
                        state = _STATE_ERROR;
                        break;
                    }
                }
            } 
            else 
            {
                while(fscanf(fp, " %lE %lE", 
                    &(*coord_arr_ansys)[i][0], &(*coord_arr_ansys)[i][1]) == 2)
                {
                    ++i;

                    if(i > ansys_elements)
                    {
                        Message("Warning (readCoordinatesAnsysAllOut()): "
                        "Problem while reading file %s," 
                        "count of read lines not equal to assumed element "
                        "count please check file or function "
                        "getAnsysElementCountFromFileLength()\n", filename);
                        state = _STATE_ERROR;
                        break;
                    }
                }
            }
            #endif
           fclose(fp);
        }
    }
    (*size_coord_arr_ansys) = ansys_elements;
    

    return state;
}
