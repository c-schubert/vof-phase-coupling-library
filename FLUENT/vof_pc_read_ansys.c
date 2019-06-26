/*
C-functions for reading the ANSYS APDL script coupling files.

License (MIT):

Copyright (c) 2016-2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include "vof_pc_read_ansys.h"

int readJouleHeatAndLorentzForcesFromAnsysOut(
                                                char filename[],
                                                real (**lorentz_force)[ND_ND],
                                                real **joule_heat,
                                                int length_arr
                                             )
{
FILE *fp;
int state = _STATE_OK;
int i = 0;
*lorentz_force = NULL;
*joule_heat = NULL;
*lorentz_force = (real (*)[ND_ND]) calloc(ND_ND * length_arr, sizeof(real));
*joule_heat = (real *) calloc(length_arr, sizeof(real));

if(*lorentz_force == NULL || *joule_heat == NULL)
{
    Message("Error (readJouleHeatAndLorentzForcesFromAnsysOut()):  Memory "
            "allocation error!Not enough Memory? \n");
    state = _STATE_ERROR;
}
else
{
    if ((fp = fopen(filename, "r")) == NULL)
    {
        Message("Error (readJouleHeatAndLorentzForcesFromAnsysOut()): Unable to "
                "open %s, create Ansys output first!\n", filename);
        state = _STATE_ERROR;
    }
    else
    {
        #if RP_3D
        if (IS_SINGLE_PRECISION) 
        {
            while(fscanf(fp, "%f %f %f %f", &(*joule_heat)[i], 
                                            &(*lorentz_force)[i][0], 
                                            &(*lorentz_force)[i][1],
                                            &(*lorentz_force)[i][2]
                        ) == 4)
            {
                ++i;

                if(i > length_arr)
                {
                    Message("Warning (readJouleHeatAndLorentzForcesFromAnsysOut()): "
                            "Reading %s, file has more lines than there is space "
                            "in coupling arrays", filename);
                    state = _STATE_ERROR;
                    break;
                }
            }
        } 
        else 
        {
            while(fscanf(fp, "%lf %lf %lf %lf", &(*joule_heat)[i], 
                                            &(*lorentz_force)[i][0], 
                                            &(*lorentz_force)[i][1],
                                            &(*lorentz_force)[i][2]
                        ) == 4)
            {
                ++i;

                if(i > length_arr)
                {
                    Message("Warning (readJouleHeatAndLorentzForcesFromAnsysOut()): "
                            "Reading %s, file has more lines than there is space "
                            "in coupling arrays", filename);
                    break;
                }
            }
        }
        #else
        if (IS_SINGLE_PRECISION) 
        {
            while(fscanf(fp, "%f %f %f", &(*joule_heat)[i], 
                                            &(*lorentz_force)[i][0], 
                                            &(*lorentz_force)[i][1]
                        ) == 3)
            {
                ++i;

                if(i > length_arr)
                {
                    Message("Warning (readJouleHeatAndLorentzForcesFromAnsysOut()): "
                            "Reading %s, file has more lines than there is space "
                            "in coupling arrays", filename);
                    break;
                }
            }
        } 
        else 
        {
            while(fscanf(fp, "%lf %lf %lf", &(*joule_heat)[i], 
                                            &(*lorentz_force)[i][0], 
                                            &(*lorentz_force)[i][1]
                        ) == 3)
            {
                ++i;

                if(i > length_arr)
                {
                    Message("Warning (readJouleHeatAndLorentzForcesFromAnsysOut()): "
                            "Reading %s, file has more lines than there is space "
                            "in coupling arrays", filename);
                    break;
                }
            }
        }
        #endif
    }
}

return state;
}

int readCoordinatesAnsysAllOut(
                                    char filename[],
                                    real (**coord_arr_ansys)[ND_ND],
                                    int *size_coord_arr_ansys
                                )
{
    FILE *fp;
    int state = _STATE_OK;
    real dummy;
    int i = 0;
    int ansys_elements = 0;


    state = countLinesOfFile(filename, &ansys_elements);

    if (state != _STATE_ERROR)
    {
        *coord_arr_ansys = NULL;
        *coord_arr_ansys = (real (*)[ND_ND]) calloc(ND_ND * ansys_elements, sizeof(real));
    }
    else
    {
        Message("Error (readCoordinatesAnsysAllOut()):  Error in "
                "getAnsysElementCountFromFileLength().\n");
    }

    if(*coord_arr_ansys == NULL)
    {
        Message("Error (readCoordinatesAnsysAllOut()):  Memory "
                "allocation error!Not enough Memory? \n");
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
                while(fscanf(fp, "%f %f %f %f %f %f %f %f", 
                &(*coord_arr_ansys)[i][0], &(*coord_arr_ansys)[i][1], &(*coord_arr_ansys)[i][2],
                &dummy, &dummy, &dummy, &dummy, &dummy) == 8)
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
                while(fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf",
                &(*coord_arr_ansys)[i][0], &(*coord_arr_ansys)[i][1], &(*coord_arr_ansys)[i][2],
                &dummy, &dummy, &dummy, &dummy, &dummy) == 8)
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
                while(fscanf(fp, "%f %f %f %f %f %f %f", 
                    &(*coord_arr_ansys)[i][0], &(*coord_arr_ansys)[i][1], 
                    &dummy, &dummy, &dummy, &dummy, &dummy) == 7)
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
                while(fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", 
                    &(*coord_arr_ansys)[i][0], &(*coord_arr_ansys)[i][1], 
                    &dummy, &dummy, &dummy, &dummy, &dummy) == 7)
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
        *size_coord_arr_ansys = ansys_elements;
        fclose(fp);
        }
    }

    return state;
}
