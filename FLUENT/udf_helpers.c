/*
License (MIT):

Copyright (c) 2016-2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include "udf_helpers.h"


DEFINE_ON_DEMAND(f_parallelInfo)
{
    /*Prints out parallelisation info*/
#if RP_NODE 
    int computeNodeCount = compute_node_count;
    int idLastNode = node_last;
    int myID = myid;

    Message("myid: %i \n", myID);
    Message0("last node: %i \n", idLastNode);
    Message0("compute node count: %i \n", computeNodeCount);
#endif
}


int countLinesOfFile(char filename[], int *line_count)
{
/*
    returns line count of given file (filename)
*/
    FILE *fp;
    int state = _STATE_OK;
    int ch;

    if ((fp = fopen(filename, "r")) == NULL)
    {
        Message("Error (countLinesOfFile()): Unable to open %s!\n", filename);
        state = _STATE_ERROR;
    }
    else
    {
        Message("Info (countLinesOfFile()): Counting lines of file %s\n", filename);
        (*line_count) = 0;
        while ((ch = fgetc(fp)) != EOF)
        {
            ch = fgetc(fp);
            if(ch == '\n')
            {
                (*line_count)++;
            }
        }
        fclose(fp);
    }

    return state;
}


int reorderRealArr(
                real *distribute_arr, 
                int size_distribute_arr, 
                int *receive_from_distribute_arr_idx_arr,
                real *receive_arr, 
                int size_receive_arr
            )
{
/*
    maps the values of distribute_arr, for the indices in array 
    receive_from_distribute_arr_idx_arr, to receive_arr.

    receive_arr must be allocated before!
*/
int i = 0;
int j = 0;
int state = _STATE_OK;

for(i = 1; i < size_receive_arr; ++i)
{
    j = receive_from_distribute_arr_idx_arr[i];

    if(j < size_distribute_arr)
    {
        receive_arr[i] = distribute_arr[j];
    }
    else
    {
        Message("Error (reorderRealArr ()): receive_from_distribute_arr_idx_arr[%i] "
                "greater than size_distribute_arr! \n",i);
        state = _STATE_ERROR;
    }
}

    return state;
}


int reorderRealND_ND_Arr(
                        real (*distribute_ND_ND_arr)[ND_ND], 
                        int size_distribute_arr, 
                        int *receive_from_distribute_arr_idx_arr,
                        real (*receive_ND_ND_arr)[ND_ND], 
                        int size_receive_arr
                     )
{
/*
    maps ND_ND dims of distribute_ND_ND_arr, for the indices in array 
    receive_from_distribute_arr_idx_arr, to ND_ND dims to receive_ND_ND_arr.

    receive_ND_ND_arr must be allocated before!
*/
int i = 0;
int j = 0;
int k = 0;
int state = _STATE_OK;

for(i = 1; i < size_receive_arr; ++i)
{
    j = receive_from_distribute_arr_idx_arr[i];

    if(j < size_distribute_arr)
    {
        for(k = 0; k < ND_ND; ++k)
        {
            receive_ND_ND_arr[i][k] = distribute_ND_ND_arr[j][k];
        }
    }
    else
    {
        Message("Error (reorderND_ND_Arr()): receive_from_distribute_arr_idx_arr[%i] "
                "greater than size_distribute_arr! \n",i);
        state = _STATE_ERROR;
    }
}

return state;
}




int writeRealArrToFile( 
                        char filename[],
                        real *arr,
                        int size_arr
                      )
{
/*
    Write a file with the follwing a column of arr. Matching the
    order of rows in the read grid file from coupled software.
*/
FILE *fp = NULL;
int i = 0;
int state = _STATE_OK;

if ((fp = fopen(filename, "w")) == NULL)
{
    Message("Error (writeRealArrToFile()): Unable to open %s \n", filename);
    state = _STATE_ERROR;
}
else
{
    if(size_arr > 0)
    {
        if ( arr == NULL)
        {
            Message("Error (writeRealArrToFile()): Memory access no values "
                    "available in arr!\n");
            state = _STATE_ERROR;
        }
        else
        {
            for (i = 0; i < size_arr; ++i)
            {
                if(IS_SINGLE_PRECISION)
                {
                    fprintf(fp, "%f\n", arr[i]);
                }
                else
                {
                    fprintf(fp, "%lf\n", arr[i]);
                }
            }
        }
    }
    else
    {
        Message("Warning (writeRealArrToFile()): Array size is zero!\n");
        state = _STATE_ERROR;
    }
    fclose(fp);
}

return state;
}


int writeIntegerArrToFile(
                            char filename[],
                            int *arr,
                            int size_arr
                         )
{
/*
    Write a file with the follwing a column of arr. Matching the
    order of rows in the read grid file from coupled software.
*/
int state = _STATE_OK;
FILE *fp = NULL;
int i = 0;

if ((fp = fopen(filename, "w")) == NULL)
{
    Message("\n Warning (in writeIntArrToFile()): Unable to open %s \n", filename);
    state = _STATE_ERROR;
}
else
{
    if(size_arr > 0)
    {
        if ( arr == NULL)
        {
            Message("Error (in writeIntArrToFile()): Memory access no values "
                    "available in arr!\n");
            state = _STATE_ERROR;
        }
        else
        {
            for (i = 0; i < size_arr; ++i)
            {
                fprintf(fp, "%i\n", arr[i]);
            }
        }
    }
    else
    {
        Message("Warning (in writeIntArrToFile()): Array size is zero!\n");
        state = _STATE_ERROR;
    }
    fclose(fp);
}

return state;
}
