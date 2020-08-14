/*
Utility functions for handling of sync.txt file in the exchange (XC) folder

License (MIT):

Copyright (c) 2016-2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "vof_pc_file_sync.h"


DEFINE_ON_DEMAND(debug_sleep)
{
    sleeptest();
}


void sleeptest()
{
    #if RP_HOST
    int sleeptime = (int) ceil(COUPLING_SLEEP_TIME_IN_S*1000);

    Message("Sleeping for %i ms\n", sleeptime);
    #if LINUX
    sleep((unsigned int) ceil(COUPLING_SLEEP_TIME_IN_S));
    #else
    Sleep((int) ceil(COUPLING_SLEEP_TIME_IN_S*1000));
    #endif
    #endif
}


int sync_coupling_state_from_file(char filename[], int *coupling_state)
{
    /*
        Read single value from file (filename):
        Single value has to be a float type because of ANSYS APDL...
    */
        FILE *fp = NULL;
        int state = _STATE_OK;
        real read_sync_state = -1;

        if ((fp = fopen(filename, "r")) == NULL)
        {
            Message("Warning (sync_coupling_state_from_file()): Unable to "
                    "open %s for reading\n", filename);
            state = _STATE_ERROR;
        }
        else
        {
            #if VOF_PC_DEBUG
            Message("Debug Info (sync_coupling_state_from_file()): Coupling sync "
                    "status check. \n");
            #endif

            if (IS_SINGLE_PRECISION)
            {
                if(fscanf(fp, "%f", &read_sync_state) == 1)
                {}
                else
                {
                    Message("Warning (sync_coupling_state_from_file()): Unable to "
                        "read sync status from sync file %s\n", filename);
                    state = _STATE_ERROR;
                }
                
            }
            else
            {
                if(fscanf(fp, "%lf", &read_sync_state) == 1)
                {}
                else
                {
                    Message("Warning (sync_coupling_state_from_file()): Unable to "
                            "read sync status from sync file %s\n", filename);
                    state = _STATE_ERROR;
                }
            }
            
            fclose(fp);

            if(MAC_IS_EQUAL(read_sync_state,ANSYS_READY))
            {
                *coupling_state = ANSYS_READY;
            }
            else if(MAC_IS_EQUAL(read_sync_state,FLUENT_READY))
            {
                *coupling_state =  FLUENT_READY;
            }
            else if(MAC_IS_EQUAL(read_sync_state,STOP_SIM))
            {
                *coupling_state = STOP_SIM;
            }
            else
            {
                *coupling_state = SYNC_ERROR;
            }
        }

        return state;
}

int sync_coupling_state_to_file(char filename[], int coupling_state)
{
    /*
        Writes given coupling state (state) to file (filename)
    */
    FILE* fp = NULL;
    int state = _STATE_OK;

    if ((fp = fopen(filename, "w")) == NULL) 
    {
        Message("\n Warning (sync_coupling_state_to_file()): Unable to open sync "
                "file %s for writing!\n", filename);
        state = _STATE_ERROR;
    }
    else 
    {
        Message("\nInfo (sync_coupling_state_to_file()): Coupling state written "
                "to sync file %s...\n", filename);
        fprintf(fp, "%i", coupling_state);
        fclose(fp);
    }

    return state;
}


int sync_wait_for_state(char filename[], int desired_state, int *actual_state)
{
    int i,j;
    int state = _STATE_OK;
    bool state_synced = false;
    int coupling_state = -100;

    i= 0;
    while (coupling_state != desired_state && coupling_state != STOP_SIM) 
    {
        state_synced = false;
        j = 0;
        while(!state_synced)
        {
            state = sync_coupling_state_from_file(filename, &coupling_state);
            if(state == _STATE_ERROR)
            {   
                Message("Warning (sync_wait_for_coupling()): Sync state could not "
                        "be read, trying again!\n");
                #if LINUX
                sleep((unsigned int) ceil(COUPLING_SLEEP_TIME_IN_S));
                #else
                Sleep((int) ceil(COUPLING_SLEEP_TIME_IN_S*1000));
                #endif
            }
            else
            {
                state_synced = true;
                break;
            }
            ++j;
        }

        if(i == MAX_COUPLING_TRIALS)
        {
            Message("Warning (sync_wait_for_coupling()): Number of max coupling "
                    "iterations reached, continue Fluent without ANSYS coupling "
                    "in this timestep!\n");
            state = _STATE_WARNING;
            coupling_state = desired_state;
            break;
        }
        else 
        {
            #if LINUX
            sleep((unsigned int) ceil(COUPLING_SLEEP_TIME_IN_S));
            #else
            Sleep((int) ceil(COUPLING_SLEEP_TIME_IN_S*1000));
            #endif

            if (i%20 == 0)
            {
                Message("Info (sync_wait_for_coupling()): Waiting for Ansys "
                "since %.2lf seconds ... \n", (real) (i*COUPLING_SLEEP_TIME_IN_S));
            }
            ++i;
        }
    }
    (*actual_state) = coupling_state;
    return state;
}


int sync_wait_for_coupling(char filename[])
{
    int i,j;
    int state = _STATE_OK;
    bool state_synced = false;
    int coupling_state = FLUENT_READY;

    j = 0;
    while(!state_synced)
    {
        state = sync_coupling_state_to_file(filename, FLUENT_READY);

        if(state == _STATE_ERROR)
        {
            Message("Warning (sync_wait_for_coupling()): Sync state could not "
                    "be written, trying again!\n");
            #if LINUX
            sleep((unsigned int) ceil(COUPLING_SLEEP_TIME_IN_S));
            #else
            Sleep((int) ceil(COUPLING_SLEEP_TIME_IN_S*1000));
            #endif
        }
        else
        {
            state_synced = true;
            break;
        }
        ++j;
    }

    state = sync_wait_for_state(filename, ANSYS_READY, &coupling_state);

    if(coupling_state == STOP_SIM)
    {
        Message("Info (sync_wait_for_coupling()): Coupling manually aborded. "
                "Please stop Fluent simulation");
        state = _STATE_ERROR;
    }
 return state;
}
