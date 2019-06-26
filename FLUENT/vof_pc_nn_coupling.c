/*
C-Function containing the implementation of NN mapping and coupling algorithm 
between ANSYS Fluent and ANSYS Mechanical APDL.

License (MIT):

Copyright (c) 2016-2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "vof_pc_main.h"
#include "vof_pc_nn_mapping.h"
#include "vof_pc_file_sync.h"
#include "vof_pc_read_ansys.h"


/* Define Globals */
#if RP_HOST
int *_g_ansys_to_fluent_mapping = NULL;
int *_g_fluent_to_ansys_mapping = NULL;
int _g_n_ansys_elements = -1;
#endif

int _g_n_fluent_cells = -1;

int *_g_fluent_cell_id_arr_full = NULL;
int *_g_fluent_compute_node_id_arr_full = NULL;
int *_g_cells_per_node = NULL;

/* Functions */
int initNNCoupling();
int strictCoupling();
int debug_setAnsysReady();
void freeGlobalArrays();


DEFINE_ON_DEMAND(debug_InitNNCoupling_oD)
{
    int state = initNNCoupling();
}

/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(debug_SetAnsysReady_oD)
{
    int state = debug_setAnsysReady();
}

/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(ResetSyncState_oD)
{
    #if RP_HOST
    int state = _STATE_OK;
    Message("Resetting sync state please make sure running Ansys first before "
            "running Init_NN_Coupling_oD.\n");

    state = sync_coupling_state_to_file(_SYNC_DAT_, COUPLING_INIT);

    if (state == _STATE_OK)
    {
        Message("Init sync file done!\n");
    }
    else
    {
        Message("Error initializing sync file!\n");
    }
    #endif
}

/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(debug_outputVofToAnsys_oD)
{
    int state = _STATE_OK;

    state = outputVofToAnsys();
}

/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(initCouplingWithAnsysCoords_oD)
{
    int state = _STATE_OK;
    real *joule_heat_to_fluent = NULL;
    real (*lorentz_force_to_fluent)[ND_ND] = NULL;
    real *lorentz_force_to_fluent_dir_temp = NULL;
    int i,j;

    #if RP_HOST
    real *joule_heat_from_ansys;
    real (*lorentz_force_from_ansys)[ND_ND];

    Message("Make sure that Ansys has printed out the actual element coordinates "
    "to %s for this function to work proper!. If this has not happend bevore you "
    "see this text please reinit!\n", _ANSYS_ALLOUT_DAT_);

    #endif

    state = initNNCoupling();

    #if RP_HOST

    if (state != _STATE_ERROR)
    {
        state = sync_wait_for_coupling(_SYNC_DAT_);
    }
    else
    {
        Message("Error: In coupling initialization!\n");
    }

    if (state != _STATE_ERROR)
    {
        state = readJouleHeatAndLorentzForcesFromAnsysOut(
                                                        _ANSYS_TO_FLUENT_OUT_DAT_,
                                                        &lorentz_force_from_ansys,
                                                        &joule_heat_from_ansys,
                                                        _g_n_ansys_elements
                                                    );
    }
    else
    {
        Message("Error: Syncing between ANSYS APDL and Fluent!\n");
    }
    

    if(state != _STATE_ERROR)
    {

        joule_heat_to_fluent = (real *) calloc(_g_n_fluent_cells, sizeof(real));
        lorentz_force_to_fluent = (real (*)[ND_ND]) calloc(ND_ND *_g_n_fluent_cells, sizeof(real));

        if(joule_heat_to_fluent != NULL)
        {
            state =  reorderRealArr(
                                joule_heat_from_ansys, 
                                _g_n_ansys_elements, 
                                _g_fluent_to_ansys_mapping,
                                joule_heat_to_fluent, 
                                _g_n_fluent_cells
                            );

            if(state != _STATE_ERROR && lorentz_force_to_fluent != NULL)
            {
                state =  reorderRealND_ND_Arr(
                                                lorentz_force_from_ansys, 
                                                _g_n_ansys_elements, 
                                                _g_fluent_to_ansys_mapping,
                                                lorentz_force_to_fluent, 
                                                _g_n_fluent_cells
                                            );
            }
        }
        else
        {
            state = _STATE_ERROR;
            Message("Error memory allocation!\n");
        }
        
    }
    else
    {
        Message("Error reading %s in readJouleHeatAndLorentzForcesFromAnsysOut()!\n",
                _ANSYS_TO_FLUENT_OUT_DAT_);
    }

    if(state != _STATE_ERROR)
    {
        Message("Host arrays mapped!\n");
    }

    #endif

    host_to_node_int_1(state);

    if(state != _STATE_ERROR)
    {
        state = safeDistributeMappedArrayToNodesInCellZone(
                                                joule_heat_to_fluent,
                                                _g_n_fluent_cells,
                                                Jouleheat,
                                                _g_fluent_compute_node_id_arr_full,
                                                _g_fluent_cell_id_arr_full,
                                                _g_cells_per_node,
                                                FLUID_ID
                                                );

        #if RP_HOST
        lorentz_force_to_fluent_dir_temp = (real *) calloc(_g_n_fluent_cells, sizeof(real));

        if(lorentz_force_to_fluent_dir_temp == NULL)
        {
            Message("Error memory allocation!\n");
            state = _STATE_ERROR;
        }
        #endif

        host_to_node_int_1(state);

        if(state != _STATE_ERROR)
        {
            for(i=0; i<ND_ND; ++i)
            {   
                #if RP_HOST
                for(j = 0; j<_g_n_fluent_cells; ++j)
                {
                    lorentz_force_to_fluent_dir_temp[j] = lorentz_force_to_fluent[j][i];
                }
                #endif

                if(state != _STATE_ERROR)
                { 
                    state = safeDistributeMappedArrayToNodesInCellZone(
                                                        lorentz_force_to_fluent_dir_temp,
                                                        _g_n_fluent_cells,
                                                        LFx+i,
                                                        _g_fluent_compute_node_id_arr_full,
                                                        _g_fluent_cell_id_arr_full,
                                                        _g_cells_per_node,
                                                        FLUID_ID
                                                        );
                }
            }

            if(lorentz_force_to_fluent_dir_temp != NULL)
            {   
                free(lorentz_force_to_fluent_dir_temp);
                lorentz_force_to_fluent_dir_temp = NULL;
            }
        }
        
    }

    #if RP_HOST
    if(joule_heat_from_ansys != NULL)
    {
        free(joule_heat_from_ansys);
    }
    if(lorentz_force_from_ansys != NULL)
    {
        free(lorentz_force_from_ansys);
    }
    #endif

    if(joule_heat_to_fluent != NULL)
    {
        free(joule_heat_to_fluent);
    }
    if(lorentz_force_to_fluent != NULL)
    {
        free(lorentz_force_to_fluent);
    }
}

/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(debug_freeNNCouplingGlobalArrays_oD)
{
    freeGlobalArrays();
}

/* ------------------------------------------------------------------------- */

DEFINE_EXECUTE_AT_END(strictCoupling_aE)
{
    int state = _STATE_OK;

    state = strictCoupling();

    if(state == _STATE_ERROR)
    {
        Message("Error within strictCoupling()!\n");
    }
}

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */


int initNNCoupling()
{
    int state = _STATE_OK;
    int controllSum = -1;
    int coupling_state = 0;

    #if RP_HOST
        real (*ansys_coord_arr)[ND_ND];
    #endif /* RP_HOST */

    real (*fluent_coord_arr_full)[ND_ND] = NULL;
    real *fluent_vof_arr_full = NULL;

    #if RP_HOST
        state = sync_coupling_state_from_file(_SYNC_DAT_, &coupling_state);
    #endif

    host_to_node_int_2(state, coupling_state);

    if (coupling_state == ANSYS_READY && state == _STATE_OK)
    {
        hostGetCouplingFieldsFromNodesinCellZone( 
                                                    &fluent_coord_arr_full,
                                                    &fluent_vof_arr_full, 
                                                    &_g_fluent_cell_id_arr_full,
                                                    &_g_fluent_compute_node_id_arr_full,
                                                    &_g_n_fluent_cells,
                                                    FLUID_ID
                                                );

        #if RP_HOST
        if(
            fluent_coord_arr_full == NULL ||
            fluent_vof_arr_full == NULL ||
            _g_fluent_cell_id_arr_full == NULL ||
            _g_fluent_compute_node_id_arr_full == NULL ||
            _g_n_fluent_cells < 1
        )
        {
            Message("Error Fluent Data could not be get, aborting!\n");
            state = _STATE_ERROR;
        }
        else
        {
            Message("Info: Please be sure that you have modified the header file "
                    "vof_pc_main and that you have exported the ANSYS EMAG mesh "
                    "coordinates as demanded by the instructions!\n");

            state = readCoordinatesAnsysAllOut(_ANSYS_ALLOUT_DAT_, 
                                        &ansys_coord_arr, 
                                        &_g_n_ansys_elements);

            if(
                ansys_coord_arr == NULL || 
                _g_n_ansys_elements < 1 ||
                state != _STATE_OK
            )
            {
                Message("Error Ansys Data could not be read, aborting!\n");
            }
            else
            {
                nearestNeighborMatching(
                                            fluent_coord_arr_full,
                                            _g_n_fluent_cells,
                                            ansys_coord_arr,
                                            _g_n_ansys_elements, 
                                            &_g_fluent_to_ansys_mapping, 
                                            &_g_ansys_to_fluent_mapping
                                        );

                if(
                    _g_fluent_to_ansys_mapping == NULL ||
                    _g_ansys_to_fluent_mapping == NULL ||
                    _g_n_ansys_elements < 1 ||
                    _g_n_fluent_cells < 1
                )
                {
                    Message("Error in NN Coupling, aborting!\n");
                    state = _STATE_ERROR;
                }
            }
        }
        #endif /* RP_HOST */

        host_to_node_int_1(state);

        if(state != _STATE_ERROR)
        {
            state = hostGetCellCountPerNodeInCellZone( 
                                &_g_cells_per_node,
                                &controllSum,
                                FLUID_ID
                                );

            #if RP_HOST

            state =  debugWriteMappings(
                                _g_fluent_to_ansys_mapping, 
                                _g_n_fluent_cells,
                                _g_ansys_to_fluent_mapping,
                                _g_n_ansys_elements, 
                                _FLUENT_TO_ANSYS_MAPPING_DAT_,
                                _ANSYS_TO_FLUENT_MAPPING_DAT_
                            );
            #endif

            host_to_node_int_1(state);

            if(state != _STATE_ERROR)
            {
                state = outputVofToAnsys();
            }
            else
            {
                #if RP_HOST
                    Message("Error in debugWriteMappings()!\n");
                #endif
            }

            #if RP_HOST
             state = sync_coupling_state_to_file(_SYNC_DAT_, FLUENT_READY);

            if(state != _STATE_ERROR)
            {
                Message("Everything OK!. Ansys Coupling Loop should be started now! "
                        "Then Fluent Coupling Can be started!\n");
            }
            #endif
            
        }

    }
    else
    {
        #if RP_HOST
            if(state != _STATE_OK)
                Message("Error: Error reading sync file %s!\n", _SYNC_DAT_);
            else
                Message("Error: Wrong sync status: %i, status should be 1!\n"
                        , coupling_state);
        #endif
    }

    if(state == _STATE_ERROR)
    {
        Message("Free globary arrays due to previous error!\n");
        freeGlobalArrays();
    }

    #if RP_HOST
    if(ansys_coord_arr != NULL)
    {
        free(ansys_coord_arr);
    }
    #endif

    if(fluent_coord_arr_full != NULL)
    {
        free(fluent_coord_arr_full);
    }
    if(fluent_vof_arr_full != NULL)
    {
        free(fluent_vof_arr_full);
    }

    return state;
}

/* ------------------------------------------------------------------------- */

int debug_setAnsysReady()
{
    int state = _STATE_OK;

    #if RP_HOST
    state = sync_coupling_state_to_file(_SYNC_DAT_, ANSYS_READY);

    if(state != _STATE_ERROR)
    {
        Message("Coupling state set!\n");
    }
    #endif

    return state;
}

/* ------------------------------------------------------------------------- */

int strictCoupling()
{
    int state = _STATE_OK;
    real *joule_heat_to_fluent = NULL;
    real (*lorentz_force_to_fluent)[ND_ND] = NULL;
    real *lorentz_force_to_fluent_dir_temp = NULL;
    int i,j;

    #if RP_HOST
    real *joule_heat_from_ansys;
    real (*lorentz_force_from_ansys)[ND_ND];
    #endif

    /* export fluent vof's */
    state = outputVofToAnsys();

    #if RP_HOST
    state = sync_coupling_state_to_file(_SYNC_DAT_, FLUENT_READY);

    if(state != _STATE_ERROR)
    {
        Message("Fluent Ready! Coupling state set.\n");
    }
    #endif

    /* wait for ansys */
    #if RP_HOST
    state = sync_wait_for_coupling(_SYNC_DAT_);

    if(state == _STATE_OK)
    {
    /* read ansys results */
        state = readJouleHeatAndLorentzForcesFromAnsysOut(
                                                        _ANSYS_TO_FLUENT_OUT_DAT_,
                                                        &lorentz_force_from_ansys,
                                                        &joule_heat_from_ansys,
                                                        _g_n_ansys_elements
                                                    );
        
        if(state != _STATE_ERROR)
        {
            joule_heat_to_fluent = (real *) calloc(_g_n_fluent_cells, sizeof(real));
            lorentz_force_to_fluent = (real (*)[ND_ND]) calloc(ND_ND * _g_n_fluent_cells, sizeof(real));

            if(joule_heat_to_fluent != NULL)
            {
                state =  reorderRealArr(
                                    joule_heat_from_ansys, 
                                    _g_n_ansys_elements, 
                                    _g_fluent_to_ansys_mapping,
                                    joule_heat_to_fluent, 
                                    _g_n_fluent_cells
                                );

                if(state != _STATE_ERROR && lorentz_force_to_fluent != NULL)
                {
                    state =  reorderRealND_ND_Arr(
                                                    lorentz_force_from_ansys, 
                                                    _g_n_ansys_elements, 
                                                    _g_fluent_to_ansys_mapping,
                                                    lorentz_force_to_fluent, 
                                                    _g_n_fluent_cells
                                                );
                }
            }
            else
            {
                state = _STATE_ERROR;
                Message("Error memory allocation!\n");
            }
            
        }
        else
        {
            Message("Error reading %s in readJouleHeatAndLorentzForcesFromAnsysOut()!\n",
                    _ANSYS_TO_FLUENT_OUT_DAT_);
        }
    }
    #endif

    host_to_node_int_1(state);

    /* distribute ansys results to nodes */
    if(state != _STATE_ERROR)
    {
        state = safeDistributeMappedArrayToNodesInCellZone(
                                                joule_heat_to_fluent,
                                                _g_n_fluent_cells,
                                                Jouleheat,
                                                _g_fluent_compute_node_id_arr_full,
                                                _g_fluent_cell_id_arr_full,
                                                _g_cells_per_node,
                                                FLUID_ID
                                                );

        #if RP_HOST
        lorentz_force_to_fluent_dir_temp = (real *) calloc(_g_n_fluent_cells, sizeof(real));

        if(lorentz_force_to_fluent_dir_temp == NULL)
        {
            Message("Error memory allocation!\n");
            state = _STATE_ERROR;
        }
        #endif

        host_to_node_int_1(state);

        if(state != _STATE_ERROR)
        {
            for(i=0; i<ND_ND; ++i)
            {   
                #if RP_HOST
                for(j = 0; j<_g_n_fluent_cells; ++j)
                {
                    lorentz_force_to_fluent_dir_temp[j] = lorentz_force_to_fluent[j][i];
                }
                #endif

                if(state != _STATE_ERROR)
                { 
                    state = safeDistributeMappedArrayToNodesInCellZone(
                                                        lorentz_force_to_fluent_dir_temp,
                                                        _g_n_fluent_cells,
                                                        LFx+i,
                                                        _g_fluent_compute_node_id_arr_full,
                                                        _g_fluent_cell_id_arr_full,
                                                        _g_cells_per_node,
                                                        FLUID_ID
                                                        );
                }
            }

            if(lorentz_force_to_fluent_dir_temp != NULL)
            {
                free(lorentz_force_to_fluent_dir_temp);
                lorentz_force_to_fluent_dir_temp = NULL;
            }
        }
    }

    #if RP_HOST

    if(state != _STATE_ERROR)
    {
        Message("ANSYS to Fluent success.\n");
    }

    if(joule_heat_from_ansys != NULL)
    {
        free(joule_heat_from_ansys);
    }
    if(lorentz_force_from_ansys != NULL)
    {
        free(lorentz_force_from_ansys);
    }
    #endif

    if(joule_heat_to_fluent != NULL)
    {
        free(joule_heat_to_fluent);
    }
    if(lorentz_force_to_fluent != NULL)
    {
        free(lorentz_force_to_fluent);
    }

    return state;
}

/* ------------------------------------------------------------------------- */

int outputVofToAnsys()
{
    int state = _STATE_OK;
    real (*fluent_coord_arr_full)[ND_ND] = NULL;
    real *fluent_vof_arr_full = NULL;
    int *fluent_cell_id_arr_full = NULL;
    int *fluent_compute_node_id_arr_full = NULL;
    int fluent_arr_full_size = -1;

    hostGetCouplingFieldsFromNodesinCellZone(
                                                &fluent_coord_arr_full,
                                                &fluent_vof_arr_full, 
                                                &fluent_cell_id_arr_full,
                                                &fluent_compute_node_id_arr_full,
                                                &fluent_arr_full_size,
                                                FLUID_ID
                                            );

    #if RP_HOST
    if(
        fluent_coord_arr_full == NULL ||
        fluent_vof_arr_full == NULL ||
        fluent_cell_id_arr_full == NULL ||
        fluent_compute_node_id_arr_full == NULL ||
        fluent_arr_full_size < 1
    )
    {
        Message("Error Fluent Data could not be get, aborting!\n");
        state = _STATE_ERROR;
    }
    else
    {
        state = writePropertyArrayToFileMapped(
                                        _FLUENT_TO_ANSYS_VOFOUT_DAT_,
                                        fluent_vof_arr_full,
                                        fluent_arr_full_size,
                                        _g_ansys_to_fluent_mapping,
                                        _g_n_ansys_elements
                                    );

        if(state != _STATE_OK)
        {
            state = _STATE_ERROR;
            Message("Error (outputVofToAnsys()): Error in writePropertyArrayToFileMapped()!\n");
        }
        else
        {
            Message("Info (outputVofToAnsys()): Done!\n");
        }
    }
    #endif

    if(fluent_coord_arr_full != NULL)
    {
        free(fluent_coord_arr_full);
    }
    if(fluent_vof_arr_full != NULL)
    {
        free(fluent_vof_arr_full);
    }
    if(fluent_cell_id_arr_full != NULL)
    {
        free(fluent_cell_id_arr_full);
    }
    if(fluent_compute_node_id_arr_full != NULL)
    {
        free(fluent_compute_node_id_arr_full);
    }

    return state;
}

/* ------------------------------------------------------------------------- */

void freeGlobalArrays()
{
    #if RP_HOST
    if (_g_fluent_to_ansys_mapping != NULL)
    {
        free(_g_fluent_to_ansys_mapping);
        _g_fluent_to_ansys_mapping = NULL;
    }
    if (_g_ansys_to_fluent_mapping != NULL)
    {
        free(_g_ansys_to_fluent_mapping);
        _g_ansys_to_fluent_mapping = NULL;
    }
    #endif

    if (_g_fluent_cell_id_arr_full != NULL)
    {
        free(_g_fluent_cell_id_arr_full);
        _g_fluent_cell_id_arr_full = NULL;
    }
    if (_g_fluent_compute_node_id_arr_full != NULL)
    {
        free(_g_fluent_compute_node_id_arr_full);
        _g_fluent_compute_node_id_arr_full = NULL;
    }
    if (_g_cells_per_node != NULL)
    {
        free(_g_cells_per_node);
        _g_cells_per_node = NULL;
    }

   #if RP_HOST
    Message("Global arrays set free\n");
   #endif 
}

/* ------------------------------------------------------------------------- */
