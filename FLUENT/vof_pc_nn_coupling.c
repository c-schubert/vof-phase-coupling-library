/*
C-Function containing the specific implementation of NN mapping and coupling algorithm 
between ANSYS Fluent and ANSYS Mechanical APDL. Can be modified to fit special case i.e. 
multi region coupling. 

License (MIT):

Copyright (c) 2016-2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "vof_pc_main.h"
#include "vof_pc_case.h"
#include "vof_pc_nn_mapping.h"
#include "vof_pc_file_sync.h"
#include "vof_pc_read_ansys.h"


const int _g_f_lf_udmi_vec[ND_ND] = {UDM_LFx, UDM_LFy, UDM_LFz};

const int _g_no_coupled_areas = 3;

const int _g_cell_zone_id[3] = {22,23,24};

const int _g_a2f_coupling_for_zone[3] = {1,1,1};
const int _g_f2a_coupling_for_zone[3] = {1,0,0};

const int _g_a2f_coupled_properties[3] = {JOULE_HEAT_PLUS_LORENTZ, JOULE_HEAT, JOULE_HEAT};
const int _g_f2a_coupled_properties[3] = {VOF, NONE, NONE}; /* Array with function pointers to cell marocs - TODO */ 


char _g_a_coupling_files_coords[3][250] = {_ANSYS_TO_FLUENT_MIXTURE_COORDS_OUT_DAT_, 
                                              _ANSYS_TO_FLUENT_SKIN_COORDS_OUT_DAT_,
                                              _ANSYS_TO_FLUENT_MOULD_COORDS_OUT_DAT_};
char _g_a_vol_val_files_jouleheat[3][250] =   {_ANSYS_TO_FLUENT_MIXTURE_JH_OUT_DAT_, 
                                              _ANSYS_TO_FLUENT_SKIN_JH_OUT_DAT_,
                                              _ANSYS_TO_FLUENT_MOULD_JH_OUT_DAT_};

char _g_a_vec_files_lorentzforce[3][250] =  {_ANSYS_TO_FLUENT_MIXTURE_LF_OUT_DAT_, 
                                               _DUMMY_DAT_,
                                               _DUMMY_DAT_};


char _g_f2a_files[3][250] = {_FLUENT_TO_ANSYS_VOFOUT_DAT_, _DUMMY_DAT_, _DUMMY_DAT_};

char _g_a2f_mapping_files[3][250]= {_ANSYS_TO_FLUENT_MAPPING_MIXTURE_DAT_, 
                                              _ANSYS_TO_FLUENT_MAPPING_SKIN_DAT_,
                                              _ANSYS_TO_FLUENT_MAPPING_MOULD_DAT_};
char _g_f2a_mapping_files[3][250]= {_FLUENT_TO_ANSYS_MAPPING_DAT_, _DUMMY_DAT_, _DUMMY_DAT_};

char _g_f_debug_coords_files[3][250]={_FLUENT_DEBUG_MIXTURE_COORDS_OUT_DAT_,
                                    _FLUENT_DEBUG_SKIN_COORDS_OUT_DAT_,
                                    _FLUENT_DEBUG_MOULD_COORDS_OUT_DAT_
                                    };


/*
Value of fluent cell c proptery of coupling region r:
prop_fluent[c] = prop_ansys[_g_a2f_mappings_zone_arr[r][c]] * _g_a2f_weights_zone_arr[r][c] 
Value of ansys element e propery of coupling region r:
prop_ansys[e] = prop_ansys[_g_f2a_mappings_zone_arr[r][c]] * _g_f2a_weights_zone_arr[r][c] 
*/
int **_g_a2f_mappings_zone_arr = NULL;
real **_g_a2f_weights_zone_arr = NULL;
int **_g_f2a_mappings_zone_arr = NULL;
real **_g_f2a_weights_zone_arr = NULL;
int *_g_no_a_elems_zone_arr = NULL;
int *_g_no_f_cells_zone_arr = NULL;

int **_g_f_ordered_cids_zone_arr = NULL; /* cell ids ordered from node 0 to node p for each fluent cell c*/
int **_g_f_ordered_myids_zone_arr = NULL; /*  node ids from 0 to p for each fluent cell c*/
int **_g_f_no_cells_per_node_zone_arr = NULL; /* cells per node for each node myid*/


/* Functions */
void initVofUDM();
real maxRelChangeVOFzones();
void updateOldVOFCouplingValues();
int initNNCouplingOfCellZones();
int initNNCouplingOfCellZone(
                            int **f2a_mappings_zone,
                            real **f2a_weights_zone,
                            int **a2f_mappings_zone,
                            real **a2f_weights_zone,
                            int **f_ordered_cids_zone,
                            int **f_ordered_myids_zone,
                            int **f_no_cells_per_node_zone,
                            int *no_f_cells_zone,
                            int *no_a_elems_zone,
                            int f_cell_zone_id, 
                            char ansys_zone_coord_file[],
                            char f2a_debug_mapping_file[],
                            char a2f_debug_mapping_file[]
                            );
int exchangeCellZones(int exch_state);
int exchangeVolumetricPropertyA2FZone(
                                        char ansys_vol_prop_file[],
                                        int *f2a_mapping_zone,
                                        int no_a_elems_zone,
                                        int no_f_cells_zone,
                                        int fluid_zone_id,
                                        int *f_compute_node_id_arr_full,
                                        int *f_cell_id_arr_full,
                                        int *f_cells_per_node,
                                        int udmi_idx
                                        );
int exchangeVecPropertyA2FZone(
                                char ansys_vec_prop_file[],
                                int *f2a_mapping_zone,
                                int no_a_elems_zone,
                                int no_f_cells_zone,
                                int fluid_zone_id,
                                int *f_compute_node_id_arr_full,
                                int *f_cell_id_arr_full,
                                int *f_cells_per_node,
                                const int udmis_vec[ND_ND]
                                );
int exchangeVolumetricPropertyF2AZone(  
                                    char f2a_vol_prop_file[],
                                    int *a2f_mapping_zone,
                                    int no_f_cells_zone,
                                    int no_a_elems_zone,
                                    real (*C_VAL_WRAPPER_FUN)(cell_t, Thread*),
                                    int fluid_zone_id
                                    );
void correctVolumetricPropertyA2F(
                                    real *a_e_prop,
                                    real *a_e_vol,
                                    int a_e_count,
                                    int fluid_zone_id,
                                    int udmi_idx
                                );
int strictCoupling();
int debug_setAnsysReady();
int hostWriteDebugCoords(
                            char filename[],
                            real (*coord_arr_full)[ND_ND],
                            int *cell_id_arr_full, 
                            int *compute_node_id_arr_full,
                            int arr_full_size
                        );
void freeGlobalArrays();


DEFINE_ON_DEMAND(debug_InitNNCoupling_oD)
{
    int state = initNNCouplingOfCellZones();
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

    state = exchangeCellZones(FLUENT_READY);
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(debug_resetVoF_UDM_oD)
{
	initVofUDM();
}
/* ------------------------------------------------------------------------- */
DEFINE_ON_DEMAND(debug_vofToUDM_oD)
{
	updateOldVOFCouplingValues();
}
/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(initCouplingWithAnsysCoords_oD)
{
    int state = _STATE_OK;
    int coupling_state = 0; 

    #if RP_HOST
    Message("If not allready done please start ANSYS now!\n");
    state = sync_wait_for_state(_SYNC_DAT_, ANSYS_READY, &coupling_state);
    #endif

    host_to_node_int_2(state, coupling_state);

    if (coupling_state == ANSYS_READY && state == _STATE_OK)
    {
        state = initNNCouplingOfCellZones();
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

    host_to_node_int_1(state);

    if(state != _STATE_ERROR && coupling_state == ANSYS_READY)
    {
        exchangeCellZones(FLUENT_READY);
    }
    else
    {
        state = _STATE_ERROR;
        #if RP_HOST
            Message("Error in initNNCouplingOfCellZones()!\n");
        #endif
    }

    #if RP_HOST
    if(state != _STATE_ERROR)
    {
        state = sync_coupling_state_to_file(_SYNC_DAT_, FLUENT_READY);
        if(state != _STATE_ERROR)
        {
        Message("Everything OK!. Ansys Coupling Loop should be started now! "
                "Then Fluent Coupling can be started!\n");
        }   
    }
    #endif

    #if RP_HOST
    if (state != _STATE_ERROR)
    {
        state = sync_wait_for_coupling(_SYNC_DAT_);
    }
    #endif

    if(state != _STATE_ERROR)
    {
        exchangeCellZones(ANSYS_READY);
    }

    initVofUDM();
    updateOldVOFCouplingValues();
}


/* ------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(debug_freeNNCouplingGlobalArrays_oD)
{
    freeGlobalArrays();
}

/* ------------------------------------------------------------------------- */
DEFINE_ON_DEMAND(debug_ReadAnsys_oD)
{
    exchangeCellZones(ANSYS_READY);
}
/* ------------------------------------------------------------------------- */
DEFINE_ON_DEMAND(debug_WriteVofToAnsys_oD)
{
    exchangeCellZones(FLUENT_READY);
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

DEFINE_EXECUTE_AT_END(looseCoupling_aE)
{
    int state = _STATE_OK;
    real maxRelChangeVof = maxRelChangeVOFzones();

    #if RP_HOST
    Message("Max rel change is: %lf\n", maxRelChangeVof);
    #endif

    if(maxRelChangeVof > VOF_MAX_REL_CHANGE)
    {
        state = strictCoupling();
        updateOldVOFCouplingValues();
    }
    else
    {   
        #if RP_HOST
        Message("No coupling this time!\n");
        #endif
    }

    if(state == _STATE_ERROR)
    {
        Message("Error within strictCoupling()!\n");
    }
}
/* ------------------------------------------------------------------------- */


DEFINE_ON_DEMAND(Debug_Export_Fluent_Coords)
{
    int state = _STATE_OK;
    int controllSum = -1;
    int *f_no_cells_per_node_zone=NULL;
    real (*f_coord_arr_full)[ND_ND] = NULL;
    int *f_ordered_cids_zone = NULL;
    int *f_ordered_myids_zone = NULL;
    int no_f_cells_zone = 0;
    int arr_full_size = -1;
    int ir;

    for(ir=0; ir<_g_no_coupled_areas;++ir)
    {
        controllSum = -1;
        no_f_cells_zone = 0;

        state = hostGetCellCountPerNodeInCellZone( 
                                                &f_no_cells_per_node_zone,
                                                &controllSum,
                                                _g_cell_zone_id[ir]
                                                );

        hostGetOrderingArraysFromNodesInCellZone( 
                                            &f_ordered_cids_zone,
                                            &f_ordered_myids_zone,
                                            &no_f_cells_zone,
                                            _g_cell_zone_id[ir]
                                            );  

        hostGetCellCoordsFromNodesInCellZone( 
                                        &f_coord_arr_full,
                                        no_f_cells_zone,
                                        _g_cell_zone_id[ir]
                                        );

        
        state = hostWriteDebugCoords (
                            _g_f_debug_coords_files[ir],
                            f_coord_arr_full,
                            f_ordered_cids_zone,
                            f_ordered_myids_zone,
                            no_f_cells_zone
                        );

        #if RP_HOST
        if(f_coord_arr_full != NULL)
        {
            free(f_coord_arr_full);
            f_coord_arr_full = NULL;
        }
        if(f_ordered_cids_zone != NULL)
        {
            free(f_ordered_cids_zone);
            f_ordered_cids_zone = NULL;
        }
        if(f_ordered_myids_zone != NULL)
        {
            free(f_ordered_myids_zone);
            f_ordered_myids_zone = NULL;
        }
        if(f_no_cells_per_node_zone != NULL)
        {
            free(f_no_cells_per_node_zone);
            f_no_cells_per_node_zone = NULL;
        }
        #endif
    }
    

    #if RP_HOST
        Message("Export Done!\n");
    #endif
}
/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */

real maxRelChangeVOFzones()
{
    real relChange = 0;

    #if RP_NODE
    cell_t c;
    Thread *t;
    Thread **pt;
    Domain *domain = Get_Domain(1); 
    real vof_diff = 0;
    int i;
  
    for (i = 0; i < _g_no_coupled_areas; i++)
    {
        if (_g_f2a_coupled_properties[i] == VOF)
        {
            t = Lookup_Thread(domain, _g_cell_zone_id[i]);
            pt = THREAD_SUB_THREADS(t);

            begin_c_loop_int(c, t) 
            {

                vof_diff = fabs(C_VOF(c,pt[COUPLING_PHASE_FRAC_IDX]) - C_UDMI(c,t, UDM_VOF_old));
                relChange = MAX(relChange,vof_diff);
            }end_c_loop_int(c, t)
        }
    }

    relChange = PRF_GRHIGH1(relChange);
    #endif

    node_to_host_real_1(relChange);

    return relChange;
}
/* ------------------------------------------------------------------------- */

void initVofUDM()
{
    #if RP_NODE
    cell_t c;
    Thread *t;
    Domain *domain = Get_Domain(1); 

    thread_loop_c(t, domain) /*loops over all cell threads in domain*/
    {
        begin_c_loop_int(c, t) 
        {
            C_UDMI(c,t, UDM_VOF_old) = 0.0 ;
        }end_c_loop_int(c, t)
    }
    #endif
}
/* ------------------------------------------------------------------------- */

void updateOldVOFCouplingValues()
{
    #if RP_NODE
    cell_t c;
    Thread *t;
    Thread **pt;
    Domain *domain = Get_Domain(1); 
    int i;
  
    for (i = 0; i < _g_no_coupled_areas; i++)
    {
        if (_g_f2a_coupled_properties[i] == VOF)
        {
            t = Lookup_Thread(domain, _g_cell_zone_id[i]);
            pt = THREAD_SUB_THREADS(t);

            begin_c_loop_int(c, t) 
            {
                C_UDMI(c,t, UDM_VOF_old) = C_VOF(c,pt[COUPLING_PHASE_FRAC_IDX]);

            }end_c_loop_int(c, t)
        }
    }
    #endif

}
/* ------------------------------------------------------------------------- */

int initNNCouplingOfCellZones()
{
    int ir;
    int state = _STATE_OK;

    _g_a2f_mappings_zone_arr = (int **) malloc(_g_no_coupled_areas*sizeof(int *));
    _g_a2f_weights_zone_arr = (real **) malloc(_g_no_coupled_areas*sizeof(real *));
    _g_f2a_mappings_zone_arr = (int **) malloc(_g_no_coupled_areas*sizeof(int *));
    _g_f2a_weights_zone_arr = (real **) malloc(_g_no_coupled_areas*sizeof(real *));

    _g_no_a_elems_zone_arr = (int *) malloc(_g_no_coupled_areas*sizeof(int));
    _g_no_f_cells_zone_arr = (int *) malloc(_g_no_coupled_areas*sizeof(int));

    _g_f_ordered_cids_zone_arr = (int **) malloc(_g_no_coupled_areas*sizeof(int *));
    _g_f_ordered_myids_zone_arr = (int **) malloc(_g_no_coupled_areas*sizeof(int *));
    _g_f_no_cells_per_node_zone_arr = (int **) malloc(_g_no_coupled_areas*sizeof(int *));

    for(ir = 0; ir < _g_no_coupled_areas; ++ir)
    {
        if(state != _STATE_ERROR)
        {
            state = initNNCouplingOfCellZone(
                                        &_g_f2a_mappings_zone_arr[ir],
                                        &_g_f2a_weights_zone_arr[ir],
                                        &_g_a2f_mappings_zone_arr[ir],
                                        &_g_a2f_weights_zone_arr[ir],
                                        &_g_f_ordered_cids_zone_arr[ir],
                                        &_g_f_ordered_myids_zone_arr[ir],
                                        &_g_f_no_cells_per_node_zone_arr[ir],
                                        &_g_no_a_elems_zone_arr[ir],
                                        &_g_no_f_cells_zone_arr[ir],
                                        _g_cell_zone_id[ir],
                                        _g_a_coupling_files_coords[ir],
                                        _g_f2a_mapping_files[ir],
                                        _g_a2f_mapping_files[ir]
                                        );
        }
        else
        {
            Message("Error initNNCouplingOfCellZones(): initializing state for "
                    "coupling region: %i\n", ir);
        }
    }

    #if RP_NODE
    /* not necessary on nodes */
    free(_g_a2f_mappings_zone_arr);
    free(_g_a2f_weights_zone_arr);
    free(_g_f2a_mappings_zone_arr);
    free(_g_f2a_weights_zone_arr);

    free(_g_no_a_elems_zone_arr);
    free(_g_no_f_cells_zone_arr);

    free(_g_f_ordered_cids_zone_arr);
    free(_g_f_ordered_myids_zone_arr);
    free(_g_f_no_cells_per_node_zone_arr);
    #endif

    return state;
}
/* ------------------------------------------------------------------------- */


int initNNCouplingOfCellZone(
                            int **f2a_mappings_zone,
                            real **f2a_weights_zone,
                            int **a2f_mappings_zone,
                            real **a2f_weights_zone,
                            int **f_ordered_cids_zone,
                            int **f_ordered_myids_zone,
                            int **f_no_cells_per_node_zone,
                            int *no_a_elems_zone,
                            int *no_f_cells_zone,
                            int f_cell_zone_id, 
                            char ansys_zone_coord_file[],
                            char f2a_debug_mapping_file[],
                            char a2f_debug_mapping_file[]
                            )
{
    int state = _STATE_OK;
    int controllSum = -1;

    #if RP_HOST
    real (*a_coord_arr)[ND_ND] = NULL;
    #endif /* RP_HOST */

    real (*f_coord_arr_full)[ND_ND] = NULL;
    real (*test_arr) = NULL;

    state = hostGetCellCountPerNodeInCellZone( 
                                                f_no_cells_per_node_zone,
                                                &controllSum,
                                                f_cell_zone_id
                                                );
    
        
    hostGetOrderingArraysFromNodesInCellZone( 
                                            f_ordered_cids_zone,
                                            f_ordered_myids_zone,
                                            no_f_cells_zone,
                                            f_cell_zone_id
                                            );
    #if RP_HOST
    Message("Cells in zone: %i \n", (*no_f_cells_zone));
    #endif

    hostGetCellCoordsFromNodesInCellZone( 
                                            &f_coord_arr_full,
                                            (*no_f_cells_zone),
                                            f_cell_zone_id
                                            );

    #if RP_HOST
    if(
        state == _STATE_ERROR ||
        f_coord_arr_full == NULL ||
        f_ordered_cids_zone == NULL ||
        f_ordered_myids_zone == NULL ||
        (*no_f_cells_zone) < 1
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

        state = readCoordinatesFromAnsysOut(ansys_zone_coord_file, 
                                            &a_coord_arr, 
                                            no_a_elems_zone);

        if(
            a_coord_arr == NULL || 
            (*no_a_elems_zone) < 1 ||
            state == _STATE_ERROR
        )
        {
            Message("Error Ansys Data could not be read, aborting!\n");
        }
        else
        {
            nearestNeighborMatching(
                                    f_coord_arr_full,
                                    (*no_f_cells_zone),
                                    a_coord_arr,
                                    (*no_a_elems_zone), 
                                    f2a_mappings_zone, 
                                    f2a_weights_zone, 
                                    a2f_mappings_zone,
                                    a2f_weights_zone
                                );

            if(
                f2a_mappings_zone == NULL ||
                a2f_mappings_zone == NULL ||
                (*no_a_elems_zone) < 1 ||
                (*no_f_cells_zone) < 1
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
        #if RP_HOST
        state =  debugWriteMappings(
                            (*f2a_mappings_zone), 
                            (*no_f_cells_zone),
                            (*a2f_mappings_zone),
                            (*no_a_elems_zone), 
                            f2a_debug_mapping_file,
                            a2f_debug_mapping_file
                        );
        #endif
    }

    if(state == _STATE_ERROR)
    {
        Message("Free globary arrays due to previous error!\n");
        freeGlobalArrays();
    }

    #if RP_HOST
    if(a_coord_arr != NULL)
    {
        free(a_coord_arr);
    }
    #endif

    if(f_coord_arr_full != NULL)
    {
        free(f_coord_arr_full);
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

int exchangeCellZones(int exch_state) 
{
    int state = _STATE_OK;
    int ir;

    for(ir = 0; ir < _g_no_coupled_areas; ++ir)
    {   
        if(
            (exch_state == ANSYS_READY) &&
            _g_a2f_coupling_for_zone[ir] && 
            (state != _STATE_ERROR)
          )
        {
            if (_g_a2f_coupled_properties[ir] == JOULE_HEAT_PLUS_LORENTZ ||
                _g_a2f_coupled_properties[ir] == JOULE_HEAT)
            {
                #if RP_HOST
                    Message("Exchanging Joule heat for zone %i\n", _g_cell_zone_id[ir]);
                #endif
                state = exchangeVolumetricPropertyA2FZone(                  
                                        _g_a_vol_val_files_jouleheat[ir],
                                        _g_f2a_mappings_zone_arr[ir],
                                        _g_no_a_elems_zone_arr[ir],
                                        _g_no_f_cells_zone_arr[ir],
                                        _g_cell_zone_id[ir],
                                        _g_f_ordered_myids_zone_arr[ir],
                                        _g_f_ordered_cids_zone_arr[ir],
                                        _g_f_no_cells_per_node_zone_arr[ir],
                                        UDM_JH
                                        );
            }

            if (_g_a2f_coupled_properties[ir] == JOULE_HEAT_PLUS_LORENTZ)
            {
                #if RP_HOST
                    Message("Exchanging Lorentz-Forces for zone %i\n", _g_cell_zone_id[ir]);
                #endif
                state = exchangeVecPropertyA2FZone(
                                                _g_a_vec_files_lorentzforce[ir],
                                                _g_f2a_mappings_zone_arr[ir],
                                                _g_no_a_elems_zone_arr[ir],
                                                _g_no_f_cells_zone_arr[ir],
                                                _g_cell_zone_id[ir],
                                                _g_f_ordered_myids_zone_arr[ir],
                                                _g_f_ordered_cids_zone_arr[ir],
                                                _g_f_no_cells_per_node_zone_arr[ir],
                                                _g_f_lf_udmi_vec
                                                    );
            }
        }
        if(
            (exch_state == FLUENT_READY) &&
            _g_f2a_coupling_for_zone[ir] && 
            (state != _STATE_ERROR)
           )
        {
            if (_g_f2a_coupled_properties[ir] == VOF)
            {
                state = exchangeVolumetricPropertyF2AZone(
                                                        _g_f2a_files[ir],
                                                       _g_a2f_mappings_zone_arr[ir],
                                                        _g_no_f_cells_zone_arr[ir],
                                                        _g_no_a_elems_zone_arr[ir],
                                                        get_c_vof,
                                                        _g_cell_zone_id[ir]
                                                            );
            }
        }
    }

    return state;
}

/* ------------------------------------------------------------------------- */

int exchangeVolumetricPropertyA2FZone(
                                        char ansys_vol_prop_file[],
                                        int *f2a_mapping_zone,
                                        int no_a_elems_zone,
                                        int no_f_cells_zone,
                                        int fluid_zone_id,
                                        int *f_compute_node_id_arr_full,
                                        int *f_cell_id_arr_full,
                                        int *f_cells_per_node,
                                        int udmi_idx
                                        )
{
    int state = _STATE_OK;
    real *vol_prop_to_fluent = NULL;

    real *vol_prop_from_ansys = NULL;
    real *elem_vol_from_ansys = NULL;

    #if RP_HOST
    state = readElemValueAndVolumeFromAnsysOut(
                                            ansys_vol_prop_file,
                                            &vol_prop_from_ansys,
                                            &elem_vol_from_ansys,
                                            no_a_elems_zone
                                        );
    if(state != _STATE_ERROR )
    {
        vol_prop_to_fluent = (real *) calloc(no_f_cells_zone, sizeof(real));       

        if(state != _STATE_ERROR && vol_prop_to_fluent != NULL)
        {
            state =  reorderRealArr(
                                    vol_prop_from_ansys, 
                                    no_a_elems_zone, 
                                    f2a_mapping_zone,
                                    vol_prop_to_fluent, 
                                    no_f_cells_zone
                                );
        }
        else
        {
            state = _STATE_ERROR;
            Message("Error exchangeVolumetricPropertyA2FZone()!\n");
        }
    }
    else
    {
        Message("Error reading %s in readElemValueAndVolumeFromAnsysOut()!\n",
                ansys_vol_prop_file);
    }

    if(state != _STATE_ERROR)
    {
        Message("Host arrays for fluid zone id: %i!\n", fluid_zone_id);
    }
    #endif

    host_to_node_int_1(state);

    if(state != _STATE_ERROR)
    {
        state = safeDistributeMappedArrayToNodesInCellZone(
                                                vol_prop_to_fluent,
                                                no_f_cells_zone,
                                                udmi_idx,
                                                f_compute_node_id_arr_full,
                                                f_cell_id_arr_full,
                                                f_cells_per_node,
                                                fluid_zone_id
                                                );
    }

    if(state != _STATE_ERROR)
    {
        correctVolumetricPropertyA2F(   
                                    vol_prop_from_ansys,
                                    elem_vol_from_ansys,
                                    no_a_elems_zone,
                                    fluid_zone_id,
                                    udmi_idx
                                    );
    }

    #if RP_HOST
    if(vol_prop_from_ansys != NULL)
    {
        free(vol_prop_from_ansys);
    }
    if(elem_vol_from_ansys != NULL)
    {
        free(elem_vol_from_ansys);
    }
    #endif

    if(vol_prop_to_fluent != NULL)
    {
        free(vol_prop_to_fluent);
    }

    return state;
}



int exchangeVecPropertyA2FZone(
                                char ansys_vec_prop_file[],
                                int *f2a_mapping_zone,
                                int no_a_elems_zone,
                                int no_f_cells_zone,
                                int fluid_zone_id,
                                int *f_compute_node_id_arr_full,
                                int *f_cell_id_arr_full,
                                int *f_cells_per_node,
                                const int udmis_vec[ND_ND]
                                )
{
    int state = _STATE_OK;
    real (*vec_prop_to_fluent)[ND_ND] = NULL;
    real *tmp_1D_prop_to_fluent = NULL;
    int i,j = 0;

    #if RP_HOST
        real (*vec_prop_from_ansys)[ND_ND] = NULL;

        state = readElemValueVecFromAnsysOut(
                                        ansys_vec_prop_file,
                                        &vec_prop_from_ansys,
                                        no_a_elems_zone
                                    );

        if(state != _STATE_ERROR )
        {
            vec_prop_to_fluent = (real (*)[ND_ND]) 
                                calloc(ND_ND * no_f_cells_zone, sizeof(real));     

            if(state != _STATE_ERROR && vec_prop_to_fluent != NULL)
            {
                state =  reorderRealND_ND_Arr(
                                        vec_prop_from_ansys, 
                                        no_a_elems_zone, 
                                        f2a_mapping_zone,
                                        vec_prop_to_fluent, 
                                        no_f_cells_zone
                                    );
            }
            else
            {
                state = _STATE_ERROR;
                Message("Error exchangeVolumetricPropertyA2FZone()!\n");
            }
        }
        else
        {
            Message("Error reading %s in readElemValueVecFromAnsysOut()!\n",
                    ansys_vec_prop_file);
        }

    if(state != _STATE_ERROR)
    {
        Message("Host arrays for fluid zone id: %i!\n", fluid_zone_id);
    }
    #endif

    host_to_node_int_1(state);

    if(state != _STATE_ERROR)
    {
        #if RP_HOST
        tmp_1D_prop_to_fluent = (real *) calloc(no_f_cells_zone, sizeof(real));

        if(tmp_1D_prop_to_fluent == NULL)
        {
            Message("Error exchangeVecPropertyA2FZone(): memory allocation!\n");
            state = _STATE_ERROR;
        }
        #endif

      host_to_node_int_1(state);

        if(state != _STATE_ERROR)
        {   
            for(i=0; i<ND_ND; ++i)
            { 

                #if RP_HOST
                for(j = 0; j<no_f_cells_zone; ++j)
                {
                    tmp_1D_prop_to_fluent[j] = vec_prop_to_fluent[j][i];
                }
                #endif

                if(state != _STATE_ERROR)
                { 
                    state = safeDistributeMappedArrayToNodesInCellZone(
                                                tmp_1D_prop_to_fluent,
                                                no_f_cells_zone,
                                                udmis_vec[i],
                                                f_compute_node_id_arr_full,
                                                f_cell_id_arr_full,
                                                f_cells_per_node,
                                                fluid_zone_id
                                                                    );
                }
            }

            if(tmp_1D_prop_to_fluent != NULL)
            {   
                free(tmp_1D_prop_to_fluent);
                tmp_1D_prop_to_fluent = NULL;
            }
        }
    }

    #if RP_HOST
    if(vec_prop_from_ansys != NULL)
    {
        free(vec_prop_from_ansys);
    }
    #endif

    if(vec_prop_to_fluent != NULL)
    {
        free(vec_prop_to_fluent);
    }

    return state;
}

int exchangeVolumetricPropertyF2AZone(  
                                    char f2a_vol_prop_file[],
                                    int *a2f_mapping_zone,
                                    int no_f_cells_zone,
                                    int no_a_elems_zone,
                                    real (*C_VAL_WRAPPER_FUN)(cell_t, Thread*),
                                    int fluid_zone_id
                                    )
{
    int state = _STATE_OK;
    real *f2a_vol_property = NULL;
    int f2a_vol_property_arr_size = 0;

    hostGetOrderedFieldValueArrayFromNodesInCellZone(
                                                &f2a_vol_property, 
                                                C_VAL_WRAPPER_FUN,
                                                &f2a_vol_property_arr_size,
                                                fluid_zone_id
                                            );

    #if RP_HOST
    if(f2a_vol_property == NULL)
    {
        Message("Error exchangeVolumetricPropertyF2AZone(): Fluent Data for zone id %i could "
                "not be get, aborting!\n",fluid_zone_id);
        state = _STATE_ERROR;
    }
    else
    {
        state = writePropertyArrayToFileMapped(
                                        f2a_vol_prop_file,
                                        f2a_vol_property,
                                        no_f_cells_zone,
                                        a2f_mapping_zone,
                                        no_a_elems_zone
                                    );

        if(state != _STATE_OK)
        {
            state = _STATE_ERROR;
            Message("Error (exchangeVolumetricPropertyF2AZone()): Error in" 
                    " for zone id %i writePropertyArrayToFileMapped()!\n", fluid_zone_id);
        }
        else
        {
            Message("Info (exchangeVolumetricPropertyF2AZone()): For zone id %i "
                    ",done!\n",fluid_zone_id);
        }
    }
    #endif

    if(f2a_vol_property != NULL)
    {
        free(f2a_vol_property);
    }

    return state;
}

void correctVolumetricPropertyA2F(
                                    real *a_e_prop,
                                    real *a_e_vol,
                                    int a_e_count,
                                    int fluid_zone_id,
                                    int udmi_idx
                                )
{
    /* Correction of conservative property 
    */
    real a_vol_weighted_prop_sum = 0;

    #if RP_HOST
    int i = 0;
    #endif

    #if RP_NODE
    real f_vol_weighted_prop_sum = 0;
    real corr_fac = 0;

    cell_t c;
    Thread *t;
    Domain *domain = Get_Domain(1); 
    #endif

    #if RP_HOST
    for(i=0;i<a_e_count;++i)
    {   
        a_vol_weighted_prop_sum = a_vol_weighted_prop_sum + a_e_prop[i]*a_e_vol[i];
    }
    Message("Host: Ansys sent heat: %lf W\n", a_vol_weighted_prop_sum);
    #endif

    host_to_node_real_1(a_vol_weighted_prop_sum);

    #if RP_NODE
    t = Lookup_Thread(domain, fluid_zone_id);

    begin_c_loop_int(c, t) 
    {
        f_vol_weighted_prop_sum += C_UDMI(c,t, udmi_idx)*C_VOLUME(c,t);
    }end_c_loop_int(c, t)

    f_vol_weighted_prop_sum = PRF_GRSUM1(f_vol_weighted_prop_sum);
    corr_fac = a_vol_weighted_prop_sum/f_vol_weighted_prop_sum;

    Message0("Fluent transfered heat to: %lf W\n", f_vol_weighted_prop_sum);
    Message0("Ansys sent heat: %lf W\n", a_vol_weighted_prop_sum);
    Message0("Correcting heat with factor: %lf\n",corr_fac);

    begin_c_loop_int(c, t) 
    {
        C_UDMI(c,t, udmi_idx) = corr_fac * C_UDMI(c,t, udmi_idx);
    }end_c_loop_int(c, t)
    #endif
}
/* ------------------------------------------------------------------------- */


int strictCoupling()
{
    int state = _STATE_OK;

    if(state != _STATE_ERROR)
    {
        state = exchangeCellZones(FLUENT_READY);
    }

    #if RP_HOST
    /* wait for ansys */
    if(state != _STATE_ERROR)
    {
         state = sync_wait_for_coupling(_SYNC_DAT_);
    }
    #endif /*RP_HOST*/ 

    if(state != _STATE_ERROR)
    {
        state = exchangeCellZones(ANSYS_READY);
    }

    #if RP_HOST
    if(state != _STATE_ERROR)
    {
        Message("ANSYS <-> Fluent coupling successfull!\n");
    }
    #endif

    return state;
}
/* ------------------------------------------------------------------------- */



int hostWriteDebugCoords(
                            char filename[],
                            real (*coord_arr_full)[ND_ND],
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
    if(arr_full_size > 0)
    {
        if ( cell_id_arr_full == NULL 
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
                    fprintf(fp, "%i %i %f %f %f\n", cell_id_arr_full[i],
                             compute_node_id_arr_full[i],  
                             coord_arr_full[i][0], coord_arr_full[i][1], 
                             coord_arr_full[i][2]);
                else
                    fprintf(fp, "%i %2i % #15.7lE % #15.7lE % #15.7lE\n", cell_id_arr_full[i], 
                            compute_node_id_arr_full[i], 
                            coord_arr_full[i][0], coord_arr_full[i][1], 
                            coord_arr_full[i][2]);
            }
            #else
            for (i = 0; i < arr_full_size; ++i)
            {
                if(IS_SINGLE_PRECISION)
                    fprintf(fp, "%i %i %f %f\n", cell_id_arr_full[i], 
                             compute_node_id_arr_full[i],  
                             coord_arr_full[i][0], coord_arr_full[i][1]);
                else
                    fprintf(fp, "%i %i %lf %lf\n", cell_id_arr_full[i], 
                             compute_node_id_arr_full[i], 
                             coord_arr_full[i][0], coord_arr_full[i][1]);
            }
            #endif /* RP_3D */
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
/* ------------------------------------------------------------------------- */

void freeGlobalArrays()
{
    int ir = 0;

    #if RP_HOST
    if (_g_f2a_mappings_zone_arr != NULL)
    {
        for(ir=0; ir<_g_no_coupled_areas; ++ir)
        {
            if(_g_f2a_mappings_zone_arr[ir] != NULL)
            {
                free(_g_f2a_mappings_zone_arr[ir]);
            }
        }
        free(_g_f2a_mappings_zone_arr);
    }

    if (_g_f2a_weights_zone_arr != NULL)
    {
        for(ir=0; ir<_g_no_coupled_areas; ++ir)
        {
            if(_g_f2a_weights_zone_arr[ir] != NULL)
            {
                free(_g_f2a_weights_zone_arr[ir]);
            }
        }
        free(_g_f2a_weights_zone_arr);
    }

    if (_g_a2f_mappings_zone_arr != NULL)
    {
        for(ir=0; ir<_g_no_coupled_areas; ++ir)
        {
            if(_g_a2f_mappings_zone_arr[ir] != NULL)
            {
                free(_g_a2f_mappings_zone_arr[ir]);
            }
        }
        free(_g_a2f_mappings_zone_arr);
    }

    if (_g_a2f_weights_zone_arr != NULL)
    {
        for(ir=0; ir<_g_no_coupled_areas; ++ir)
        {
            if(_g_a2f_weights_zone_arr[ir] != NULL)
            {
                free(_g_a2f_weights_zone_arr[ir]);
            }
        }
        free(_g_a2f_weights_zone_arr);
    }
    #endif

    if (_g_f_no_cells_per_node_zone_arr != NULL)
    {
        for(ir=0; ir<_g_no_coupled_areas; ++ir)
        {
            if(_g_f_no_cells_per_node_zone_arr[ir] != NULL)
            {
                free(_g_f_no_cells_per_node_zone_arr[ir]);
            }
        }
        free(_g_f_no_cells_per_node_zone_arr);
    }
    
    if (_g_f_ordered_cids_zone_arr != NULL)
    {
        for(ir=0; ir<_g_no_coupled_areas; ++ir)
        {
            if(_g_f_ordered_cids_zone_arr[ir] != NULL)
            {
                free(_g_f_ordered_cids_zone_arr[ir]);
            }
        }
        free(_g_f_ordered_cids_zone_arr);
    }
    
    if (_g_f_ordered_myids_zone_arr != NULL)
    {
        for(ir=0; ir<_g_no_coupled_areas; ++ir)
        {
            if(_g_f_ordered_myids_zone_arr[ir] != NULL)
            {
                free(_g_f_ordered_myids_zone_arr[ir]);
            }
        }
        free(_g_f_ordered_myids_zone_arr);
    }

    if (_g_no_f_cells_zone_arr != NULL)
    {
        free(_g_no_f_cells_zone_arr);
    }

    if (_g_no_a_elems_zone_arr != NULL)
    {
        free(_g_no_a_elems_zone_arr);
    }

   #if RP_HOST
    Message("Global arrays set free\n");
   #endif 
}

/* ------------------------------------------------------------------------- */
