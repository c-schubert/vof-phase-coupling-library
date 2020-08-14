/*
Functions to get fields and node distribution information from nodes to host

License (MIT):

Copyright (c) 2016-2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "vof_pc_fluent_get_fields.h" 

int hostGetCellCountPerNodeInCellZone( 
                                      int **cells_per_node,
                                      int *cell_sum_over_nodes,
                                      int cellZoneID
                                      )
{  
 /*
   Writes cell count for each compute node in cell zone 
   in domain d with ID cellZoneID to the array cells_per_node) and 
   returns sum of cell count for validation purposes.
 */

int state = _STATE_OK;
int cell_count_in_node = 0;
int node_id;
int pe;

#if RP_HOST
int i = 0;
#endif

#if !RP_HOST
cell_t c;
Domain *d = Get_Domain(1);
Thread *t;
int return_cell_count = -1;

#if RP_NODE
cell_count_in_node = 0;
t = Lookup_Thread(d, cellZoneID);
cell_count_in_node = THREAD_N_ELEMENTS_INT(t);
return_cell_count = cell_count_in_node;
node_id = myid;

#if VOF_PC_DEBUG
Message("Info (hostGetCellCountPerNodeInCellZone()):cell_count_in_node: %i, "
        "myid: %i.\n", cell_count_in_node, node_id);
#endif

PRF_GSYNC();
pe = (I_AM_NODE_ZERO_P) ? node_host : node_zero;
/*Sent data from nodes to node0 or from node0 to host*/
PRF_CSEND_INT(pe, &cell_count_in_node, 1, myid);
PRF_CSEND_INT(pe, &node_id, 1, myid);


/* node_0 now collect data sent by other compute nodes */
/*  and sends it to the host */
if (I_AM_NODE_ZERO_P)
{  
  /* pe only acts as a counter in this loop */
  compute_node_loop_not_zero (pe) 
  {
    /* Receive data */
    PRF_CRECV_INT(pe, &cell_count_in_node, 1, pe);
    PRF_CRECV_INT(pe, &node_id, 1, pe);

    /* send data */
    PRF_CSEND_INT(node_host, &cell_count_in_node, 1, myid);
    PRF_CSEND_INT(node_host, &node_id, 1, myid);
  }
}
#endif /* RP_NODE in !RP_HOST*/

return state;
#endif /* !RP_HOST*/

#if RP_HOST
*cells_per_node = NULL;

*cells_per_node = calloc(compute_node_count, sizeof(int));

if (*cells_per_node == NULL)
{
  Message("Error (hostGetCellCountPerNodeInCellZone()): Memory allocation error "
          "in hostGetCellCountPerNodeInCellZone()\n");
  state = _STATE_ERROR;
}
else
{  
  for (i = 0; i<compute_node_count; ++i)
  {
    (*cells_per_node)[i] = 0;
  }
}

/* pe only acts as a counter in this loop */
(*cell_sum_over_nodes) = 0;
compute_node_loop (pe) 
{ 
  PRF_CRECV_INT(node_zero, &cell_count_in_node, 1, node_zero);
  PRF_CRECV_INT(node_zero, &node_id, 1, node_zero);

  if (state == _STATE_OK)
  {
    if (node_id >= 0 && node_id < compute_node_count)
    {
    (*cells_per_node)[node_id] = cell_count_in_node;
    }
    else
    {
      Message("Error (hostGetCellCountPerNodeInCellZone()): Index out of bounds "
              "in hostGetCellCountPerNodeInCellZone\n");
      state = _STATE_ERROR;
    }
  }
  
  (*cell_sum_over_nodes) += cell_count_in_node;
}

return state;

#endif /* RP_HOST*/
}


real get_c_vof(cell_t c, Thread *t)
{
  Thread **pt;
  pt = THREAD_SUB_THREADS(t);
  return C_VOF(c, pt[COUPLING_PHASE_FRAC_IDX]);
}

real get_coord(cell_t c,Thread *t,int dim)
{
  real x[ND_ND];
  C_CENTROID(x,c,t);

  return x[dim];
}

real get_x_coord(cell_t c,Thread *t)
{
  return get_coord(c,t,0);
}

real get_y_coord(cell_t c,Thread *t)
{
  return get_coord(c,t,1);
}

real get_z_coord(cell_t c,Thread *t)
{
  return get_coord(c,t,2);
}


void hostGetCellCoordsFromNodesInCellZone(
                                          real (**coord_arr_full)[ND_ND],
                                          int length_arrs_full,
                                          int cell_zone
                                          )
{
  real *x_arr_full=NULL;
  real *y_arr_full=NULL;
  real *z_arr_full=NULL;
  int arr_control_size = 0;
  int i = 0;
  
  hostGetOrderedFieldValueArrayFromNodesInCellZone(
                                                  &x_arr_full,
                                                  get_x_coord,
                                                  &arr_control_size,
                                                  cell_zone
                                                );
  hostGetOrderedFieldValueArrayFromNodesInCellZone(
                                                  &y_arr_full,
                                                  get_y_coord,
                                                  &arr_control_size,
                                                  cell_zone
                                                );
  hostGetOrderedFieldValueArrayFromNodesInCellZone(
                                                  &z_arr_full,
                                                  get_z_coord,
                                                  &arr_control_size,
                                                  cell_zone
                                                  );                  

  #if RP_HOST
  if( arr_control_size == length_arrs_full)
  {
    (*coord_arr_full) = (real (*)[ND_ND]) calloc(ND_ND * length_arrs_full, sizeof(real));

    if ((*coord_arr_full) != NULL)
    {
      for(i = 0; i<length_arrs_full;++i)
      {
        (*coord_arr_full)[i][0] = x_arr_full[i];
        (*coord_arr_full)[i][1] = y_arr_full[i];
        (*coord_arr_full)[i][2] = z_arr_full[i];
      }
    }
    else
    {
      Message("Error hostGetCellCoordsFromNodesInCellZone(): Allocating Memory Problem!\n");
    }
  }
  else
  {
    Message("Error hostGetCellCoordsFromNodesInCellZone(): Arrays have different size!\n");
    Message("length_arrs_full: %i\narr_control_size: %i\n",length_arrs_full, arr_control_size);
  }


  
  #endif

  if(x_arr_full!= NULL)
  {
    free(x_arr_full);
  }

  if(y_arr_full!= NULL)
  {
    free(y_arr_full);
  }

  if(z_arr_full!= NULL)
  {
    free(z_arr_full);
  }

}


void hostGetOrderingArraysFromNodesInCellZone(
                                                int **cid_arr_full, 
                                                int **myid_arr_full,
                                                int *length_arrs_full,
                                                int cell_zone
                                              )
{
/*
  Function to get information from vof coordinates and node distribution info 
  to host.
  IF PARTIONING HAS CHANGED RUN THIS FUNCTION AGAIN!

  TODO: add error states
*/

int sum_size_full = 0;
int i,j = 0;
int size = 0; 

*cid_arr_full = NULL;
*myid_arr_full = NULL;

int compute_node_id;
int *cell_id_arr_node = NULL;
int pe;

#if RP_HOST
int sum_size_nodes = 0;
Message("Receiving ordering of field values for coupling operations...\n");
#endif

#if !RP_HOST 
cell_t c;
Thread *t;
Domain *domain = Get_Domain(1);
Thread **pt;
real xc[ND_ND];

t = Lookup_Thread(domain, cell_zone);


#if RP_NODE /* in !RP_HOST*/
size = THREAD_N_ELEMENTS_INT(t);
cell_id_arr_node = (int *) calloc(size, sizeof(int));

i = 0;
begin_c_loop_int(c, t) 
{
  pt = THREAD_SUB_THREADS(t);
  cell_id_arr_node[i] = (int) c;
  ++i;
}
end_c_loop_int(c, t)

compute_node_id = myid;

/* Set pe to destination node */
/* If on node_0 send data to host */
/* Else send to node_0 because */
/*  compute nodes connect to node_0 & node_0 to host */
/* Order of sending and receiving is very important!*/

PRF_GSYNC();
sum_size_full = PRF_GISUM1(size);

if (I_AM_NODE_ZERO_P)
{  
 PRF_CSEND_INT(node_host, &sum_size_full, 1, myid);
}


pe = (I_AM_NODE_ZERO_P) ? node_host : node_zero;
/*Sent data from nodes to node0 or from node0 to host*/
PRF_CSEND_INT(pe, &size, 1, myid);
PRF_CSEND_INT(pe, &compute_node_id, 1, myid);
PRF_CSEND_INT(pe, cell_id_arr_node, size, myid);

/* free array on nodes once data sent */
free(cell_id_arr_node);

/* node_0 now collect data sent by other compute nodes */
/*  and sends it to the host */
if (I_AM_NODE_ZERO_P)
{  
 /* pe only acts as a counter in this loop */
 compute_node_loop_not_zero (pe) 
 {
   PRF_CRECV_INT(pe, &size, 1, pe);
   
   cell_id_arr_node = (int *) calloc(size, sizeof(int));

   /* Receive data */
   PRF_CRECV_INT(pe, &compute_node_id, 1, pe);
   PRF_CRECV_INT(pe, cell_id_arr_node, size, pe);

   /* send data */
   PRF_CSEND_INT(node_host, &size, 1, myid);
   PRF_CSEND_INT(node_host, &compute_node_id, 1, myid);
   PRF_CSEND_INT(node_host, cell_id_arr_node, size, myid);

   free((char *)cell_id_arr_node);
 }
}
#endif /* RP_NODE in !RP_HOST */
#endif


#if RP_HOST
PRF_CRECV_INT(node_zero, &sum_size_full, 1, node_zero);

if(sum_size_full > 0)
{  
 *cid_arr_full = (int *) calloc(sum_size_full, sizeof(int));
 *myid_arr_full = (int *) calloc(sum_size_full, sizeof(int));

 if ( *cid_arr_full == NULL || *myid_arr_full == NULL)
 {
   Message("Error at Memory Allocation of Variables in host_getDebugCoordinates!");
   
   free(*cid_arr_full);
   free(*myid_arr_full);
 }
 else
 {
   sum_size_nodes = 0;
   /* pe only acts as a counter in this loop */
   compute_node_loop (pe) 
   { 
     PRF_CRECV_INT(node_zero, &size, 1, node_zero);

     cell_id_arr_node = (int *) calloc(size, sizeof(int));


     /* Receive data */
     PRF_CRECV_INT(node_zero, &compute_node_id, 1, node_zero);
     PRF_CRECV_INT(node_zero, cell_id_arr_node, size, node_zero);

     sum_size_nodes += size;

     for(i=(sum_size_nodes-size); i<sum_size_nodes; ++i)
     {	
       if (i >= sum_size_full)
       {
         Message("Warning: Possible memory overflow error in function blocked!\n");
       }
       else
       {

         if ((i - (sum_size_nodes-size)) >= size)
         {
           Message("Warning: Possible pointer out of memory error in function blocked!\n");
         }
         else
         {
           (*cid_arr_full)[i] = cell_id_arr_node[i - (sum_size_nodes-size)];
           (*myid_arr_full)[i] = compute_node_id;
         }
       }
     }
     free(cell_id_arr_node);/* free array on nodes once data sent */
   }
 }
}
(*length_arrs_full) = sum_size_full;

#endif /* RP_HOST */
}


void hostGetOrderedFieldValueArrayFromNodesInCellZone(
                                                      real **val_arr_full,
                                                      real (*C_VAL_WRAPPER_FUN)(cell_t, Thread*),
                                                      int *length_arrs_full, 
                                                      int cell_zone
                                                    )
{
/*
  Function to get information from vof coordinates and node distribution info 
  to host.
  IF PARTIONING HAS CHANGED RUN THIS FUNCTION AGAIN!

  TODO: add error states
*/

int sum_size_full = 0;
int i = 0;
int size = 0; 
int pe;

real *val_arr_node = NULL;

*val_arr_full = NULL;

#if RP_HOST
int sum_size_nodes = 0;
Message("Receiving FLUENT cell field for coupling operations...\n");
#endif

#if !RP_HOST 
cell_t c;
Thread *t;
Domain *domain = Get_Domain(1);
Thread **pt;
real xc[ND_ND];

t = Lookup_Thread(domain, cell_zone);


#if RP_NODE /* in !RP_HOST*/
/* Each Node loads up its data passing array */
size = THREAD_N_ELEMENTS_INT(t);
val_arr_node = (real *) calloc(size, sizeof(real));

i = 0;
begin_c_loop_int(c, t) 
{
  val_arr_node[i] = (*C_VAL_WRAPPER_FUN)(c,t);
  ++i;
}
end_c_loop_int(c, t)

/* Set pe to destination node */
/* If on node_0 send data to host */
/* Else send to node_0 because */
/*  compute nodes connect to node_0 & node_0 to host */
/* Order of sending and receiving is very important!*/

PRF_GSYNC();
sum_size_full = PRF_GISUM1(size);

if (I_AM_NODE_ZERO_P)
{  
 PRF_CSEND_INT(node_host, &sum_size_full, 1, myid);
}


pe = (I_AM_NODE_ZERO_P) ? node_host : node_zero;
/*Sent data from nodes to node0 or from node0 to host*/
PRF_CSEND_INT(pe, &size, 1, myid);
PRF_CSEND_REAL(pe, val_arr_node, size, myid);

/* free array on nodes once data sent */
free(val_arr_node);

/* node_0 now collect data sent by other compute nodes */
/*  and sends it to the host */
if (I_AM_NODE_ZERO_P)
{  
 /* pe only acts as a counter in this loop */
 compute_node_loop_not_zero (pe) 
 {
   PRF_CRECV_INT(pe, &size, 1, pe);
   val_arr_node = (real *) calloc(size, sizeof(real));

   /* Receive data */
   PRF_CRECV_REAL(pe, val_arr_node, size, pe);

   /* send data */
   PRF_CSEND_INT(node_host, &size, 1, myid);
   PRF_CSEND_REAL(node_host, val_arr_node, size, myid);

   free((char *)val_arr_node);
 }
}
#endif /* RP_NODE in !RP_HOST */
#endif

#if RP_HOST
PRF_CRECV_INT(node_zero, &sum_size_full, 1, node_zero);

if(sum_size_full > 0)
{  
 *val_arr_full = (real *) calloc(sum_size_full, sizeof(real));

 if (  *val_arr_full == NULL)
 {
   Message("Error at Memory Allocation of Variables in host_getDebugCoordinates!");
   free(*val_arr_full);
 }
 else
 {
   sum_size_nodes = 0;
   /* pe only acts as a counter in this loop */
   compute_node_loop (pe) 
   { 
     PRF_CRECV_INT(node_zero, &size, 1, node_zero);
     val_arr_node = (real *) calloc(size, sizeof(real));

     /* Receive data */
     PRF_CRECV_REAL(node_zero, val_arr_node, size, node_zero);

     sum_size_nodes += size;

     for(i=(sum_size_nodes-size); i<sum_size_nodes; ++i)
     {	
       if (i >= sum_size_full)
       {
         Message("Warning: Possible memory overflow error in function blocked!\n");
       }
       else
       {

         if ((i - (sum_size_nodes-size)) >= size)
         {
           Message("Warning: Possible pointer out of memory error in function blocked!\n");
         }
         else
         {
           (*val_arr_full)[i] = val_arr_node[i - (sum_size_nodes-size)];
         }
       }
     }

     free(val_arr_node);/* free array on nodes once data sent */
   }
   (*length_arrs_full) = sum_size_nodes;
 }
}
#endif /* RP_HOST */
}
