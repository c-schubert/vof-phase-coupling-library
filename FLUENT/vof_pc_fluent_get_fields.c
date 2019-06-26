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



void hostGetCouplingFieldsFromNodesinCellZone(
                                                real (**coord_arr_full)[ND_ND],
                                                real **vof_arr_full, 
                                                int **cell_id_arr_full, 
                                                int **compute_node_id_arr_full,
                                                int *arr_full_size,
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

*vof_arr_full = NULL;
*cell_id_arr_full = NULL;
*coord_arr_full = NULL;
*compute_node_id_arr_full = NULL;

real (*coord_arr_node)[ND_ND] = NULL;
real *vof_arr_node = NULL;
int compute_node_id;
int *cell_id_arr_node = NULL;
int pe;

#if RP_HOST
int sum_size_nodes = 0;
Message("Receiving FLUENT cell coordinates and VOF field for coupling operations...\n");
#endif


#if !RP_HOST 
cell_t c;
Thread *t;
Domain *domain = Get_Domain(1); // mixture domain
Thread **pt;
real xc[ND_ND];

t = Lookup_Thread(domain, cell_zone);


#if RP_NODE /* in !RP_HOST*/
/* Each Node loads up its data passing array */

/*
On a compute-node, THREAD_N_ELEMENTS(t) contains the combined
number of elements in the compute-node's interior region and
all exterior layers of the compute-node

THREAD_N_ELEMENTS_INT(t) - only the interior number of elements

THREAD_N_ELEMENTS_EXT(t) - only the exterior number of elements

see Exterior Thread Storage in Parallel Considerations in Fluent Customization Guide
*/
size = THREAD_N_ELEMENTS_INT(t);

vof_arr_node = (real *) calloc(size, sizeof(real));
cell_id_arr_node = (int *) calloc(size, sizeof(int));
coord_arr_node = (real (*)[ND_ND]) calloc(ND_ND * size, sizeof(real));

i = 0;
begin_c_loop_int(c, t) 
{
  pt = THREAD_SUB_THREADS(t);
  vof_arr_node[i] = C_VOF(c, pt[0]);
  cell_id_arr_node[i] = (int) c;
  C_CENTROID(coord_arr_node[i], c, t);
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
PRF_CSEND_REAL(pe, vof_arr_node, size, myid);
PRF_CSEND_REAL(pe, coord_arr_node[0], ND_ND * size, myid);

/* free array on nodes once data sent */
free(vof_arr_node);
free(cell_id_arr_node);
free(coord_arr_node);

/* node_0 now collect data sent by other compute nodes */
/*  and sends it to the host */
if (I_AM_NODE_ZERO_P)
{  
 /* pe only acts as a counter in this loop */
 compute_node_loop_not_zero (pe) 
 {
   PRF_CRECV_INT(pe, &size, 1, pe);
   
   cell_id_arr_node = (int *) calloc(size, sizeof(int));
   vof_arr_node = (real *) calloc(size, sizeof(real));
   coord_arr_node = (real (*)[ND_ND]) calloc(ND_ND * size, sizeof(real));

   /* Receive data */
   PRF_CRECV_INT(pe, &compute_node_id, 1, pe);
   PRF_CRECV_INT(pe, cell_id_arr_node, size, pe);
   PRF_CRECV_REAL(pe, vof_arr_node, size, pe);
   PRF_CRECV_REAL(pe, coord_arr_node[0], ND_ND * size, pe);

   /* send data */
   PRF_CSEND_INT(node_host, &size, 1, myid);
   PRF_CSEND_INT(node_host, &compute_node_id, 1, myid);
   PRF_CSEND_INT(node_host, cell_id_arr_node, size, myid);
   PRF_CSEND_REAL(node_host, vof_arr_node, size, myid);
   PRF_CSEND_REAL(node_host, coord_arr_node[0], ND_ND * size, myid);

   free((char *)vof_arr_node);
   free((char *)cell_id_arr_node);
   free((char *)coord_arr_node);
 }
}
#endif /* RP_NODE in !RP_HOST */
#endif



#if RP_HOST
PRF_CRECV_INT(node_zero, &sum_size_full, 1, node_zero);
// Message("sum_size_full: %i\n", sum_size_full);

if(sum_size_full > 0)
{  
 *vof_arr_full = (real *) calloc(sum_size_full, sizeof(real));
 *cell_id_arr_full = (int *) calloc(sum_size_full, sizeof(int));
 *coord_arr_full = (real (*)[ND_ND]) calloc(ND_ND * sum_size_full, sizeof(real));
 *compute_node_id_arr_full = (int *) calloc(sum_size_full, sizeof(int));

 if (  *vof_arr_full == NULL 
   || *cell_id_arr_full == NULL 
   || *coord_arr_full == NULL 
   || *compute_node_id_arr_full == NULL
   )
 {
   Message("Error at Memory Allocation of Variables in host_getDebugCoordinates!");
   
   free(*coord_arr_full);
   free(*vof_arr_full);
   free(*cell_id_arr_full);
   free(*compute_node_id_arr_full);
 }
 else
 {
   sum_size_nodes = 0;
   /* pe only acts as a counter in this loop */
   compute_node_loop (pe) 
   { 
     PRF_CRECV_INT(node_zero, &size, 1, node_zero);

     vof_arr_node = (real *) calloc(size, sizeof(real));
     cell_id_arr_node = (int *) calloc(size, sizeof(int));
     coord_arr_node = (real (*)[ND_ND]) calloc(ND_ND * size, sizeof(real));

     /* Receive data */
     PRF_CRECV_INT(node_zero, &compute_node_id, 1, node_zero);
     PRF_CRECV_INT(node_zero, cell_id_arr_node, size, node_zero);
     PRF_CRECV_REAL(node_zero, vof_arr_node, size, node_zero);
     PRF_CRECV_REAL(node_zero, coord_arr_node[0], ND_ND * size, node_zero);

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
           (*vof_arr_full)[i] = vof_arr_node[i - (sum_size_nodes-size)];
           (*cell_id_arr_full)[i] = cell_id_arr_node[i - (sum_size_nodes-size)];

           for(j = 0; j < ND_ND; ++j)
           {
             (*coord_arr_full)[i][j] = coord_arr_node[i - (sum_size_nodes-size)][j];
           }

           (*compute_node_id_arr_full)[i] = compute_node_id;
         }
       }
     }

     free(vof_arr_node);/* free array on nodes once data sent */
     free(cell_id_arr_node);/* free array on nodes once data sent */
     free(coord_arr_node);/* free array on nodes once data sent */
   }
 }
}
*arr_full_size = sum_size_full;

#endif /* RP_HOST */
}
