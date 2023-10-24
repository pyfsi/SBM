#include "udf.h"
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdbool.h>

/* dynamic memory allocation for 1D and 2D arrays */
#define DECLARE_MEMORY(name, type) type *name = NULL

#define DECLARE_MEMORY_N(name, type, dim) type *name[dim] = {NULL}

#define RELEASE_MEMORY(name)                                        \
if (NNULLP(name)) {                                                 \
    free(name);                                                     \
    name = NULL;                                                    \
}

#define RELEASE_MEMORY_N(name, dim)                                 \
for (_d = 0; _d < dim; _d++) {                                      \
    RELEASE_MEMORY(name[_d]);                                       \
}

#define ASSIGN_MEMORY(name, size, type)                             \
if (size) {                                                         \
    if (NNULLP(name)) {                                             \
        name = (type *)realloc(name, size * sizeof(type));          \
    } else {                                                        \
        name = (type *)malloc(size * sizeof(type));                 \
    }                                                               \
    if (NULLP(name)) {                                              \
        Error("\nUDF-error: Memory assignment failed for name.");   \
        exit(1);                                                    \
    }                                                               \
}

#define ASSIGN_MEMORY_N(name, size, type, dim)                      \
for (_d = 0; _d < dim; _d++) {                                      \
    ASSIGN_MEMORY(name[_d], size, type);                            \
}

/* sending and receiving arrays in parallel */
#define PRF_CSEND_INT_N(to, name, n, tag, dim)                      \
for (_d = 0; _d < dim; _d++) {                                      \
    PRF_CSEND_INT(to, name[_d], n, tag);                            \
}

#define PRF_CSEND_REAL_N(to, name, n, tag, dim)                     \
for (_d = 0; _d < dim; _d++) {                                      \
    PRF_CSEND_REAL(to, name[_d], n, tag);                           \
}

#define PRF_CRECV_INT_N(from, name, n, tag, dim)                    \
for (_d = 0; _d < dim; _d++) {                                      \
    PRF_CRECV_INT(from, name[_d], n, tag);                          \
}

#define PRF_CRECV_REAL_N(from, name, n, tag, dim)                   \
for (_d = 0; _d < dim; _d++) {                                      \
    PRF_CRECV_REAL(from, name[_d], n, tag);                         \
}

/* global variables */
#define mnpf 4 /*max nodes per face hardcoded */
#define n_slices 20000
int _d; /* don't use in UDFs! */
int n_threads;
int n_faces;
DECLARE_MEMORY(thread_ids, int);
DECLARE_MEMORY_N(ids, int, mnpf);
DECLARE_MEMORY_N(vof_w, int, n_slices);
int iteration = 0;
int timestep = 0;


  /*----------------------*/
 /* get_inlet_thread_ids */
/*----------------------*/

DEFINE_ON_DEMAND(get_inlet_thread_ids) {
    /* read in thread thread ids, should be called early on */
    if (myid == 0) {printf("\n\nStarted UDF get_inlet_thread_ids.\n"); fflush(stdout);}

#if !RP_NODE
    int k;
    FILE *file;
    file = fopen("inlets.txt", "r");
    fscanf(file, "%i", &n_threads);
#endif /* !RP_NODE */

    host_to_node_int_1(n_threads);
    ASSIGN_MEMORY(thread_ids, n_threads, int);

#if !RP_NODE
    for (k = 0; k < n_threads; k++) {
        fscanf(file, "%*s %i", &thread_ids[k]);
    }
    fclose(file);
#endif /* !RP_NODE */

    host_to_node_int(thread_ids, n_threads);

    if (myid == 0) {printf("\nFinished UDF get_inlet_thread_ids.\n"); fflush(stdout);}

}


  /*-------------------------*/
 /* store_faces_normals_ids */
/*-------------------------*/

DEFINE_ON_DEMAND(store_faces_normals_ids) {
    if (myid == 0) {printf("\n\nStarted UDF store_faces_normals_id.\n"); fflush(stdout);}

    int thread, i_f, d, compute_node;
    DECLARE_MEMORY_N(face_coords, real, ND_ND);
	DECLARE_MEMORY_N(face_areas, real, ND_ND);
    DECLARE_MEMORY_N(face_ids, int, mnpf);

#if !RP_HOST
    Domain *domain;
    Thread *face_thread;
	Node *node;
    face_t face;
    int node_number, j;
    real centroid[ND_ND];
	real area[ND_ND];
#endif /* !RP_HOST */

#if !RP_NODE
    char file_faces_name[ ] = "faces.dat";
    FILE *file_faces = NULL;
#endif /* !RP_NODE */


    for (thread=0; thread<n_threads; thread++) {

#if !RP_NODE
        if (NULLP(file_faces = fopen(file_faces_name, "w"))) {
            Error("\nUDF-error: Unable to open %s for writing\n", file_faces_name);
            exit(1);
        }

#if RP_2D
        fprintf(file_faces, "%27s %27s %27s %27s      %10s\n", "x-coordinate", "y-coordinate", "x-area", "y-area", "unique-ids");
#else /* RP_2D */
        fprintf(file_faces, "%27s %27s %27s %27s %27s %27s      %10s\n", "x-coordinate", "y-coordinate", "z-coordinate", "x-area", "y-area", "z-area",  "unique-ids");
#endif /* RP_2D */
#endif /* !RP_NODE */

#if !RP_HOST
        domain = Get_Domain(1);
        face_thread = Lookup_Thread(domain, thread_ids[thread]);
        n_faces = THREAD_N_ELEMENTS_INT(face_thread);

        ASSIGN_MEMORY_N(face_coords, n_faces, real, ND_ND);
        ASSIGN_MEMORY_N(face_ids, n_faces, int, mnpf);
		ASSIGN_MEMORY_N(face_areas, n_faces, real, ND_ND);

        i_f = 0;
        begin_f_loop(face, face_thread) {
            if (i_f >= n_faces) {Error("\nIndex %i >= array size %i.", i_f, n_faces);}

            F_CENTROID(centroid, face, face_thread); /*centroid of face returned from F_CENTROID */
			F_AREA(area,face,face_thread); /* face normal vector returned from F_AREA */

            for (d = 0; d < ND_ND; d++) {
                face_coords[d][i_f] = centroid[d];
				face_areas[d][i_f] = area[d];
            }

            for (j = 0; j < mnpf; j++) {
                face_ids[j][i_f] = -1;
            }

            j = 0;
            f_node_loop(face, face_thread, node_number) {
                node = F_NODE(face, face_thread, node_number);

                if (j >= mnpf) {Error("\nIndex %i >= array size %i.", j, mnpf);}
                face_ids[j][i_f] = NODE_DM_ID(node);
                j++;
            }
            i_f++;
        } end_f_loop(face, face_thread);
#endif /* !RP_HOST */

#if RP_NODE
        compute_node = (I_AM_NODE_ZERO_P) ? node_host : node_zero;

        PRF_CSEND_INT(compute_node, &n_faces, 1, myid);

        PRF_CSEND_REAL_N(compute_node, face_coords, n_faces, myid, ND_ND);
        PRF_CSEND_REAL_N(compute_node, face_areas, n_faces, myid, ND_ND);
        PRF_CSEND_INT_N(compute_node, face_ids, n_faces, myid, mnpf);

        RELEASE_MEMORY_N(face_coords, ND_ND);
        RELEASE_MEMORY_N(face_areas, ND_ND);
        RELEASE_MEMORY_N(face_ids, mnpf);

        if(I_AM_NODE_ZERO_P){
            compute_node_loop_not_zero(compute_node) {
                PRF_CRECV_INT(compute_node, &n_faces, 1, compute_node);

                ASSIGN_MEMORY_N(face_coords, n_faces, real, ND_ND);
                ASSIGN_MEMORY_N(face_areas, n_faces, real, ND_ND);
                ASSIGN_MEMORY_N(face_ids, n_faces, int, mnpf);

                PRF_CRECV_REAL_N(compute_node, face_coords, n_faces, compute_node, ND_ND);
                PRF_CRECV_REAL_N(compute_node, face_areas, n_faces, compute_node, ND_ND);
                PRF_CRECV_INT_N(compute_node, face_ids, n_faces, compute_node, mnpf);

                PRF_CSEND_INT(node_host, &n_faces, 1, compute_node);

                PRF_CSEND_REAL_N(node_host, face_coords, n_faces, compute_node, ND_ND);
                PRF_CSEND_REAL_N(node_host, face_areas, n_faces, compute_node, ND_ND);
                PRF_CSEND_INT_N(node_host, face_ids, n_faces, compute_node, mnpf);

                RELEASE_MEMORY_N(face_coords, ND_ND);
                RELEASE_MEMORY_N(face_areas, ND_ND);
                RELEASE_MEMORY_N(face_ids, mnpf);
            }
        }
#endif /* RP_NODE */

#if RP_HOST
        compute_node_loop(compute_node) {
            PRF_CRECV_INT(node_zero, &n_faces, 1, compute_node);

            ASSIGN_MEMORY_N(face_coords, n_faces, real, ND_ND);
            ASSIGN_MEMORY_N(face_areas, n_faces, real, ND_ND);
            ASSIGN_MEMORY_N(face_ids, n_faces, int, mnpf);

            PRF_CRECV_REAL_N(node_zero, face_coords, n_faces, compute_node, ND_ND);
            PRF_CRECV_REAL_N(node_zero, face_areas, n_faces, compute_node, ND_ND);
            PRF_CRECV_INT_N(node_zero, face_ids, n_faces, compute_node, mnpf);
#endif /* RP_HOST */

#if !RP_NODE
            for (i_f = 0; i_f < n_faces; i_f++) {
                for (d = 0; d < ND_ND; d++) {
                    fprintf(file_faces, "%27.17e ", face_coords[d][i_f]);
				}
                for (d = 0; d < ND_ND; d++) {
                    fprintf(file_faces, "%27.17e ", face_areas[d][i_f]);
                }
                for (d = 0; d < mnpf; d++) {
                    fprintf(file_faces, " %10d", face_ids[d][i_f]);
                }
                fprintf(file_faces, "\n");
            }

            RELEASE_MEMORY_N(face_coords, ND_ND);
            RELEASE_MEMORY_N(face_areas, ND_ND);
            RELEASE_MEMORY_N(face_ids, mnpf);
#endif /* !RP_NODE */

#if RP_HOST
        } /* close compute_node_loop */
#endif /* RP_HOST */

#if !RP_NODE
        fclose(file_faces);
#endif /* !RP_NODE */

    } /* close loop over threads */

    if (myid == 0) {printf("\nFinished UDF store_faces_normals_id.\n"); fflush(stdout);}
}

  /*-------------------*/
 /* read_vof_face_ids */
/*-------------------*/

DEFINE_ON_DEMAND(read_vof_face_ids) {

#if RP_NODE
int compute_node;
#endif /*RP_NODE*/

#if !RP_NODE
    int i, j;
    char file_ids[] = "face_ids.dat";
    char file_vof[] = "inlet_VOFw.dat";
    FILE *fp_ids;
    FILE *fp_vof;

    if (NULLP(fp_ids = fopen(file_ids, "r"))) {
        Error("\nUDF-error: Unable to open %s for reading\n", file_ids);
        exit(1);
    }
    fscanf(fp_ids, "# %i", &n_faces);

    if (NULLP(fp_vof = fopen(file_vof, "r"))) {
        Error("\nUDF-error: Unable to open %s for reading\n", file_vof);
        exit(1);
    }

    ASSIGN_MEMORY_N(ids, n_faces, int, mnpf);
    ASSIGN_MEMORY_N(vof_w, n_faces, int, n_slices);

    for (i=0; i<n_faces; i++){
        for (j=0; j<mnpf; j++){
            fscanf(fp_ids, "%i", &ids[j][i]);
            printf("%i ", ids[j][i]);
        }
        printf("\n");
        for (j=0; j<n_slices; j++){
            fscanf(fp_vof, "%i", &vof_w[j][i]);
        }
    }
    fclose(fp_ids);
    fclose(fp_vof);

#endif /*!RP_NODE*/

    /* Send data to nodes */
    host_to_node_int_1(n_faces);

#if RP_HOST
    PRF_CSEND_INT_N(node_zero, ids, n_faces, myid, mnpf); /* send from host to node0 */
    PRF_CSEND_INT_N(node_zero, vof_w, n_faces, myid, mnpf); /* send from host to node0 */
#endif /*RP_HOST*/

#if RP_NODE
    ASSIGN_MEMORY_N(ids, n_faces, int, mnpf);
    ASSIGN_MEMORY_N(vof_w, n_faces, int, n_slices);

    if(I_AM_NODE_ZERO_P){
        PRF_CRECV_INT_N(node_host, ids, n_faces, node_host, mnpf);
        PRF_CRECV_INT_N(node_host, vof_w, n_faces, node_host, mnpf);
        compute_node_loop_not_zero(compute_node){
        PRF_CSEND_INT_N(compute_node, ids, n_faces, myid, mnpf);
        PRF_CSEND_INT_N(compute_node, vof_w, n_faces, myid, mnpf);
        }
    }
    else {
        PRF_CRECV_INT_N(node_zero, ids, n_faces, node_zero, mnpf);
        PRF_CRECV_INT_N(node_zero, vof_w, n_faces, node_zero, mnpf);
    }
#endif /*RP_NODE*/

    Message0("\nFinished UDF read_vof_faces_ids.\n");
}

DEFINE_ON_DEMAND(test){

#if !RP_HOST
    int i_f, i_n, j, thread, node_number;
    int i_f_min = 0;
    int match_counter =0 ;
    int loop_counter = 0;
    bool match_found, all_matched;
    Domain *domain = Get_Domain(1);
    Thread *face_thread;
    face_t face;
    Node *node;
    real centroid[ND_ND];
	real area[ND_ND];
    real x[ND_ND]; /* this will hold the position vector */
    int node_ids[mnpf];

    for (thread=0; thread<n_threads; thread++) { /* not in udf_profile because there you get the thread*/
        face_thread = Lookup_Thread(domain, thread_ids[thread]);
        begin_f_loop(face,face_thread)
        {
            /* store node ids in array*/
            for (i_n=0; i_n<mnpf; i_n++)
                node_ids[i_n] = -1; /* if not all cells same number of nodes, fill with -1*/
            i_n = 0;

            f_node_loop(face, face_thread, node_number) {
                node = F_NODE(face, face_thread, node_number);
                node_ids[i_n++] = NODE_DM_ID(node); /* store ID and post-increment iterator*/
                /* Message("Node ID: %i\n", node_ids[i_n-1]); */
            }
            /*compare ids to file to find correct line (assumed different order) */
            for (i_f=i_f_min; i_f<n_faces; i_f++) {
                loop_counter++;
                all_matched = true;
                i_n = 0;
                while (all_matched && (i_n < mnpf)) {
                    match_found = false;
                    j = 0;
                    while (!match_found && (j < mnpf)) {  /* compare them without assuming ordering */
                        match_found = (ids[i_n][i_f] == node_ids[j++]); /* compare and post-increment j */
                    }
                    if (!match_found) {
                        all_matched = false;
                    }
                    else {
                    i_n++;
                    }
                }
                if (all_matched) {
                    match_counter++;
                    break;  /* current i_f will be number in file */
                }
            }
            if (i_f == i_f_min + 1) {
                i_f_min++; /* for efficiency */
            }
            /* do SBM profile stuff */
        }
        end_f_loop(f,t)
        Message("\nTotal matches: %i in a total of %i loop executions\n", match_counter, loop_counter);
    }

#endif /*!RP_HOST*/
}