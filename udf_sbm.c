#include "udf.h"
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

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
int _d; /* don't use in UDFs! */
int n_threads;
DECLARE_MEMORY(thread_ids, int);
int iteration = 0;
int timestep = 0;


  /*----------------*/
 /* get_inlet_thread_ids */
/*----------------*/

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

    int thread, n_faces, i_f, d, compute_node;
    DECLARE_MEMORY_N(face_coords, real, ND_ND);
    DECLARE_MEMORY_N(face_ids, int, mnpf);
	DECLARE_MEMORY_N(face_areas, real, ND_ND);

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

            F_CENTROID(centroid, face, face_thread);
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

