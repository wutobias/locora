#define _USE_MATH_DEFINES
#include <math.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/arrayobject.h"
#include <stdio.h>
#include <string.h>
#include "Vec.h"

PyObject *_crd_systems_ext_axis_paral(PyObject *self, PyObject *args) {

    int n_crds;
    int crds_i, crds_j;            // Coordinate iterator
    double vec_q1[3], vec_q2[3];   // Coordinate vectors for query coordinates
    double vec_q1_avg[3];          // Average orientation of vector vec_q1
    vec_q1_avg[0] = 0.;
    vec_q1_avg[1] = 0.;
    vec_q1_avg[2] = 0.;

    PyArrayObject *query_coords_array, *target_coords_array; // Numpy array objects of query coordinates and
                                                             // and target coordinates.

    PyArrayObject *ext_ref_coords_array; // external references coordinates
    double ext_ref_vec[3];      // external reference vector
    int has_ref = 0; // bool that tells us if should use external reference or not

    double ref[3];         // A reference vector
    double diff_vec[3];    // Difference vector for construction of connection vector

    double paral_vec[3];   // The final vector. Will be copied to target_coords_array in the end
    paral_vec[0] = 0.;
    paral_vec[1] = 0.;
    paral_vec[2] = 0.;

    int n_comb = 0;  // Total number of different combinations considered in calculation.

    int verbose = 0; // verboes mode on/off
    
    //Here we retrieve the data from python
    if (!PyArg_ParseTuple(args, "O!O!i|O!i",
        &PyArray_Type, &query_coords_array,
        &PyArray_Type, &target_coords_array,
        &has_ref,
        &PyArray_Type, &ext_ref_coords_array,
        &verbose
        ))
    {
        return NULL;
    }

    n_crds = PyArray_DIM(query_coords_array, 0);

    for (crds_i=0; crds_i<n_crds; crds_i++){

        vec_q1[0] = *(float *) PyArray_GETPTR2(query_coords_array, crds_i, 0);
        vec_q1[1] = *(float *) PyArray_GETPTR2(query_coords_array, crds_i, 1);
        vec_q1[2] = *(float *) PyArray_GETPTR2(query_coords_array, crds_i, 2);

        for (crds_j=crds_i+1; crds_j<n_crds; crds_j++){

            if (verbose){
                printf("Iteration (%d %d)\n", crds_i, crds_j);
            }

            vec_q2[0] = *(float *) PyArray_GETPTR2(query_coords_array, crds_j, 0);
            vec_q2[1] = *(float *) PyArray_GETPTR2(query_coords_array, crds_j, 1);
            vec_q2[2] = *(float *) PyArray_GETPTR2(query_coords_array, crds_j, 2);

            vec_q1_avg[0] += vec_q1[0];
            vec_q1_avg[1] += vec_q1[1];
            vec_q1_avg[2] += vec_q1[2];

            diff_vec[0] = vec_q2[0] - vec_q1[0];
            diff_vec[1] = vec_q2[1] - vec_q1[1];
            diff_vec[2] = vec_q2[2] - vec_q1[2];

            if (verbose){
                printf("Vec1: (%f %f %f)\n", vec_q1[0], vec_q1[1], vec_q1[2]);
                printf("Vec2: (%f %f %f)\n", vec_q2[0], vec_q2[1], vec_q2[2]);
                printf("Diff: (%f %f %f)\n", diff_vec[0], diff_vec[1], diff_vec[2]);
            }

            norm3(diff_vec);

            if (verbose){
                printf("Norm Diff: (%f %f %f)\n", diff_vec[0], diff_vec[1], diff_vec[2]);
            }

            if ( crds_i == 0 && crds_j == 1){
                ref[0] = diff_vec[0];
                ref[1] = diff_vec[1];
                ref[2] = diff_vec[2];

            }

            else{
                if (dotp(ref, diff_vec) < 0.) {

                    diff_vec[0] *= -1.;
                    diff_vec[1] *= -1.;
                    diff_vec[2] *= -1.;
                }
            }

            if (verbose){
                printf("Ref  Diff: (%f %f %f)\n", diff_vec[0], diff_vec[1], diff_vec[2]);
            }

            paral_vec[0] += diff_vec[0];
            paral_vec[1] += diff_vec[1];
            paral_vec[2] += diff_vec[2];

            n_comb++;

        }

    }

    paral_vec[0] /= (double) n_comb;
    paral_vec[1] /= (double) n_comb;
    paral_vec[2] /= (double) n_comb;

    norm3(paral_vec);

    if (has_ref){

        vec_q1_avg[0] /= (double) n_comb;
        vec_q1_avg[1] /= (double) n_comb;
        vec_q1_avg[2] /= (double) n_comb;

        ext_ref_vec[0] = *(float *) PyArray_GETPTR1(ext_ref_coords_array, 0);
        ext_ref_vec[1] = *(float *) PyArray_GETPTR1(ext_ref_coords_array, 1);
        ext_ref_vec[2] = *(float *) PyArray_GETPTR1(ext_ref_coords_array, 2);

        ext_ref_vec[0] = ext_ref_vec[0] - vec_q1_avg[0];
        ext_ref_vec[1] = ext_ref_vec[1] - vec_q1_avg[1];
        ext_ref_vec[2] = ext_ref_vec[2] - vec_q1_avg[2];

        if (dotp(paral_vec, ext_ref_vec) < 0.){
            paral_vec[0] *= -1.;
            paral_vec[1] *= -1.;
            paral_vec[2] *= -1.;
        }

    }

    *(double *)PyArray_GETPTR1(target_coords_array, 0) = paral_vec[0];
    *(double *)PyArray_GETPTR1(target_coords_array, 1) = paral_vec[1];
    *(double *)PyArray_GETPTR1(target_coords_array, 2) = paral_vec[2];

    return Py_BuildValue("i", 1);

}

PyObject *_crd_systems_ext_axis_ortho(PyObject *self, PyObject *args) {

    int n_crds;                            // Number of 3-dim coordinates in query_coords array
    int crds_i, crds_j, crds_k;            // Coordinate iterator
    double vec_q1[3], vec_q2[3], vec_q3[3];// Coordinate vectors for query coordinates
    double vec_q1_avg[3];
    vec_q1_avg[0] = 0.;
    vec_q1_avg[1] = 0.;
    vec_q1_avg[2] = 0.;

    PyArrayObject *query_coords_array, *target_coords_array; // Numpy array objects of query coordinates and
                                                             // and target coordinates.

    PyArrayObject *ext_ref_coords_array; // external references coordinates
    double ext_ref_vec[3];      // external reference vector
    int has_ref = 0;            // bool that tells us if should use external reference or not

    double ref[3];                         // A reference vector
    double plane_vec1[3], plane_vec2[3], plane_vec3[3]; // Vectors for construction of plane

    double ortho_vec[3]; // The final vector. Will be copied to target_coords_array in the end
    ortho_vec[0] = 0.;
    ortho_vec[1] = 0.;
    ortho_vec[2] = 0.;

    int n_comb = 0;   // Total number of different combinations considered in calculation.

    int verbose = 0; // verboes mode on/off
    
    //Here we retrieve the data from python
    if (!PyArg_ParseTuple(args, "O!O!i|O!i",
        &PyArray_Type, &query_coords_array,
        &PyArray_Type, &target_coords_array,
        &has_ref,
        &PyArray_Type, &ext_ref_coords_array,
        &verbose
        ))
    {
        return NULL;
    }

    n_crds = PyArray_DIM(query_coords_array, 0);

    for (crds_i=0; crds_i<n_crds; crds_i++){

            vec_q1[0] = *(float *) PyArray_GETPTR2(query_coords_array, crds_i, 0);
            vec_q1[1] = *(float *) PyArray_GETPTR2(query_coords_array, crds_i, 1);
            vec_q1[2] = *(float *) PyArray_GETPTR2(query_coords_array, crds_i, 2);

            for (crds_j=crds_i+1; crds_j<n_crds; crds_j++){

                vec_q2[0] = *(float *) PyArray_GETPTR2(query_coords_array, crds_j, 0);
                vec_q2[1] = *(float *) PyArray_GETPTR2(query_coords_array, crds_j, 1);
                vec_q2[2] = *(float *) PyArray_GETPTR2(query_coords_array, crds_j, 2);

                for (crds_k=crds_j+1; crds_k<n_crds; crds_k++){

                    if (verbose){
                        printf("Iteration (%d %d %d)\n", crds_i, crds_j, crds_k);
                    }

                    vec_q3[0] = *(float *) PyArray_GETPTR2(query_coords_array, crds_k, 0);
                    vec_q3[1] = *(float *) PyArray_GETPTR2(query_coords_array, crds_k, 1);
                    vec_q3[2] = *(float *) PyArray_GETPTR2(query_coords_array, crds_k, 2);

                    plane_vec1[0] = vec_q3[0] - vec_q1[0];
                    plane_vec1[1] = vec_q3[1] - vec_q1[1];
                    plane_vec1[2] = vec_q3[2] - vec_q1[2];

                    plane_vec2[0] = vec_q2[0] - vec_q1[0];
                    plane_vec2[1] = vec_q2[1] - vec_q1[1];
                    plane_vec2[2] = vec_q2[2] - vec_q1[2];

                    vec_q1_avg[0] += vec_q1[0];
                    vec_q1_avg[1] += vec_q1[1];
                    vec_q1_avg[2] += vec_q1[2];

                    crossp(plane_vec1, plane_vec2, plane_vec3);
                    norm3(plane_vec3);

                    if ( crds_j == 1 && crds_k == 2){
                        ref[0] = plane_vec3[0];
                        ref[1] = plane_vec3[1];
                        ref[2] = plane_vec3[2];
                    }

                    else{
                        if (dotp(ref, plane_vec3) < 0.){

                            plane_vec3[0] *= -1.;
                            plane_vec3[1] *= -1.;
                            plane_vec3[2] *= -1.;
                        }
                    }

                    ortho_vec[0] += plane_vec3[0];
                    ortho_vec[1] += plane_vec3[1];
                    ortho_vec[2] += plane_vec3[2];

                    n_comb++;
                }
        }
    }

    ortho_vec[0] /= (double) n_comb;
    ortho_vec[1] /= (double) n_comb;
    ortho_vec[2] /= (double) n_comb;

    norm3(ortho_vec);

    if (has_ref){

        vec_q1_avg[0] /= (double) n_comb;
        vec_q1_avg[1] /= (double) n_comb;
        vec_q1_avg[2] /= (double) n_comb;

        ext_ref_vec[0] = *(float *) PyArray_GETPTR1(ext_ref_coords_array, 0);
        ext_ref_vec[1] = *(float *) PyArray_GETPTR1(ext_ref_coords_array, 1);
        ext_ref_vec[2] = *(float *) PyArray_GETPTR1(ext_ref_coords_array, 2);

        ext_ref_vec[0] = ext_ref_vec[0] - vec_q1_avg[0];
        ext_ref_vec[1] = ext_ref_vec[1] - vec_q1_avg[1];
        ext_ref_vec[2] = ext_ref_vec[2] - vec_q1_avg[2];

        if (dotp(ortho_vec, ext_ref_vec) < 0.){
            ortho_vec[0] *= -1.;
            ortho_vec[1] *= -1.;
            ortho_vec[2] *= -1.;
        }

    }

    *(double *)PyArray_GETPTR1(target_coords_array, 0) = ortho_vec[0];
    *(double *)PyArray_GETPTR1(target_coords_array, 1) = ortho_vec[1];
    *(double *)PyArray_GETPTR1(target_coords_array, 2) = ortho_vec[2];

    return Py_BuildValue("i", 1);

}


static PyMethodDef _crd_systems_ext_methods[] = {
    {
        "axis_ortho",
        (PyCFunction)_crd_systems_ext_axis_ortho,
        METH_VARARGS,
        "Takes a set of coordinates and calculates\n"
        "the mean normal vector of all possible\n"
        "planes constructed any uniq triple of vectors.",
    },
    {
        "axis_paral",
        (PyCFunction)_crd_systems_ext_axis_paral,
        METH_VARARGS,
        "Takes a set of coordinates and calculates\n"
        "the mean connecting vector between all\n"
        "possible uniq tuple of vectors.",
    },
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
    #define MOD_ERROR_VAL NULL
    #define MOD_SUCCESS_VAL(val) val
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
    #define MOD_DEF(ob, name, doc, methods) \
            static struct PyModuleDef extmoduledef = { \
              PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
            ob = PyModule_Create(&extmoduledef);
#else
    #define MOD_ERROR_VAL
    #define MOD_SUCCESS_VAL(val)
    #define MOD_INIT(name) void init##name(void)
    #define MOD_DEF(ob, name, doc, methods) \
            ob = Py_InitModule3(name, methods, doc);
#endif

/* Initialization function for this module
 */

MOD_INIT(_crd_systems_ext)
{
    PyObject *m;

    MOD_DEF(m, "_crd_systems_ext", "Perform auxilliary calculations for.\n"
                                   "coordinate system Initialization.\n", _crd_systems_ext_methods)
    
    if (m == NULL)
        return MOD_ERROR_VAL;

    import_array();

    return MOD_SUCCESS_VAL(m);
}