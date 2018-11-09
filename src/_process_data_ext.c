#define _USE_MATH_DEFINES
#include <math.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/arrayobject.h"
#include <stdio.h>
#include <string.h>
#include "Vec.h"

static PyObject * 
_lifetime_distr_majmin_ext(PyObject *self, PyObject *args) {

    PyArrayObject *input_array, *lifetime_array; // input_array contains water oxygen atom ids for each frame
                                                 // lifetime_array will store the average lifetimes
    PyArrayObject *frames_array; // frames_array containes the frame idx of each water oxygen atom present in input_array
    PyArrayObject *vec_array;   // vec_array contains orientation vector for each water molecule

    int *wat_id; // oxygen atom idx of water molecule that should be analyzed
    int *wat_0, *wat_k; // water oxygen atom idx at starting point 0 and after lag time k
    int *frame_0, *frame_k; // frame idx at starting point 0 and after lag time k;
    double vec_0[3], vec_k[3]; // orientation vector to analyze at starting point 0 and after lag time k
    double vec_prod; // vector product between vec_0 and vec_k
    int i,j,k; // counting variables
    int lag; // lag time
    int frame_kl; // last frame idx that was populated before current frame frame_0, frame_k
    int tol; // tolerance time for transient escape. Solvent molecules that leave some hydration state shorter than this time,
             // will be treated as if they did not leave the hydration state.
    int order; // order of the lagrange polynom used for orientation calculation
    int length; // length of the input_array

    int verbose=0; //verbose on/off

    //Here we retrieve the data from python
    if (!PyArg_ParseTuple(args, "O!O!O!ii|i",
        &PyArray_Type, &input_array,
        &PyArray_Type, &frames_array,
        &PyArray_Type, &vec_array,
        &tol,
        &order,
        &verbose
        ))
    {
        return NULL;
    }

    length = PyArray_DIM(input_array, 0);
    if (length != PyArray_DIM(frames_array, 0)){
        PyErr_Format(PyExc_ValueError,
                     "Both frames_array and input_array must have same length.\n"
                     );
        return NULL;
    }

    if (order>2){
        if (verbose){
            printf("order>2 is not allowed. Forcing order=2.\n");
        }
        order=2;
    }

    if (order<1){
        if (verbose){
            printf("order<1 is not allowed. Forcing order=1.\n");
        }
        order=1;
    }

    if (verbose){
        printf("Performing order=%d calculation.\n", order);
    }

    // Get uniq water ids
    if (verbose){
        printf("Getting unique water ids...\n");
    }
    int uniq_list[length];
    int next_i_frwrd=0;
    int next_i_bckwrds=length-1;
    int found_wat=0;
    int maxframe=0;
    for (i=0; i<length; i++){
        uniq_list[i] = -9999;
    }
    for (i=0; i<length; i++){
        wat_0 = PyArray_GETPTR1(input_array, i);
        frame_0 = (int *) PyArray_GETPTR1(frames_array, i);
        if (*frame_0>maxframe){
            maxframe = 0;
            maxframe += *frame_0;
        }
        found_wat = 0;
        for (j=0; j<length; j++){
            if (*wat_0 == uniq_list[j]){
                uniq_list[next_i_bckwrds] = -1;
                next_i_bckwrds--;
                found_wat = 1;
                break;
            }
        }
        if (!found_wat){
            uniq_list[next_i_frwrd] = *wat_0;
            next_i_frwrd++;
        }
    }
    if (verbose){
        printf("List of unique water ids: ");
        for (i=0; i<length; i++){
            if (uniq_list[i]>-1){
                printf("%d ", uniq_list[i]);
            }
        }
        printf("\n");
    }

    maxframe += tol;
    npy_intp dims[1];
    dims[0] = maxframe;
    lifetime_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (lifetime_array == NULL) {
        printf("Creating lifetime_array of shape (%d,) failed.\n", length);
        return NULL; /*raises exception*/
    }
    for (i=0; i<maxframe; i++){
        *(double *) PyArray_GETPTR1(lifetime_array, i) = 0.;
    }

    if (verbose){
        printf("Starting lifetime calculation routine...\n");
    }
    // Calculate lifetimes
    for (k=0; k<length; k++){
        wat_id = &uniq_list[k];
        if (*wat_id == -1) break;

        for (i=0; i<length; i++){
            wat_0   = (int *) PyArray_GETPTR1(input_array, i);
            frame_0 = (int *) PyArray_GETPTR1(frames_array, i);
            if (*wat_0 != *wat_id) continue;
            frame_kl = *frame_0;

            vec_0[0] = *(double *) PyArray_GETPTR2(vec_array, i, 0);
            vec_0[1] = *(double *) PyArray_GETPTR2(vec_array, i, 1);
            vec_0[2] = *(double *) PyArray_GETPTR2(vec_array, i, 2);

            for (j=i; j<length; j++){
                wat_k   = (int *) PyArray_GETPTR1(input_array, j);
                frame_k = (int *) PyArray_GETPTR1(frames_array, j);

                if (verbose) {
                    printf("Analyzing frames (%d,%d) and water molecules (%d,%d).\n", *frame_0, *frame_k, *wat_0, *wat_k);
                }

                if (*wat_k != *wat_id) continue;
                if ((*frame_k-frame_kl)>tol+1) break; // frame_k-frame_kl+1: time in between leaving and re-entering (transient escape time)
                lag = *frame_k - *frame_0;          // frame_k-frame_0: actual time difference between frames
                frame_kl = *frame_k;
                
                vec_k[0] = *(double *) PyArray_GETPTR2(vec_array, j, 0);
                vec_k[1] = *(double *) PyArray_GETPTR2(vec_array, j, 1);
                vec_k[2] = *(double *) PyArray_GETPTR2(vec_array, j, 2);

                vec_prod = dotp(vec_0, vec_k);
                if (order==2){
                    *(double *) PyArray_GETPTR1(lifetime_array, lag) += (3.*pow(vec_prod,2)-1)/2.;
                }
                if (order==1){
                    *(double *) PyArray_GETPTR1(lifetime_array, lag) += vec_prod;
                }
                if (verbose) {
                    printf("Set lifetime_array[%d] to %f.\n", lag, *(double *) PyArray_GETPTR1(lifetime_array, lag));
                }
            }
        }
    }

    if (verbose){
        printf("Finished lifetime calculation routine.\n");
    }

    return PyArray_Return(lifetime_array);

}

static PyObject * 
_lifetime_distr_ext(PyObject *self, PyObject *args) {

    PyArrayObject *input_array;    // input_array contains water oxygen atom ids for each frame
    PyArrayObject *lifetime_array; // lifetime_array will store the average lifetimes
    PyArrayObject *frames_array;   // frames_array containes the frame idx of each water oxygen atom present in input_array

    int *wat_id; // oxygen atom idx of water molecule that should be analyzed
    int *wat_0, *wat_k; // water oxygen atom idx at starting point 0 and after lag time k
    int *frame_0, *frame_k; // frame idx at starting point 0 and after lag time k;
    int i,j,k; // counting variables
    int lag; // lag time
    int frame_kl; // last frame idx that was populated before current frame frame_0, frame_k
    int tol; // tolerance time for transient escape. Solvent molecules that leave some hydration state shorter than this time,
             // will be treated as if they did not leave the hydration state.
    int length; // length of the input_array

    int verbose=0; //verbose on/off

    //Here we retrieve the data from python
    if (!PyArg_ParseTuple(args, "O!O!i|i",
        &PyArray_Type, &input_array,
        &PyArray_Type, &frames_array,
        &tol,
        &verbose
        ))
    {
        return NULL;
    }

    length = PyArray_DIM(input_array, 0);
    if (length != PyArray_DIM(frames_array, 0)){
        PyErr_Format(PyExc_ValueError,
                     "Both frames_array and input_array must have same length.\n"
                     );
        return NULL;
    }

    // Get uniq water ids
    if (verbose){
        printf("Getting unique water ids...\n");
    }
    int uniq_list[length];
    int next_i_frwrd=0;
    int next_i_bckwrds=length-1;
    int found_wat=0;
    int maxframe=0;
    for (i=0; i<length; i++){
        uniq_list[i] = -9999;
    }
    for (i=0; i<length; i++){
        wat_id = (int *) PyArray_GETPTR1(input_array, i);
        frame_0 = (int *) PyArray_GETPTR1(frames_array, i);
        if (*frame_0>maxframe){
            maxframe = 0;
            maxframe += *frame_0;
        }
        found_wat = 0;
        for (j=0; j<length; j++){
            if (*wat_id == uniq_list[j]){
                uniq_list[next_i_bckwrds] = -1;
                next_i_bckwrds--;
                found_wat = 1;
                break;
            }
        }
        if (!found_wat){
            uniq_list[next_i_frwrd] = 0;
            uniq_list[next_i_frwrd] += *wat_id;
            next_i_frwrd++;
        }
    }
    if (verbose){
        printf("List of unique water ids: ");
        for (i=0; i<length; i++){
            if (uniq_list[i]>-1){
                printf("%d ", uniq_list[i]);
            }
        }
        printf("\n");
    }

    maxframe += tol;
    npy_intp dims[1];
    dims[0] = maxframe;
    lifetime_array = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_INT);
    if (lifetime_array == NULL) {
        printf("Creating lifetime_array of shape (%d,) failed.\n", length);
        return NULL; /*raises exception*/
    }
    for (i=0; i<maxframe; i++){
        *(int *) PyArray_GETPTR1(lifetime_array, i) = 0;
    }

    if (verbose){
        printf("Starting lifetime calculation routine...\n");
    }

    // Calculate lifetimes
    tol+=1; //Because we want lag that is lower of equal to tol
    for (k=0; k<length; k++){
        wat_id = &uniq_list[k];
        if (*wat_id == -1) break;

        for (i=0; i<length; i++){
            wat_0   = (int *) PyArray_GETPTR1(input_array, i);
            frame_0 = (int *) PyArray_GETPTR1(frames_array, i);
            if (*wat_0 != *wat_id) continue;
            frame_kl = *frame_0;

            for (j=i; j<length; j++){
                wat_k   = (int *) PyArray_GETPTR1(input_array, j);
                frame_k = (int *) PyArray_GETPTR1(frames_array, j);

                if (verbose) {
                    printf("Iteration step (%d,%d): ", i,j);
                    printf("Analyzing frames (%d,%d) and water molecules (%d,%d).\n", *frame_0, *frame_k, *wat_0, *wat_k);
                }

                if (*wat_k != *wat_id) continue;
                if (*frame_k-frame_kl>tol) break; // frame_k-frame_kl+1: time in between leaving and re-entering (transient escape time)
                lag = *frame_k - *frame_0;        // frame_k-frame_0: actual time difference between frames
                frame_kl = *frame_k;
                *(int *)PyArray_GETPTR1(lifetime_array, lag) += 1;
                if (verbose) {
                    printf("Set lifetime_array[%d] to %d.\n", lag, *(int *) PyArray_GETPTR1(lifetime_array, lag));
                }
            }
        }
    }

    if (verbose){
        printf("Finished lifetime calculation routine.\n");
    }

    return PyArray_Return(lifetime_array);

}


static PyMethodDef _utils_ext_methods[] = {
    {
        "lifetime_distr_ext",
        (PyCFunction)_lifetime_distr_ext,
        METH_VARARGS,
        "Calculates lifetime distributions of binary time series.\n"
        "Suitable for analyzing translation diffusion of water molecules.",
    },
    {
        "lifetime_distr_majmin_ext",
        (PyCFunction)_lifetime_distr_majmin_ext,
        METH_VARARGS,
        "Calculate first order/second order lifetime distribution of vector orientation\n"
        "overlayed with a binary time series. Suitable for analyzing orientational diffusion\n"
        "of water molecules confined in certain region, which is described by the binary\n"
        "time series.",
    },

    { NULL, NULL, 0, NULL }
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

MOD_INIT(_utils_ext)
{
    PyObject *m;

    MOD_DEF(m, "_utils_ext", "Different routines for analyzing\n"
                             "diffusion properties of solvent molecules.\n", _utils_ext_methods)
    
    if (m == NULL)
        return MOD_ERROR_VAL;

    import_array();

    return MOD_SUCCESS_VAL(m);
}