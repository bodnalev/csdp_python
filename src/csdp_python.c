#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <Python.h>
#include "declarations.h"


PyObject* c_solve_sdp(PyObject* self, PyObject* args) {
	PyObject *block_sizes_obj;
	PyObject *a_obj;
	PyObject *mat_inds_obj;
	PyObject *mat_vals_obj;
	if (!PyArg_ParseTuple(args, "OOOO", &block_sizes_obj, &a_obj, &mat_inds_obj, &mat_vals_obj)){
		return NULL;
	}
	int block_num = (int) PyObject_Length(block_sizes_obj);
	int k = (int) PyObject_Length(a_obj);
	int rows = (int) PyObject_Length(mat_vals_obj);
	
	int* block_sizes = (int*)malloc(block_num * sizeof(int));
	for (int i = 0; i < block_num; i++) {
		block_sizes[i] = (int)PyLong_AsLong(PyList_GetItem(block_sizes_obj, i));
	}
	
	double* a = (double*)malloc((k + 1) * sizeof(double)); /* idk why the list is shifted by 1 */
	for (int i = 1; i <= k; i++) {
		a[i] = (double)PyFloat_AsDouble(PyList_GetItem(a_obj, i-1));
	}
	
	int* mat_inds = (int*)malloc(4 * rows * sizeof(int));
	for (int i = 0; i < 4 * rows; i++) {
		mat_inds[i] = (int)PyLong_AsLong(PyList_GetItem(mat_inds_obj, i));
	}
	
	double* mat_vals = (double*)malloc(rows * sizeof(double));
	for (int i = 0; i < rows; i++) {
		mat_vals[i] = (double)PyFloat_AsDouble(PyList_GetItem(mat_vals_obj, i));
	}
	
	
	int ret;
	int n;
	struct blockmatrix C;
	struct constraintmatrix *constraints;
	struct blockmatrix X,Z;
	double *y;
	double pobj,dobj;
	double rval;
	
	ret=from_sparse_data(k,block_num,block_sizes,rows,mat_inds,mat_vals, &n, &C, &constraints, 1);
	if (ret != 0)
	{
		printf("Something is wrong with the problem\n");
		exit(201);
	};
	initsoln(n,k,C,a,constraints,&X,&y,&Z);
	ret=easy_sdp(n,k,C,a,constraints,0.0,&X,&y,&Z,&pobj,&dobj);
	if (ret != 0)
	{
		printf("Something is wrong with trying to solve the problem\n");
		exit(201);
	};
	rval = pobj;
	free_prob(n,k,C,a,constraints,X,y,Z);
	free(block_sizes);
	free(mat_inds);
	free(mat_vals);
	fflush(stdout);
	
	if (PyErr_Occurred()) {
		PyObject *exc_type, *exc_value, *exc_tb;
		PyErr_Fetch(&exc_type, &exc_value, &exc_tb);
		PyErr_Print();
		Py_XDECREF(exc_type);
		Py_XDECREF(exc_value);
		Py_XDECREF(exc_tb);
		return NULL;
    }
	
	return Py_BuildValue("d", rval);
	}


PyMethodDef methods[] = {
	{"solve_sdp", c_solve_sdp, METH_VARARGS, "Run the ccsdp solver"},
	{NULL, NULL, 0, NULL}
};


struct PyModuleDef module = {
	PyModuleDef_HEAD_INIT,
	"csdpy",
	NULL,
	-1,
	methods
};

PyMODINIT_FUNC PyInit_csdpy(void) {
	return PyModule_Create(&module);
}