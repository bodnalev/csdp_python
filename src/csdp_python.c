#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <Python.h>
#include "declarations.h"

PyObject* convert_blockrec_to_python_list(struct blockrec *block) {
	double value = 0.0;
	PyObject *result = PyList_New(block->blocksize);
	if (!result) {
		PyErr_SetString(PyExc_MemoryError, "Failed to create Python list");
		return NULL;
	}
	for (int i = 1; i <= block->blocksize; ++i) {
		if (block->blockcategory==MATRIX){
			PyObject *sublist = PyList_New(block->blocksize);
			if (!sublist) {
				Py_DECREF(result);
				PyErr_SetString(PyExc_MemoryError, "Failed to create Python list");
				return NULL;
			}
			for (int j = 1; j <= block->blocksize; ++j) {
				value = block->data.mat[ijtok(i,j,block->blocksize)];
				PyObject *py_value = PyFloat_FromDouble(value);
				PyList_SET_ITEM(sublist, j - 1, py_value);
			}
			PyList_SET_ITEM(result, i - 1, sublist);
		}
		else{
			value = block->data.vec[i];
			PyObject *py_value = PyFloat_FromDouble(value);
			PyList_SET_ITEM(result, i - 1, py_value);
		}
	}
	return result;
}

PyObject* convert_blockmatrix_to_python_list(struct blockmatrix *matrix) {
	PyObject *result = PyList_New(matrix->nblocks);
	if (!result) {
		PyErr_SetString(PyExc_MemoryError, "Failed to create Python list");
		return NULL;
	}
	for (int i = 1; i <= matrix->nblocks; ++i) {
		PyObject *block_list = convert_blockrec_to_python_list(&matrix->blocks[i]);
		if (!block_list) {
			Py_DECREF(result);
			return NULL;
		}
		PyList_SET_ITEM(result, i - 1, block_list);
	}
	return result;
}

PyObject* c_solve_sdp(PyObject* self, PyObject* args) {
	PyObject *block_sizes_obj;
	PyObject *a_obj;
	PyObject *mat_inds_obj;
	PyObject *mat_vals_obj;
	
	PyObject *result_obj;
	PyObject *X_obj;
	PyObject *y_obj;
	PyObject *Z_obj;
	PyObject *pobj_obj;
	PyObject *dobj_obj;
	PyObject *ret_obj;
	
	int *block_sizes;
	double *a;
	int *mat_inds;
	double *mat_vals;
	
	int ret;
	int n;
	struct blockmatrix C;
	struct constraintmatrix *constraints;
	struct blockmatrix X,Z;
	double *y;
	double pobj,dobj;
	
	if (!PyArg_ParseTuple(args, "OOOO", &block_sizes_obj, &a_obj, &mat_inds_obj, &mat_vals_obj)){
		PyErr_SetString(PyExc_ValueError, "Couldn't parse the input");
		return NULL;
	}
	
	int block_num = (int) PyObject_Length(block_sizes_obj);
	int k = (int) PyObject_Length(a_obj);
	int rows = (int) PyObject_Length(mat_vals_obj);
	
	
	block_sizes = (int*)malloc(block_num * sizeof(int));
	if (!block_sizes){
		PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
		return NULL;
	}
	for (int i = 0; i < block_num; i++) {
		block_sizes[i] = (int)PyLong_AsLong(PyList_GetItem(block_sizes_obj, i));
	}
	
	
	a = (double*)malloc((k + 1) * sizeof(double));
	if (!a){
		PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
		free(block_sizes);
		return NULL;
	}
	for (int i = 1; i <= k; i++) {
		a[i] = (double)PyFloat_AsDouble(PyList_GetItem(a_obj, i-1));
	}
	
	
	mat_inds = (int*)malloc(4 * rows * sizeof(int));
	if (!mat_inds){
		PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
		free(block_sizes);
		free(a);
		return NULL;
	}
	for (int i = 0; i < 4 * rows; i++) {
		mat_inds[i] = (int)PyLong_AsLong(PyList_GetItem(mat_inds_obj, i));
	}
	
	
	mat_vals = (double*)malloc(rows * sizeof(double));
	if (!mat_vals){
		PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
		free(block_sizes);
		free(a);
		free(mat_inds);
		return NULL;
	}
	for (int i = 0; i < rows; i++) {
		mat_vals[i] = (double)PyFloat_AsDouble(PyList_GetItem(mat_vals_obj, i));
	}
	
	
	
	ret = from_sparse_data(k,block_num,block_sizes,rows,mat_inds,mat_vals, &n, &C, &constraints, 1);
	
	if (ret != 0)
	{
		PyErr_SetString(PyExc_ValueError, "Couldn't create an SDP from the data");
		free(block_sizes);
		free(a);
		free(mat_inds);
		free(mat_vals);
		return NULL;
	}
	
	initsoln(n,k,C,a,constraints,&X,&y,&Z);
	ret = easy_sdp(n,k,C,a,constraints,0.0,&X,&y,&Z,&pobj,&dobj);
	
	result_obj = PyDict_New();
	
	X_obj = convert_blockmatrix_to_python_list(&X);
	Z_obj = convert_blockmatrix_to_python_list(&Z);
	y_obj = PyList_New(k);
	for (int i=1; i<=k; i++){
		PyObject *py_value = PyFloat_FromDouble(y[i]);
		PyList_SET_ITEM(y_obj, i - 1, py_value);
	}
	
	pobj_obj = PyFloat_FromDouble(pobj);
	dobj_obj = PyFloat_FromDouble(dobj);
	ret_obj = PyLong_FromLong(ret);
	
	if (!X_obj || !Z_obj || !y_obj || !pobj_obj || !dobj_obj || !ret_obj || !result_obj){
		free(block_sizes);
		free(mat_inds);
		free(mat_vals);
		free_prob(n,k,C,a,constraints,X,y,Z);
		Py_XDECREF(X_obj);
		Py_XDECREF(Z_obj);
		Py_XDECREF(y_obj);
		Py_XDECREF(pobj_obj);
		Py_XDECREF(dobj_obj);
		Py_XDECREF(ret_obj);
		Py_XDECREF(result_obj);
	}
	
	PyDict_SetItemString(result_obj, "code", ret_obj);
	PyDict_SetItemString(result_obj, "primal", pobj_obj);
	PyDict_SetItemString(result_obj, "dual", dobj_obj);
	PyDict_SetItemString(result_obj, "X", X_obj);
	PyDict_SetItemString(result_obj, "y", y_obj);
	PyDict_SetItemString(result_obj, "Z", Z_obj);
	
	free(block_sizes);
	free(mat_inds);
	free(mat_vals);
	free_prob(n,k,C,a,constraints,X,y,Z);
	
	fflush(stdout);
	return result_obj;
}

PyMethodDef methods[] = {
	{"solve_sdp", c_solve_sdp, METH_VARARGS, "Run the csdp solver"},
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