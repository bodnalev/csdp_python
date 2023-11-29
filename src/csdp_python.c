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
		PyObject *sublist = PyList_New(block->blocksize);
		if (!sublist) {
			Py_DECREF(result);
			PyErr_SetString(PyExc_MemoryError, "Failed to create Python list");
			return NULL;
		}
		for (int j = 1; j <= block->blocksize; ++j) {
			value = 0.0;
			if (block->blockcategory==MATRIX){
				value = block->data.mat[ijtok(i,j,block->blocksize)];
			}
			else{
				if (i==j){
					value = block->data.vec[i];
				}
			}
			PyObject *py_value = PyFloat_FromDouble(value);
			if (!py_value) {
				Py_DECREF(result);
				Py_DECREF(sublist);
				PyErr_SetString(PyExc_MemoryError, "Failed to create Python float");
				return NULL;
			}
			PyList_SET_ITEM(sublist, j - 1, py_value);
		}
		PyList_SET_ITEM(result, i - 1, sublist);
	}
	return result;
}

PyObject* convert_blockmatrix_to_python_list(struct blockmatrix *matrix) {
	PyObject *result = PyList_New(matrix->nblocks);
	if (!result) {
		PyErr_SetString(PyExc_MemoryError, "Failed to create Python list");
		return NULL;
	}
	for (int i = 0; i < matrix->nblocks; ++i) {
		PyObject *block_list = convert_blockrec_to_python_list(&matrix->blocks[i]);
		if (!block_list) {
			Py_DECREF(result);
			return NULL;
		}
		PyList_SET_ITEM(result, i, block_list);
	}
	return result;
}

PyObject* c_solve_sdp(PyObject* self, PyObject* args) {
	PyObject *block_sizes_obj;
	PyObject *a_obj;
	PyObject *mat_inds_obj;
	PyObject *mat_vals_obj;
	
	if (!PyArg_ParseTuple(args, "OOOO", &block_sizes_obj, &a_obj, &mat_inds_obj, &mat_vals_obj)){
		PyErr_SetString(PyExc_ValueError, "CSDPY couldn't parse the input");
		return NULL;
	}
	
	int block_num = (int) PyObject_Length(block_sizes_obj);
	int k = (int) PyObject_Length(a_obj);
	int rows = (int) PyObject_Length(mat_vals_obj);
	
	
	
	int* block_sizes = (int*)malloc(block_num * sizeof(int));
	if (!block_sizes){
		PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
		return NULL;
	}
	for (int i = 0; i < block_num; i++) {
		block_sizes[i] = (int)PyLong_AsLong(PyList_GetItem(block_sizes_obj, i));
	}
	
	
	
	double* a = (double*)malloc((k + 1) * sizeof(double)); /* idk why the list is shifted by 1 */
	if (!a){
		PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
		return NULL;
	}
	for (int i = 1; i <= k; i++) {
		a[i] = (double)PyFloat_AsDouble(PyList_GetItem(a_obj, i-1));
	}
	
	
	int* mat_inds = (int*)malloc(4 * rows * sizeof(int));
	if (!mat_inds){
		PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
		return NULL;
	}
	for (int i = 0; i < 4 * rows; i++) {
		mat_inds[i] = (int)PyLong_AsLong(PyList_GetItem(mat_inds_obj, i));
	}
	
	
	
	double* mat_vals = (double*)malloc(rows * sizeof(double));
	if (!mat_vals){
		PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
		return NULL;
	}
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
	
	ret=from_sparse_data(k,block_num,block_sizes,rows,mat_inds,mat_vals, &n, &C, &constraints, 1);
	
	free(block_sizes);
	free(mat_inds);
	free(mat_vals);
	
	if (ret != 0)
	{
		PyErr_SetString(PyExc_ValueError, "Couldn't create an SDP from the data");
		free_prob(n,k,C,a,constraints,X,y,Z);
		return NULL;
	}
	
	initsoln(n,k,C,a,constraints,&X,&y,&Z);
	
	ret=easy_sdp(n,k,C,a,constraints,0.0,&X,&y,&Z,&pobj,&dobj);
	
	PyObject* result_obj = PyDict_New();
	if (!result_obj){
		PyErr_SetString(PyExc_MemoryError, "Failed to create Python dictionary");
		free_prob(n,k,C,a,constraints,X,y,Z);
		return NULL;
	}
	
	
	PyObject* X_obj = convert_blockmatrix_to_python_list(&X);
	if (!X_obj){
		PyErr_SetString(PyExc_ValueError, "Couldn't convert the result to python objects");
		free_prob(n,k,C,a,constraints,X,y,Z);
		Py_DECREF(result_obj);
		return NULL;
	}
	
	
	PyObject* Z_obj = convert_blockmatrix_to_python_list(&Z);
	if (!Z_obj){
		PyErr_SetString(PyExc_ValueError, "Couldn't convert the result to python objects");
		free_prob(n,k,C,a,constraints,X,y,Z);
		Py_DECREF(result_obj);
		Py_DECREF(X_obj);
		return NULL;
	}
	
	
	PyObject* y_obj = PyList_New(k);
	if (!y_obj){
		PyErr_SetString(PyExc_ValueError, "Couldn't convert the result to python objects");
		free_prob(n,k,C,a,constraints,X,y,Z);
		Py_DECREF(result_obj);
		Py_DECREF(X_obj);
		Py_DECREF(Z_obj);
		return NULL;
	}
	for (int i=1; i<=k; i++){
		PyObject *py_value = PyFloat_FromDouble(y[i]);
		if (!py_value){
			PyErr_SetString(PyExc_ValueError, "Couldn't convert the result to python objects");
			free_prob(n,k,C,a,constraints,X,y,Z);
			Py_DECREF(result_obj);
			Py_DECREF(X_obj);
			Py_DECREF(Z_obj);
			Py_DECREF(y_obj);
			return NULL;
		}
		PyList_SET_ITEM(y_obj, i - 1, py_value);
	}
	
	PyObject* pobj_obj = PyFloat_FromDouble(pobj);
	if (!pobj_obj){
		PyErr_SetString(PyExc_ValueError, "Couldn't convert the result to python objects");
		free_prob(n,k,C,a,constraints,X,y,Z);
		Py_DECREF(result_obj);
		Py_DECREF(X_obj);
		Py_DECREF(Z_obj);
		Py_DECREF(y_obj);
		return NULL;
	}
	PyObject* dobj_obj = PyFloat_FromDouble(dobj);
	if (!dobj_obj){
		PyErr_SetString(PyExc_ValueError, "Couldn't convert the result to python objects");
		free_prob(n,k,C,a,constraints,X,y,Z);
		Py_DECREF(result_obj);
		Py_DECREF(X_obj);
		Py_DECREF(Z_obj);
		Py_DECREF(y_obj);
		Py_DECREF(pobj_obj);
		return NULL;
	}
	PyObject* ret_obj = PyLong_FromLong(ret);
	if (!ret_obj){
		PyErr_SetString(PyExc_ValueError, "Couldn't convert the result to python objects");
		free_prob(n,k,C,a,constraints,X,y,Z);
		Py_DECREF(result_obj);
		Py_DECREF(X_obj);
		Py_DECREF(Z_obj);
		Py_DECREF(y_obj);
		Py_DECREF(pobj_obj);
		Py_DECREF(dobj_obj);
		return NULL;
	}
	
	PyDict_SetItemString(result_obj, "code", ret_obj);
	PyDict_SetItemString(result_obj, "primal", pobj_obj);
	PyDict_SetItemString(result_obj, "dual", dobj_obj);
	PyDict_SetItemString(result_obj, "X", X_obj);
	PyDict_SetItemString(result_obj, "y", y_obj);
	PyDict_SetItemString(result_obj, "Z", Z_obj);
	
	
	free_prob(n,k,C,a,constraints,X,y,Z);
	/*Py_DECREF(X_obj);
	Py_DECREF(Z_obj);
	Py_DECREF(y_obj);
	Py_DECREF(pobj_obj);
	Py_DECREF(dobj_obj);
	Py_DECREF(ret_obj);*/
	fflush(stdout);
	
	return result_obj;
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