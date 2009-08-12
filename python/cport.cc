// MOODS.xs - algorithm for finding PWM matches from sequences
// Copyright (C) 2009 Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <Python.h>
#include "pssm_algorithms.hpp"

charArray convertSequence(const char *sequence) {
	charArray c_seq;
	int lenght = strlen(sequence);
	for(int i = 0; i < lenght; i++) {
		char toput = 5;
		switch (sequence[i])
		{
			case 'a':
			case 'A': toput = 0; break;
			case 'c':
			case 'C': toput = 1; break;
			case 'g':
			case 'G': toput = 2; break;
			case 't':
			case 'T': toput = 3; break;
			case 'n':
			case 'N': toput = (int) (4 * rand() / (1.0 + RAND_MAX)); break;
			default:
				break;
		}
		if(toput != 5) {
			c_seq.push_back(toput);
		}
	}
	return c_seq;
}

scoreArray atoDoubleArray(PyObject *o) {
	scoreArray t;
	if(!PyList_Check(o)) {
		return t;
	}
	Py_ssize_t length = PyList_Size(o);
	PyObject *tmp;
	for(int i=0; i< length; i++) {
		tmp = PyList_GET_ITEM(o, i);
		t.push_back(PyFloat_AsDouble(tmp));
	}
	return t;
}
intArray atoIntArray(PyObject *o) {
	intArray t;
	if(!PyList_Check(o)) {
		return t;
	}
	Py_ssize_t length = PyList_Size(o);
	PyObject *tmp;
	for(int i=0; i< length; i++) {
		tmp = PyList_GET_ITEM(o, i);
		t.push_back((int)PyFloat_AsDouble(tmp));
	}
	return t;
}
intMatrix atoIntMatrix(PyObject *o) {
	intMatrix t;
	if(!PyList_Check(o)) {
		return t;
	}
	Py_ssize_t length = PyList_Size(o);
	PyObject *tmp;
	for(int i=0; i< length; i++) {
		tmp = PyList_GET_ITEM(o, i);
		t.push_back(atoIntArray(tmp));
	}
	return t;
}
scoreMatrix atoDoubleMatrix(PyObject *o) {
	scoreMatrix t;
	if(!PyList_Check(o)) {
		return t;
	}
	Py_ssize_t length = PyList_Size(o);
	PyObject *tmp;
	for(int i=0; i< length; i++) {
		tmp = PyList_GET_ITEM(o, i);
		t.push_back(atoDoubleArray(tmp));
	}
	return t;
}

static PyObject *
_search(PyObject *self, PyObject *args)
{
    const char *sequence;
    charArray c_seq;
    PyObject *py_matrices;
    PyObject *py_thresholds;
    const char *algorithm;
    int q;
    PyObject *py_absolute_threshold;
    PyObject *py_bg;
    PyObject *py_combine;
    bool absolute_threshold;
    double ps = 0.1;
    bool combine = false;

    std::vector<scoreMatrix> matrices;
    scoreArray thresholds;
    scoreArray bg;

    if (!PyArg_ParseTuple(args, "sOOOsiOO", &sequence, &py_matrices, &py_thresholds, &py_bg, &algorithm, &q, &py_absolute_threshold, &py_combine))
        return NULL;

    c_seq = convertSequence(sequence);
    bg = atoDoubleArray(py_bg);
    absolute_threshold = (bool) PyObject_IsTrue(py_absolute_threshold);
    combine = (bool) PyObject_IsTrue(py_combine);
    thresholds = atoDoubleArray(py_thresholds);

    if(!PyList_Check(py_matrices)) {
		return NULL;
	}
	int num_matrices = (int) PyList_Size(py_matrices);
	if(num_matrices != thresholds.size()) {
		//TODO Raise exception
		return NULL;
	}
	for(int i=0; i< num_matrices; i++) {
		if(absolute_threshold) {
			matrices.push_back(atoDoubleMatrix(PyList_GET_ITEM(py_matrices, i)));
		}
		else {
			matrices.push_back(counts2LogOdds(atoIntMatrix(PyList_GET_ITEM(py_matrices, i)), bg, ps));
			thresholds[i] = tresholdFromP(matrices.back(), bg, thresholds[i]);
		}
		if(matrices[i].size() != 4) {
			//TODO Raise exception
			return NULL;
		}
	}
	if(matrices.size() == 1)
		combine = false;
	std::vector<matchArray> matches;
	if(strcmp(algorithm, "lf")== 0 && combine) {
		matches = multipleMatrixLookaheadFiltrationDNA(q, c_seq, matrices, bg, thresholds);
	}
	else {
		for(int i = 0; i < num_matrices; i++) {
			matchArray result;
			if(maxScore(matrices[i]) >= thresholds[i]) {
				if(strcmp(algorithm, "naive") == 0) {
					result = naiveAlgorithm(c_seq, matrices[i], thresholds[i]);
				}
				else if(strcmp(algorithm, "pla") == 0) {
					result = permutatedLookAhead(c_seq, matrices[i], bg, thresholds[i]);
				}
				else if(strcmp(algorithm, "supera") == 0 ||  matrices[0].size() > (3*q)) {
					result = naiveSuperalphabetAlgorithmDNA(q, c_seq, matrices[i], thresholds[i]);
				}
				else if(strcmp(algorithm, "lf") == 0) {
					result = lookaheadFiltrationDNA(q, c_seq, matrices[i], bg, thresholds[i]);
				}
				else {
					//TODO Raise unsupported algorithm exception
				}
			}
			matches.push_back(result);
		}
	}
	PyObject *results = PyList_New(matches.size());
	for(int i = 0; i < matches.size(); i++) {
		PyObject *new_match_list = PyList_New(matches[i].size());
		for(int j=0; j < matches[i].size(); j++) {
			PyList_SET_ITEM(new_match_list, j, Py_BuildValue("Ld", matches[i][j].position, matches[i][j].score));
		}
		PyList_SET_ITEM(results, i, new_match_list);
	}
    return results;
}

static PyMethodDef SpamMethods[] = {
    {"search",  _search, METH_VARARGS,
     "Search for a pwm match."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

extern "C" {
PyMODINIT_FUNC
init_cmodule(void)
{
	(void) Py_InitModule("_cmodule", SpamMethods);
}

}
