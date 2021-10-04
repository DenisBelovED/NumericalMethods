#ifndef PY_SSIZE_T_CLEAN
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdio.h>
#endif

/*
PyObject* tanh_impl(PyObject*, PyObject* o)
{
    double x = PyFloat_AsDouble(o);
    double tanh_x = sinh(x) / cosh(x);
    return PyFloat_FromDouble(tanh_x);
}
*/

PyObject* half_division(PyObject*, PyObject* args)
{
    double r, c, left_edge, right_edge, epsilon;
    double (*f)(double) = [](double x) { return pow(x, 5) + pow(x, 2) - 5; };

    if (PyArg_ParseTuple(args, "ddd", &left_edge, &right_edge, &epsilon))
        printf("Get args: [%.9f, %.9f] eps=%.9f\n", left_edge, right_edge, epsilon);
    else
    {
        printf("Wrong args, call example: half_division(1.2, 3.4, 0.001)\n");
        return nullptr;
    }
    if (left_edge >= right_edge)
    {
        printf("Wrong input range: left edge shoud < right edge\n");
        return nullptr;
    }
    if ((f(left_edge) * f(right_edge) >= 0) || (epsilon <= 0))
    {
        printf("Wrong initial conditions, must be f(a) * f(b) < 0 on [a, b], epsilon > 0\n");
        return nullptr;
    }

    PyObject* list = PyList_New(0);
    PyList_Append(list, PyTuple_Pack(2,  PyFloat_FromDouble(left_edge), PyFloat_FromDouble(f(left_edge))));
    PyList_Append(list, PyTuple_Pack(2, PyFloat_FromDouble(right_edge), PyFloat_FromDouble(f(right_edge))));
    while (right_edge - left_edge > epsilon)
    {
        c = (right_edge + left_edge) / 2;
        r = f(c);
        PyList_Append(list, PyTuple_Pack(2, PyFloat_FromDouble(c), PyFloat_FromDouble(r)));
        if (f(right_edge) * r < 0)
            left_edge = c;
        else
            right_edge = c;
    }
    return list;
}

static PyMethodDef compile_methods[] = {
    // The first property is the name exposed to Python, fast_tanh
    // The second is the C++ function with the implementation
    // METH_O means it takes a single PyObject argument
    // { "fast_tanh", (PyCFunction)tanh_impl, METH_O, nullptr },
    { "half_division", (PyCFunction)half_division, METH_VARARGS, nullptr},
    // Terminate the array with an object containing nulls.
    { nullptr, nullptr, 0, nullptr }
};

static PyModuleDef compile_module = {
    PyModuleDef_HEAD_INIT,
    "NumericalMethods",                              // Module name to use with Python import statements
    "Provides implementation numerical methods",     // Module description
    0,
    compile_methods                                  // Structure that defines the methods of the module
};

PyMODINIT_FUNC PyInit_NumericalMethods() {
    return PyModule_Create(&compile_module);
}