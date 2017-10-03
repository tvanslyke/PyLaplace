#include <python3.5/Python.h>
#include "structmember.h"
#include "laplace.h"
#include <string>

static auto py_object_as_callable(PyObject* arg)
{
	if(not PyCallable_Check(arg))
		throw std::runtime_error("Attempt to create C++ lambda function from non-callable Python object.");
	Py_INCREF(arg);
	return [=](double d)->double
	{
		PyObject* num = PyFloat_FromDouble(d);
		if(!arg)
			throw std::runtime_error("Failed to create Python float from C double.");
		PyObject* result = PyObject_CallFunctionObjArgs(arg, num, NULL);
		if(!result)
			throw std::runtime_error("Bad result recieved from callable Python object.");
		Py_DECREF(result);
		d = PyFloat_AsDouble(num);
		if(d == -1.0 and PyErr_Occurred())
			throw std::runtime_error("Error occurred while converting Python float to C double.");
		Py_DECREF(result);
		return d;
	};
}

extern "C"{

struct PyLaplaceObject
{
	PyObject_HEAD
	LaplaceTransform<double> lt_;
};

static void PyLaplaceObject_dealloc(PyObject* self)
{
	PyLaplaceObject* slf = (PyLaplaceObject*)self;
	slf->lt_.~LaplaceTransform<double>();
}
static PyObject* PyLaplaceObject_new(PyTypeObject* type, PyObject* args, PyObject* kwargs)
{
	PyObject* self = type->tp_alloc(type, 0);
	if(!self)
	{
		PyErr_NoMemory();
		return nullptr;
	}
	return self;
}

static int PyLaplaceObject_init(PyObject* self, PyObject* args, PyObject* kwargs)
{
	std::string s_kw{"s"};
	std::string order_kw{"order"};
	char* keywords[3]{s_kw.data(), order_kw.data(), nullptr};
	PyObject* function;
	unsigned PY_LONG_LONG order{10};
	if(not PyArg_ParseTupleAndKeywords(args, kwargs, "O|K", keywords, &function, &order))
		return -1;
	if(order < 1)
	{
		Py_DECREF(function);
		Py_DECREF(self);
		PyErr_SetString(PyExc_TypeError, "'order' must be greater than 0 in laplace_transform().");
		return -1;
	}
	try
	{
		auto func = py_object_as_callable(function);
		new(&((PyLaplaceObject*)(self))->lt_) LaplaceTransform<double>(laplace_transform(func, order));
	}
	catch(const std::exception & e)
	{
		Py_DECREF(self);
		Py_DECREF(function);
		PyErr_SetString(PyExc_RuntimeError, e.what());
		return -1;
	}
	return 0;
}

static PyObject* PyLaplaceObject_call(PyObject* self, PyObject* arg)
{
	PyLaplaceObject* slf = (PyLaplaceObject*)self;
	std::complex<double> num{0.0, 0.0};
	try{
		if(PyComplex_Check(arg))
		{
			num.real(PyComplex_RealAsDouble(arg));
			num.imag(PyComplex_ImagAsDouble(arg));
			num = slf->lt_(num);
			return PyComplex_FromDoubles(num.real(), num.imag());
		}
		PyObject* real = PyNumber_Float(arg);
		if(!real)
			return nullptr;
		num.real(PyFloat_AsDouble(real));
		Py_DECREF(real);
		if(num.real() == -1.0 and PyErr_Occurred())
			return nullptr;
		num = slf->lt_(num);
		return PyComplex_FromDoubles(num.real(), num.imag());	
	}
	catch(const std::exception& e)
	{
		PyErr_SetString(PyExc_RuntimeError, e.what());
		return nullptr;
	}
}
static PyObject* PyLaplaceObject_call_func(PyObject* self, PyObject* args, PyObject* kwargs)
{
	char* kws[1] = {nullptr};
	PyObject* x;
	if(!PyArg_ParseTuple(args, "O", &x))
		return nullptr;

	return PyLaplaceObject_call(self, x);
}

static PyMethodDef PyLaplaceObject_methods[] = 
{
	{"__call__", (PyCFunction)PyLaplaceObject_call, METH_O,
       		"Evaluate the LaplaceTransform object at a point 's' in the complex plane." 
	},
	{nullptr}  /* Sentinel */
};



static PyTypeObject PyLaplaceObject_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "laplace.LaplaceTransform",/* tp_name */
    sizeof(PyLaplaceObject),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)PyLaplaceObject_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    PyLaplaceObject_call_func, /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "Laplace Transform of a function.",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    (PyMethodDef*)(&PyLaplaceObject_methods),   /* tp_methods */
    0,		               /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)PyLaplaceObject_init,      /* tp_init */
    0,                         /* tp_alloc */
    PyLaplaceObject_new,                 /* tp_new */
};


static PyModuleDef laplacemodule = {
	PyModuleDef_HEAD_INIT,
	"laplace",
	"Provides a LaplaceTransform object to numerically compute the laplace transform of arbitrary Python functions.",
	-1,
	NULL, NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC
PyInit_laplace(void)
{
	PyObject* m;

	if (PyType_Ready(&PyLaplaceObject_Type) < 0)
		return nullptr;

	m = PyModule_Create(&laplacemodule);
	if (m == nullptr)
		return nullptr;

	Py_INCREF(&PyLaplaceObject_Type);
	PyModule_AddObject(m, "LaplaceTransform", (PyObject *)&PyLaplaceObject_Type);
	return m;
}

} /* extern "C" */
