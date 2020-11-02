# kate: syntax python
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.unordered_set cimport unordered_set


cdef extern from "hapcolor.h":
	void compute_hapcolor(vector[int]*, vector[int]*, vector[int]*, vector[string]*, float, vector[int]*, vector[int]*, vector[float]*, char*, char*) except + 



