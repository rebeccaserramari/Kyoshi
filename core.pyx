# distutils: language=c++
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.unordered_set cimport unordered_set
cimport cpp

from collections import namedtuple
from cython.operator cimport dereference as deref



def hmm(haplo1, haplo2, E_whole, panel, mutation, recomb_list, posterior_h1, posterior_h2):
	filename_byte_string_post_h1 = posterior_h1.encode("UTF-8")
	filename_byte_string_post_h2 = posterior_h2.encode("UTF-8")
	cdef char* fname_posterior_h1 = filename_byte_string_post_h1
	cdef char* fname_posterior_h2 = filename_byte_string_post_h2
	cdef float mutation_parameter = int(mutation)
	c_recomblist = new vector[float]()
	for rec in recomb_list:
		c_recomblist.push_back(rec)
	E = new vector[int]()
	for boundary in E_whole:
		E.push_back(boundary)
	cdef vector[int]* haplo1_vec = new vector[int]()
	cdef vector[int]* haplo2_vec = new vector[int]()
	for hap in haplo1:
		haplo1_vec.push_back(hap)
	for hap in haplo2:
		haplo2_vec.push_back(hap)
	cdef vector[string]* ref = new vector[string]()
	for row in panel:
		ref.push_back(row.encode('UTF-8'))
	cdef vector[int]* path1 = new vector[int]()
	cdef vector[int]* path2 = new vector[int]()
	cpp.compute_hapcolor(haplo1_vec, haplo2_vec, E, ref, mutation_parameter, path1, path2, c_recomblist, fname_posterior_h1, fname_posterior_h2)
	del ref
	del haplo1_vec
	del haplo2_vec
	del E
	del c_recomblist
	return deref(path1), deref(path2)

