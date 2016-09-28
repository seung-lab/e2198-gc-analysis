
#defining NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

# @cython.boundscheck(False) # shouldn't need this, we round down partial windows
# @cython.wraparound(False)
@cython.nonecheck(False)

cpdef np.ndarray[DTYPE_t, ndim=3] sliding_window( np.ndarray[DTYPE_t,ndim=3] data3d, 
						  function, 
						  tuple filter_size, 
						  tuple stride):


	cdef tuple data_shape = (data3d.shape[0], data3d.shape[1], data3d.shape[2])

	cdef np.ndarray[np.uint32_t, ndim=1] ending = np.zeros((3,),dtype=np.uint32)
	cdef np.ndarray[np.uint32_t, ndim=1] beginning = np.zeros((3,),dtype=np.uint32)
	#ROUNDS DOWN TO CLIP PARTIAL WINDOWS
	cdef np.ndarray[np.uint32_t, ndim=1] final_shape = output_shape( data_shape, filter_size, stride )

	cdef np.ndarray[DTYPE_t, ndim=3] output = np.empty( final_shape, dtype=DTYPE )
	cdef np.ndarray[DTYPE_t, ndim=3] data_window

	cdef np.uint32_t z,y,x

	for z in xrange(final_shape[0]):
		for y in xrange(final_shape[1]):
			for x in xrange(final_shape[2]):

				# at least one of the factors should be the same
				# type as beginning, otherwise you'll see an error
				# after compilation about unexpected types
				beginning[0] = stride[0] * z
				beginning[1] = stride[1] * y
				beginning[2] = stride[2] * x

				ending[0] = beginning[0] + filter_size[0]
				ending[1] = beginning[1] + filter_size[1]
				ending[2] = beginning[2] + filter_size[2]

				data_window = data3d[
					beginning[0]:ending[0],
					beginning[1]:ending[1],
					beginning[2]:ending[2]]

				output[z,y,x] = function( data_window )

	return output


cdef np.ndarray[np.uint32_t,ndim=1] output_shape( tuple data_shape, tuple filter_size, tuple stride ):

	cdef tuple shape_w_o_stride = (
		data_shape[0] - filter_size[0] + 1,
		data_shape[1] - filter_size[1] + 1,
		data_shape[2] - filter_size[2] + 1)

	return np.floor( shape_w_o_stride / np.array(stride).astype('float64') ).astype('uint32')


cpdef np.ndarray[DTYPE_t, ndim=3] expand_volume( np.ndarray[DTYPE_t,ndim=3] data3d, tuple stride ):

	#this has to be the same type as beginning and ending
	cdef np.ndarray[np.uint32_t, ndim=1] data_shape = np.array((data3d.shape[0], data3d.shape[1], data3d.shape[2]), dtype=np.uint32)

	cdef tuple final_shape = (data_shape[0] * stride[0],
							  data_shape[1] * stride[1],
							  data_shape[2] * stride[2])

	cdef np.ndarray[DTYPE_t, ndim=3] output = np.empty( final_shape, dtype=DTYPE )

	cdef np.ndarray[np.uint32_t, ndim=1] beginning = np.zeros((3,),dtype=np.uint32)
	cdef np.ndarray[np.uint32_t, ndim=1] ending    = np.zeros((3,),dtype=np.uint32)

	cdef np.uint32_t z,y,x	

	for z in range(data_shape[0]):
		for y in range(data_shape[1]):
			for x in range(data_shape[2]):

				beginning[0] = stride[0] * z
				beginning[1] = stride[1] * y
				beginning[2] = stride[2] * x

				ending[0] = beginning[0] + stride[0]
				ending[1] = beginning[1] + stride[1]
				ending[2] = beginning[2] + stride[2]

				output[
					beginning[0]:ending[0],
					beginning[1]:ending[1],
					beginning[2]:ending[2]] = data3d[z,y,x]

	return output
