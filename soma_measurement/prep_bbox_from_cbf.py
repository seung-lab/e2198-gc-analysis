#!/usr/bin/env python

import h5py
import u
import numpy as np
from omni2h5 import ChunkedBinary
from prep_chunked import prep_chunked_file
from sys import argv

output_filename = argv[1]

# fname = xxx
# shape = (6528,10496,2432)
#shape = (128,256,2432)
# chunk_shape = (128,128,128)
# dtype = 'uint8'
#fname = xxx
shape = (2048,2048,896)
chunk_shape = (128,128,128)
dtype = 'uint32'

#size of chunk to process at once
input_chunk_shape = np.array([ 128, 128, 128])

beginning  = np.array([ 0,   0,   0  ])
ending     = np.array([ 896, 2048, 896])

processing_stride = np.array([ 1,  1,  1])
processing_filter = np.array([ 1, 1, 1])


def process_chunk(data):
	return data
	# return u.minfilter( data, (5,9,9), (5,5,5) )


print "Reading data"
data = ChunkedBinary( fname, shape, chunk_shape, dtype )

output = prep_chunked_file( data, process_chunk,
	beginning, ending,
	input_chunk_shape,
	processing_filter,
	processing_stride )

u.write_chunked_h5( output, output_filename )
