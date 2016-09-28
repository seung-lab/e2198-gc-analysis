#!/usr/bin/env python

import os, h5py, timeit
import numpy as np
import u
from ChunkedBinary import ChunkedBinary

#NOTE - the last time I used this code, the ChunkedBinary class was written within
# the same file, take care when running again

fname = xxx
shape = (6528,10496,2432)
#shape = (128,256,2432)
chunk_shape = (128,128,128)
dtype = 'uint8'


def main(output_filename):

    f = chunkedBinary(fname, shape, chunk_shape, dtype)

    out = h5py.File(output_filename,'w')
    d = out.create_dataset('/main', shape, chunks=chunk_shape, compression='gzip', compression_opts=6)

    bounds = f.compute_bounds()
    i = 1
    num_bounds = len(bounds)

    start = timeit.default_timer()
    for beg, end in bounds:

        print "chunk %d of %d" % (i, num_bounds)


        d[ beg[0]:end[0],
           beg[1]:end[1],
           beg[2]:end[2]] = f.read_next_chunk()

        i += 1

    end = timeit.default_timer()
    print "transferred in %f seconds" % (end-start)
    out.close()
    f.close()

if __name__ == '__main__':

    main('cbf_bigger_test.h5')

