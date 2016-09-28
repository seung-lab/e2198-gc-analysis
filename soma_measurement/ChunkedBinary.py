#!/usr/bin/env python

import numpy as np
import u, os

class ChunkedBinary( object ):

    def __init__(self, fname, shape, chunk_shape, dtype):
        self.fname = fname
        self.shape = shape
        self.chunk_shape = chunk_shape
        self.dtype = dtype

        self.chunk_size = reduce( lambda x,y: x*y, chunk_shape )

        self.f = open(fname, 'rb')
        self.chunk_id = 0

        self.bounds = None

        if dtype in ('uint8','int8'):
            self.dtype_size=1
        elif dtype in ('uint16','int16'):
            self.dtype_size=2
        elif dtype in ('uint32','float32'):
            self.dtype_size=4
        elif dtype in ('uint64','float64','double'):
            self.dtype_size=8


    def __getitem__(self, index):
        sl = self.__prep_slice__(index)

        if self.bounds is None:
            self.compute_bounds()

        to_read = []
        for i, bound in enumerate(self.bounds):

            overlap = self.__in_range__(bound, sl)
            if overlap:
              to_read.append( (i, overlap) )

        if len(to_read) > 0:
            return self.__read_and_map__( sl, to_read )
        else:
            return None #this should be changed at some point


    def __prep_slice__(self, sl):
        res = []
        for elem in sl:
          if isinstance(elem, slice):
            res.append((elem.start, elem.stop))
          else:
            res.append((elem, elem+1))

        return map(np.array,zip(*res))


    def compute_bounds(self, chunk_shape=None):

        if chunk_shape is None:
            chunk_shape = self.chunk_shape

        self.bounds = u.find_bounds( self.shape, chunk_shape )
        return self.bounds


    def __in_range__(self, bound, prepped_slice):
        b, e = bound
        ib, ie = prepped_slice

        if np.any(ie < b) or np.any(ib >= e): return False #also change this

        res_b = np.maximum(b, ib)
        res_e = np.minimum(e, ie)

        if np.any(res_b == res_e): return False #also change this

        rel_b = res_b - b
        rel_e = res_e - b

        out_b = res_b - ib
        out_e = res_e - ib

        return ((rel_b, rel_e), (out_b, out_e))


    def __read_and_map__(self, sl, to_read):

        output_shape = sl[1]-sl[0]
        output = np.empty(output_shape, dtype=self.dtype)

        for chunk_id, ranges in to_read:
            chunk = self.read_chunk(chunk_id)

            chunk_beg,  chunk_end  = ranges[0]
            output_beg, output_end = ranges[1]

            output[ output_beg[0]:output_end[0],
                    output_beg[1]:output_end[1],
                    output_beg[2]:output_end[2]
                  ] = chunk[
                    chunk_beg[0]:chunk_end[0],
                    chunk_beg[1]:chunk_end[1],
                    chunk_beg[2]:chunk_end[2]]

        return output


    def close(self):
        self.f.close()


    def read_next_chunk(self):

        self.chunk_id += 1
        return self.read_chunk(self.chunk_id)


    def read_chunk(self, chunk_id, verbose=False):
        
        if verbose:
            print "Reading chunk %d" % chunk_id
            
        file_pos = chunk_id * self.chunk_size * self.dtype_size 

        self.f.seek(file_pos, os.SEEK_SET)
        d = np.fromfile(self.f, dtype=self.dtype, count=self.chunk_size)
        return d.reshape(self.chunk_shape)

