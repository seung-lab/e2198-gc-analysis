#!/usr/bin/env python

import numpy as np
import u

def prep_chunked_file( chunked_file_obj, processing_fn,
		beginning, ending,
		input_chunk_shape,
		processing_filter, processing_stride ):

	print "Reading data"
	data = chunked_file_obj

	print "Defining patch bounds"
	input_bounds, output_bounds = define_bounds( beginning, ending, input_chunk_shape,
						processing_filter, processing_stride )

	output_vol_shape = output_shape( ending - beginning, processing_filter, processing_stride )
	output = np.empty(output_vol_shape, dtype=data.dtype)

	assert len(input_bounds) == len(output_bounds)

	for i in range(len(input_bounds)):
		print "%d of %d" % (i+1, len(input_bounds))
		input_beg,  input_end  = input_bounds[i]
		output_beg, output_end = output_bounds[i]

		input_beg, input_end = beginning + input_beg, beginning + input_end

		input_window = data[
			input_beg[0]:input_end[0],
			input_beg[1]:input_end[1],
			input_beg[2]:input_end[2]]

		output[
			output_beg[0]:output_end[0],
			output_beg[1]:output_end[1],
			output_beg[2]:output_end[2]] = processing_fn(input_window)

	return output

def define_bounds( beginning, ending, input_chunk_shape, 
	processing_filter, processing_stride ):

	input_vol_shape = ending - beginning

	output_chunk_shape = output_shape( input_chunk_shape, processing_filter, processing_stride )
	output_vol_shape   = output_shape( input_vol_shape,   processing_filter, processing_stride ).astype('uint32')

	# If the filter stride doesn't evenly divide the input, we need to remove some of the input
	truncated_input = truncate_input_shape( output_vol_shape, processing_filter, processing_stride )

	assert input_shape_divides_evenly( output_chunk_shape )
	output_chunk_shape = output_chunk_shape.astype('uint32')

	input_step = input_chunk_shape - processing_filter + 1
	input_bounds  = u.find_bounds( truncated_input, input_chunk_shape, input_step )
	output_bounds = u.find_bounds( output_vol_shape, output_chunk_shape, output_chunk_shape )

	return input_bounds, output_bounds
	

def output_shape( data_shape, filter_size, stride ):

	shape_w_o_stride = (
		data_shape[0] - filter_size[0] + 1,
		data_shape[1] - filter_size[1] + 1,
		data_shape[2] - filter_size[2] + 1)

	return shape_w_o_stride / np.array(stride).astype('float64')


def truncate_input_shape( output_vol_shape, filter_size, stride ):
	return output_vol_shape * stride + filter_size - 1


def input_shape_divides_evenly( output_chunk_shape ):
	return np.all((output_chunk_shape % 1) == 0)


def test_input_bounds(input_bounds):

	for bound in input_bounds:
		beg, end = bound

		if np.any((np.array(end) - processing_filter + 1) % 9 != 0):
			print bound
