#!/usr/bin/env python

import numpy as np
import h5py
import tifffile
import os, glob
import subprocess, timeit
import sliding_window
#hopefully the only function we'll need from here
from emirt.volume_util import crop

#=======================================================================
#CONSTANTS
DTYPE='float32'
DISPLAY_DTYPE = 'float32'#'uint8'

omnify_directory = xxx
omnify_global_vars = xxx

#Default watershed_params
default_gws_high = 0.9999
default_gws_low = 0.3
default_gws_dust = 100
default_gws_dust_low = 0.3

error_curve = xxx
error_metric_dirname = 'error_metrics/'
#=======================================================================

def read_bin(fname):
	d = np.fromfile(fname, dtype='double')
	shape = np.fromfile(fname + '.size', dtype='uint32')[::-1]
	return d.reshape(shape)

def write_bin(data, fname):
	shape = np.array(data.shape).astype('uint32')
	data.tofile(fname)
	shape.tofile(fname + '.size')

def read_s3_aff(fname):
	data = np.fromfile(fname, dtype='float32')
	return data.reshape((3,128,128,128))

def read_tif(fname):
	return tifffile.imread(fname)

def write_tif(data, fname):
	tifffile.imsave(fname, data)

def read_h5(fname):
	f = h5py.File(fname)
	d = f['/main'].value
	f.close()
	return d

def write_h5(data, fname):

	if os.path.exists(fname):
		os.remove(fname)

	f = h5py.File(fname)
	f.create_dataset('/main',data=data)
	f.close()

def write_chunked_h5(data, fname):

	if os.path.exists(fname):
		os.remove(fname)

	f = h5py.File(fname)
	f.create_dataset('/main',data=data, chunks=(128,128,128), compression='gzip',compression_opts=4)
	f.close()

def read_file(fname):
	if '.h5' in fname:
		return read_h5(fname)
	elif '.tif' in fname:
		return read_tif(fname)
	elif os.path.isfile( fname + '.size'):
		return read_bin(fname)
	else:
		print "Unsupported file? (no .tif or .h5)"

def write_file(data, fname):
	if '.h5' in fname:
		write_h5(data, fname)
	elif '.tif' in fname:
		write_tif(data, fname)
	else:
		print "Unsupported file? (no .tif or .h5)"

def write_file_for_display(data, fname):
    '''
    Breaks the volume into 3d chunks, and saves each as
    the display dtype. Normalizes to [0,255] if saving as
    uint8
    '''
    to_save = np.copy(data)

    if len(to_save.shape) == 4:

        for i in range(to_save.shape[0]):

            to_save_vol = to_save[i,:,:,:]

            if DISPLAY_DTYPE == 'uint8':
              # to_save_vol = (to_save_vol - np.min(to_save_vol))
              to_save_vol = (to_save_vol / float(np.max(to_save_vol))) * 255

            to_save_vol = to_save_vol.astype(DISPLAY_DTYPE)

            root, ext = os.path.splitext(fname)
            new_fname = "%s_%d%s" % (root, i, ext)

            write_file(to_save_vol, new_fname)

    else:
        write_file(to_save, fname)

def tif2h5(fname, new_fname):
	write_h5( read_tif(fname), new_fname )

def h52tif(fname, new_fname):
	write_tif( read_h5(fname), new_fname )

def cls_err(output, label, verbose=True, return_vol=False):
	assert output.shape == label.shape

	error_vol = (output > 0.5) != (label > 0.5)
	errors = np.count_nonzero( error_vol )

	cls_err = errors / float(output.size)

	if verbose:
		print errors
		print output.size
		print cls_err

	if return_vol:
		return cls_err, error_vol
	else:
		return cls_err

def sq_err(output, label, verbose=True, return_vol=False):
	assert output.shape == label.shape

	error_vol = np.square( output - label )
	err = np.sum( error_vol )
	print np.max( error_vol )

	if verbose:
		print err

	if return_vol:
		return err, error_vol
	else:
		return err

def abs_err(output, label, verbose=True, return_vol=False):
	assert output.shape == label.shape

	error_vol = np.abs( output - label )
	err = np.sum( error_vol )

	if verbose:
		print err

	if return_vol:
		return err, error_vol
	else:
		return err

def unif_transformation(d, memopt=True):

	dflat = d.flatten()

	if memopt:

		print "Double Sort"
		indices = np.argsort(np.argsort(dflat))

	else:

		print "Sort"
		dsort = np.argsort(dflat)

		print "Inverting sort"
		indices = np.empty((d.size,),dtype=np.uint32)
		indices[dsort] = range(d.size)

	print "Generating values"
	values = np.linspace(0,1,d.size)

	print "Making new array"
	newd = values[indices].reshape(d.shape)

	return newd.astype(DTYPE)

def gaussian_transformation(d, mu, sigma):

	dflat = d.flatten()
	print "Double sort"
	indices = np.argsort(np.argsort(dflat))

	print "Generating values"
	values = np.random.normal(mu,sigma,d.size)
	values = np.sort(values)

	print "Making new array"
	newd = values[indices].reshape(d.shape)

	return newd.astype(DTYPE)


def expand_volume( data, factor ):
	#assert data.dtype == np.uint32
	assert data.dtype == np.float64

	start = timeit.default_timer()
	res = sliding_window.expand_volume( data, factor )
	end = timeit.default_timer()

	print "completed in %f seconds" % (end-start)
	return res


def avgfilter( data, kernel_size, stride, dtype=DTYPE ):

	start = timeit.default_timer()
	res   = sliding_window.sliding_window( data, np.mean, kernel_size, stride )
	end   = timeit.default_timer()

	print "completed in %f seconds" % (end-start)
	return res


def minfilter( data, filter_size, stride, dtype=DTYPE ):

	start = timeit.default_timer()
	res   = sliding_window.sliding_window( data, np.min, filter_size, stride)
	end   = timeit.default_timer()

	print "completed in %f seconds" % (end-start)
	return res

def slow_minfilter( data, filter_size, stride, dtype=DTYPE ):

	def min_func( d ):
		return np.min( d ).astype(dtype)

	start = timeit.default_timer()
	res = sliding_window.slow_sliding_window( data, min_func, filter_size, stride )
	end = timeit.default_timer()

	print "completed in %f seconds" % (end-start)
	return res

def minpool( data, filter_size, dtype=DTYPE ):
	return minfilter( data, filter_size, filter_size, dtype )

def slow_minpool( data, filter_size, dtype=DTYPE ):
	return slow_minfilter( data, filter_size, filter_size, dtype )

def center_subs( data, filter_size, stride=None ):

	if stride is None:
		stride = filter_size

	#bias to the right/larger index pixel
	center_index = np.array(filter_size) / 2

	def center_func( d ):
		return d[
			center_index[0],
			center_index[1],
			center_index[2]]

	start = timeit.default_timer()
	res = sliding_window.sliding_window( data, center_func, filter_size, stride )
	end = timeit.default_timer()

	print "completed in %f seconds" % (end-start)
	return res

def consolidate_ids(vol):
	unique_ids = np.unique(vol)

	for i in range(len(unique_ids)):
		vol[ vol == unique_ids[i]] = i

	return vol

def rgb_to_id(vol, consolidate=True):

	# vol = vol.astype('uint32')
	vol = (vol[0,:,:,:] +
		   vol[1,:,:,:] * 256 +
		   vol[2,:,:,:] * 256 ** 2)

	if consolidate:
		print "Consolidating ids"
		vol = consolidate_ids(vol)

	return vol

def fill_omnify_vars(dir_name, aff_fname, channel_fname,
	gws_high = default_gws_high,
	gws_low = default_gws_low,
	gws_dust = default_gws_dust,
	gws_dust_low = default_gws_dust_low):

    with open(omnify_global_vars) as f:

    	lines = f.readlines()
    	f.close()

    lines = replace_line( lines, 'gtmp = ', dir_name,
    					add_quotations=True )
    lines = replace_line( lines, "gchann_file = ", channel_fname,
    					add_quotations=True )
    lines = replace_line( lines, "gaffin_file = ", aff_fname,
    					add_quotations=True )
    lines = replace_line( lines, "gws_high = ", gws_high )
    lines = replace_line( lines, "gws_low = ", gws_low )
    lines = replace_line( lines, "gws_dust = ", gws_dust )
    lines = replace_line( lines, "gws_dust_low = ", gws_dust_low )

    with open(omnify_global_vars, 'w+') as f:

    	f.writelines(lines)
    	f.close()

def replace_line( lines, expression, value, add_quotations=False ):

	value = str(value)

	if add_quotations:
		value = '"{}"'.format(value)

	for i in range(len(lines)):

		if expression in lines[i]:

			print "Replacing: {}".format(lines[i])
			print "with {}".format(expression + value)

			lines[i] = expression + value + '\n'

	return lines

def run_omnification():

	stack = pushd( omnify_directory )

	subprocess.call( ['python','main.py'] )

	stack = popd( stack )

def pushd(destination, stack = []):
	stack.append( os.getcwd() )
	os.chdir(destination)

	return stack

def popd(stack):
	destination = stack.pop()
	os.chdir(destination)

	return stack

def makedir(dir_name):

	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def run_error_curve(ws_filename, label_filename, out_tail,
	thr_low=0.5, thr_high=1, num_thresholds=51):
    '''Runs segerror.error_curve, saves the results to an
    error_metrics subdirectory, and returns the maximum full RFS'''

    makedir(error_metric_dirname)

    outname = error_metric_dirname + out_tail

    thr_low , thr_high, num_thresholds = map(str,(thr_low,thr_high,num_thresholds))

    command = [error_curve,
    		ws_filename,
    		label_filename,
    		thr_low, thr_high, num_thresholds,
    		'-outname', outname]

    subprocess.call(command)

    RFS_fname = "%s_RFS_%s_%s_%s" % (outname, thr_low, thr_high, num_thresholds)

    print out_tail
    print "Max Score: {}".format(max_full_score(RFS_fname))

def max_full_score(score_fname):

	#Reading input file
	shape = np.fromfile(score_fname + '.shape', dtype='int64')
	values = np.fromfile(score_fname, dtype='float64').reshape(shape)

	return np.max(values[:,0]), np.argmax(values[:,0])


def zeropad(dataset, new_size, pick_right=False):

	size_diffs = np.array(new_size) - np.array(dataset.shape)

	assert all( size_diffs >= 0 )

	starting_index = (size_diffs / 2).astype('uint32')

	#boolean mult works like 'and' operation:
	# increment the starting index if there's a choice to
	# make AND the user wants us to
	increment_index = ( (size_diffs % 2 != 0) * pick_right ).astype('uint32')
	starting_index += increment_index

	result = np.zeros(new_size, dtype=dataset.dtype)

	# Placing original data into the center of the result
	zmin = starting_index[0]; zmax = starting_index[0] + dataset.shape[0]
	ymin = starting_index[1]; ymax = starting_index[1] + dataset.shape[1]
	xmin = starting_index[2]; xmax = starting_index[2] + dataset.shape[2]

	result[zmin:zmax,ymin:ymax, xmin:xmax] = dataset

	return result

def make_it_so(dataset, new_size, pick_right=False):
	'''Engage'''
	#max along each dimension
	new_size = np.array(new_size, dtype=np.uint32)
	temp_size = np.max( (np.array(new_size),np.array(dataset.shape)),0 )
	# temp_size = np.array((500,500,500))

	temp_dataset = zeropad( dataset, temp_size, pick_right )
	return crop( temp_dataset, new_size, pick_right=pick_right )


def find_bounds( total_shape, chunk_dimensions, step_dimensions=None, force_final=True ):

	if step_dimensions is None:
		step_dimensions = chunk_dimensions

	z_bounds = bounds_1d( total_shape[0], chunk_dimensions[0], step_dimensions[0], force_final )
	y_bounds = bounds_1d( total_shape[1], chunk_dimensions[1], step_dimensions[1], force_final )
	x_bounds = bounds_1d( total_shape[2], chunk_dimensions[2], step_dimensions[2], force_final )

	bounds = []
	for z in z_bounds:
		for y in y_bounds:
			for x in x_bounds:
				bounds.append(
					(
					np.array((z[0],y[0],x[0])),
					np.array((z[1],y[1],x[1]))
					)
				)

	return bounds

def bounds_1d( vol_width, chunk_width, step_width, force_final=True ):

	bounds = []

	beginning = 0
	ending = chunk_width

	while ending < vol_width:

		bounds.append(
			( beginning, ending )
			)

		beginning += step_width
		ending    += step_width

	#last bound
	if ending == vol_width or force_final:
		bounds.append(
			( vol_width - chunk_width, vol_width )
			)

	return bounds

def tuple_mult( tup, other ):
	if isinstance(other, tuple):
		#elementwise mult
		return map( lambda x: x[0] * x[1], zip(tup, other) )
	else:
		#scalar mult
		return map( lambda x: other * x, tup )

def tuple_add( tup, other ):
	if isinstance(other, tuple):
		#elementwise add
		return map( lambda x: x[0] + x[1], zip(tup, other) )
	else:
		#scalar add
		return map( lambda x: other + x, tup )

def tuple_sub( tup, other ):
	if isinstance(other, tuple):
		#elementwise mult
		return map( lambda x: x[0] - x[1], zip(tup, other) )
	else:
		#scalar mult
		return map( lambda x: x - other, tup )

def tuple_div( tup, other ):
	if isinstance(other, tuple):
		#elementwise mult
		return map( lambda x: x[0] / x[1], zip(tup, other) )
	else:
		#scalar mult
		return map( lambda x: other / x, tup )
