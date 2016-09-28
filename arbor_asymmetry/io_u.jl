#!/usr/bin/env julia

__precompile__()

module io_u

using HDF5
using MAT

export write_bin, read_bin
export write_chunked_h5, read_h5, write_h5, write_chunked_h5
export write_map_file, read_map_file

#using MAT
function read_mat_var( filename, varname )

  f = matopen(filename)
  res = read(f, varname)
  close(f)

  res
end

#using H5 here
function read_h5( filename, read_whole_dataset=true, h5_dset_name="/main" )
  #Assumed to be an h5 file for now
  if read_whole_dataset
    d = h5read( filename, h5_dset_name );
  else
    f = h5open( filename );
    d = f[ h5_dset_name ];
  end

  return d;
end


function write_h5( dset, filename, h5_dset_name="/main" )
  if isfile(filename)
    rm(filename)
  end

  HDF5.h5write(filename, "/main", dset);
end


function write_chunked_h5( dset, output_filename, h5_dset_name="/main", chunk_shape=(128,128,128) )

  h5open( output_filename, "w" ) do f

    d = d_create( f, h5_dset_name, datatype(typeof(dset[1,1,1])), dataspace(size(dset)), "chunk", chunk_shape, "compress", 4 );
    d[:,:,:] = dset;

  end
end


"""
    write_map_file( output_filename, dicts... )
  Take an iterable of dictionaries which all have the same
  keys. Write a file in which each line takes the form

  key;val1;val2... for the number of dictionaries
"""
function write_map_file( output_filename, dicts... )

  open(output_filename, "w+") do f
    if length(dicts) > 0
    for k in keys(dicts[1])

      vals = ["$(d[k]);" for d in dicts ];
      write(f, "$k;$(vals...)\n" )

    end #for k
    end #if length
  end #open(fname)
end


"""

    read_map_file( input_filename, num_columns, sep=";" )

  Reads in map files written by write_map_file. Returns the first
  num_columns dicts stored within the file, but will break if you specify
  too many dicts to read.
"""
function read_map_file( input_filename, num_columns, sep=";" )

  dicts = [Dict() for i=1:num_columns];

  open(input_filename) do f

    for ln in eachline(f)

      fields = map(x -> eval(parse(x)), split(ln, sep));

      key = fields[1]

      for d in 1:length(dicts)
        dicts[d][key] = fields[1+d];
      end

    end#for ln

  end#do f

  dicts
end


function write_list_file( output_filename, iterable )

  open(output_filename, "w+") do f

    for v in iterable
      write(f, "$(v) ")
    end

  end
end

#module end
end
