#!/usr/people/nturner/julia/julia

import e2198_u
import io_u

#----------------------------
percentage_threshold = 0.7;
pixel_resolution = [16 16 25] * 10; # x10 for downsampling

#Whether to read the entire dataset into memory
read_whole_dataset = true;
dtype = UInt32

#Size of the 'inspection window'
# I explicitly assume that each cell is less than
# 2 * each of these dimensions in downsampled voxels
window_radius = [400,400,400];
#----------------------------

segmentation_filename = ARGS[1];
coord_filename = ARGS[2];
output_filename = ARGS[3];

println("Reading data...");
d = io_u.read_data( segmentation_filename );

println("Reading coordinate file...")
cell_id_to_coord, cell_id_to_type  = e2198_u.read_coordinate_file( coord_filename ); #map from id to coordinates


cell_id_to_soma_id = Dict{Int,dtype}();
cell_id_to_estvol  = Dict{Int,Float64}();
cell_id_to_axes    = Dict{Int,Int}();
cell_id_to_nv      = Dict{Int,Int64}();

for (cell_id, coord) in cell_id_to_coord

  println(cell_id);
  cell_id_to_soma_id[cell_id] = d[ coord... ];

  @time (cell_id_to_axes[cell_id], cell_id_to_nv[cell_id] 
          = e2198_u.extract_ellaxes_in_window( d, coord ) );

  cell_id_to_estvol = e2198_u.ellipse_volume( cell_id_to_axes[cell_id] );

end

io_u.write_map_file( output_filename, soma_id_to_type, soma_id_to_radius, soma_id_to_coord, soma_id_to_nv );

