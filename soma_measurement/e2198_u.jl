#!/usr/bin/env julia

__precompile__()

# println("NOT IN WORKING CONDITION")
module e2198_u

import io_u

export read_coordinate_file, read_gc_type_file

#----------------------------
#coord file
contains_coords_regexp = r"\[.*\]"
delim = ';'

#gc type file
type_line_regexp = r"^struct.*";
id_list_regexp   = r"\[.*\]";

#downsampling parameters
downsampling_filter = [18, 18, 18];
downsampling_stride = [10, 10, 10];
window_radius = [400,400,400];


current_cell_info = (
"/usr/people/nturner/research/retina/e2198_Ca_imaging/data/cell_info_clustering.20160623warpedSomaCorrection.mat"
)
current_classification = (
"/usr/people/nturner/Dropbox/seunglab/plotting_retina_figures/final_final_classification.csv"
)
skel_file_directory = (
"/usr/people/nturner/seungmount/research/Alex/e2198/Cells_done"
)
#----------------------------

function read_soma_coords( cell_info = current_cell_info )

  cell_info = io_u.read_mat_var( cell_info, "cell_info" );

  ids = Int[cell_info["cell_id"]...];
  soma_coords = cell_info["soma_coords_warped_mip2_zscaled"];

  Dict{Int,Any}([ ids[i] => soma_coords[i] for i in 1:length(ids) ])
end

function read_classification_csv( filename = current_classification )

  type_to_ids = Dict{AbstractString,Vector{Int}}();

  open(filename) do f
    for (i,ln) in enumerate(eachline(f))

      if i == 1 continue end

      #hard-coding the delim here since
      # I don't want this to change by my
      # owh whims for this file
      fields = split(ln, ',');

      cell_type = replace(fields[1], "\"","");
      cell_id = parse(Int, fields[2]);

      if !haskey(type_to_ids, cell_type)
        type_to_ids[cell_type] = Vector();
      end

      push!(type_to_ids[cell_type], cell_id);
    end
  end

  type_to_ids
end

#NOTE currently returns the coordinates as Float64's
function read_warped_points( cell_id )

  filename = "$(skel_file_directory)/skel_$(cell_id).mat"

  io_u.read_mat_var( filename, "p" )[:,[3,2,1]] #flipping to XYZ
end

function read_skeleton_points( cell_id )

  filename = "$(skel_file_directory)/skel_$(cell_id).mat"

  io_u.read_mat_var( filename, "n" )[:,[3,2,1]] #flipping to XYZ
end

function read_candidate_coord_file( filename, coord_index )

  id_to_coord = Dict{Int,Array{Int64,1}}();

  open(filename) do f
    for ln in eachline(f)

      fields = split(ln,delim);

      cell_id = parse(Int,fields[1]);

      coord_str = fields[coord_index]
      coords  = map(x->parse(Int,x), split(coord_str[2:end-1],','));

      id_to_coord[ cell_id ] = coords;
    end #for ln
  end #open filename


  return id_to_coord;
end


function read_coordinate_file( filename )

  id_to_coord = Dict{Int,Array{Int64,1}}();
  id_to_type = Dict{Int,AbstractString}();

  open(filename) do f
    for ln in eachline(f)

      if ! ismatch(contains_coords_regexp,ln) continue end

      fields = split(ln, delim);

      cell_type = fields[1];
      cell_id   = parse(Int,fields[2]);
      coords    = map(x->parse(Int,x), split(fields[3][2:end-2],','));

      id_to_coord[cell_id] = coords;
      id_to_type[cell_id]  = cell_type;

    end #for ln
  end #open

  return id_to_coord, id_to_type;
end

function read_gc_type_file( filename )

  id_to_type = Dict{Int, AbstractString}();

  open( filename ) do f
    for ln in eachline(f)

      if ! ismatch(type_line_regexp, ln) continue end

      typename = split(ln,"'")[4];

      cell_ids_str = match(id_list_regexp, ln).match[2:end-1]; #trimming brackets

      cell_ids = map(x->parse(Int,x), split(cell_ids_str," "));

      for id in cell_ids
        id_to_type[id] = typename;
      end

    end#for ln
  end#open

  return id_to_type;
end


"""
    find_COM_in_window( dset, coord )
"""
function find_COM_in_window( dset, coord )

  window, rel_coord = get_inspection_window( dset, coord );
  rel_COM = center_of_mass( window, rel_coord );

  return rel_COM + (coord - rel_coord); #adjusting for window selection

end


"""
    extract_ellaxes_in_window( dset, coord )
"""
function extract_ellaxes_in_window( dset, coord::Array{Int64,1},
  pixel_resolution, percentile_threshold )

  window, rel_coord = get_inspection_window( dset, coord );
  return extract_ellaxes_by_coord( window, rel_coord,
    pixel_resolution, percentile_threshold );
end


"""
    extract_radius_in_window( dset, coord )
"""
function extract_radius_in_window( dset, coord::Array{Int64,1},
  pixel_resolution, percentile_threshold )

  window, rel_coord = get_inspection_window( dset, coord );
  return extract_radius_by_coord( window, rel_coord,
    pixel_resolution, percentile_threshold );
end


#Mid level functions
"""
    get_inspection_window( dset, coord )
"""
function get_inspection_window( dset, coord::Array{Int64,1} )

  beg_candidate = coord - window_radius;
  end_candidate = coord + window_radius;

  beg_i = max( [1,1,1], beg_candidate );
  end_i = min( [size(dset)...], end_candidate );

  rel_coord = coord - beg_i;

  window = dset[
    beg_i[1]:end_i[1],
    beg_i[2]:end_i[2],
    beg_i[3]:end_i[3]];

  return window, rel_coord;
end


"""
    extract_ellvol_by_coord( seg_vol, coord )

  Computes the volume from the ellipse formed by the 70th percentile
  distances along each
"""
function extract_ellaxes_by_coord{T}( seg_vol::Array{T,3}, coord::Array{Int64,1},
  pixel_resolution, percentile_threshold )

  seg_id = seg_vol[ coord... ];

  seg_points = findn( seg_vol .== seg_id );

  axes_distances = axes_dist_from_COM( seg_points );

  #dummy code for now
  percentile_index = Int( ceil( percentile_threshold * size(axes_distances,1) ) );

  axes = axes_distances[percentile_index];

  return axes, size(a, 1);
end


"""
    extract_radius_by_coord( seg_vol, coord )
"""
function extract_radius_by_coord{T}( seg_vol::Array{T,3}, coord::Array{Int64,1},
  pixel_resolution, percentile_threshold )

  seg_id = seg_vol[ coord... ];

  seg_points = findn( seg_vol .== seg_id );

  radii = radii_from_COM( seg_points );

  percentage_index = Int( ceil( percentile_threshold * length(radii) ) );

  return radii[ percentage_index ], length(radii);
end


"""
    radii_from_COM( indices )

  Find the sorted distances for each point from COM in 3D
"""
function radii_from_COM( indices, pixel_resolution )

  spatial_diffs = diffs_from_COM( indices );

  radii = sqrt(sum( spatial_diffs .^ 2, 2));

  return sort(vec(radii));
end

"""
    axis_dist_from_COM( indices )

  Compute absolute distances from COM along each x, y, and x axis
  from the raw indices
"""
function axis_dist_from_COM( indices, pixel_resolution )

  spatial_diffs = diffs_from_COM( indices );

  axes_distances = abs( spatial_diffs );

  sort!(axes_distances, 1);

  return axes_distances;
end


"""
    axis_diffs_from_COM( indices )
  Find distances of each point from COM along each axis,
  sort the values for each axis
"""
function diffs_from_COM( indices, pixel_resolution )

  indices_a = [ indices[1] indices[2] indices[3] ];

  center_of_mass = [ mean(indices[i]) for i=1:length(indices) ];

  spatial_diffs = (indices_a .- center_of_mass' ) .* pixel_resolution;

  return spatial_diffs;
end


"""
    center_of_mass( seg_vol, coord )
"""
function center_of_mass{T}( seg_vol::Array{T,3}, coord::Array{Int64,1},
  round_result=true )

  seg_id = seg_vol[ coord... ];

  indices = findn( seg_vol .== seg_id );

  center_of_mass = Array{Float64}((length(indices),));

  for i in 1:length(indices)
    center_of_mass[i] = mean(indices[i]);
  end

  if round_result
    center_of_mass = convert(Array{Int64,1}, round(center_of_mass));
  end

  return center_of_mass;
end


"""
    ellipse_volume( axis_lengths )

  Compute the volume of the ellipse described
  by the axis lengths within each row
"""
function ellipse_volume( axis_lengths )
    return (0.75) * pi * prod(axis_lengths, 1)
end


function downsampled_coord_to_global( coord )
  return (coord .* downsampling_stride) + div(downsampling_filter,2);
end

function global_coord_to_downsampled( coord )
  return (coord - div(downsampling_filter,2)) ./ downsampling_stride;
end

function parse_typedef_file( filename )

  t = AbstractString("");

  open(filename) do f
    l = readlines(f)
    l = filter( line -> !ismatch( r"^%.*", line), l );
    t = join(l, " ");
    close(f)
  end

  t = replace(t, "...",  "");
  t = replace(t, r"\s+", " ");

  structs = split_into_structs( t );

  types = [ parse_typedef_struct(struct) for struct in structs ];

  types
end


function split_into_structs( contents )


  structs = [];
  start_index = 1;

  while true

    start_indices = search( contents, "struct(", start_index);

    if start(start_indices) == 0 break end

    after_start   = start(start_indices) + endof(start_indices);
    end_indices = search( contents, ");", after_start);

    push!( structs, contents[ after_start: (start(end_indices) - 1) ] );

    start_index = start(end_indices) + endof(end_indices);
  end

  structs
end

function parse_typedef_struct( contents )

  cell_type = Dict{AbstractString,Any}();

  fields = split( contents, "'" )[2:2:end];
  append!(fields, [split(contents,",")[end]] );

  @assert mod(length(fields), 2) == 0;

  i = 1;
  while i < length(fields)

    key = fields[i];

    if key != "cells"
      value = fields[i+1];
    else

      raw = fields[i+1];
      open_br  = search( raw, '[' );
      close_br = search( raw, ']' );

      value = map( x -> parse(Int,x), split( raw[open_br+1:close_br-1], " ") );

    end

    cell_type[key] = value;
    i += 2;

  end

  cell_type
end

end #module end
