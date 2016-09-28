#!/usr/bin/env julia

__precompile__()

module asymmetry

import MAT
import io_u

#-----------------------------------------
voxel_dimensions = [16.5 16.5 23];
zscaled_voxel_dimensions = [16.5 16.5 16.5];

skel_file_node_key = "n"
soma_coord_key = "soma_coords_warped_mip2_zscaled" # => requires nodes to be rescaled in z
cell_id_regexp = r"skel_([0-9]+).mat"

#default inputs
#default_classification_file = xxx
#default_skeleton_directory = xxx
#default_cell_info_filename = xxx
default_output_filename = "asymmetry.csv"

#NOTE skeleton coordinates are ZYX, not XYZ (python import)
#NOTE soma_warped coordinates are YZX, and Z is pre-scaled by 23/16.5

#NOTE this is in mip2 warped coords (no z scaling)
#dir_2an = [-0.5209310049005197 0.6284252444851756]; old, no zscaling
dir_2an = [-0.46187405392971975 0.6987995575075379];
#-----------------------------------------

searchdir(path,key) = filter(x->contains(x,key), readdir(path))

"""

    read_classification_file( classification_filename, delim="," )
"""
function read_classification_file( classification_filename, delim="," )

  id_to_type = Dict{Int, AbstractString}();

  open(classification_filename) do f
    for (i,ln) in enumerate(eachline(f))
    
      if i == 1 continue end #header line
      
      fields = split(ln, delim);

      cell_type = replace(fields[1],'"',"");#removing " characters
      cell_id = parse(Int, fields[2]);

      id_to_type[ cell_id ] = cell_type;
    end #for ln
  end

  id_to_type
end 

"""

    read_nodes_from_mat( filename, key="$(soma_coord_key)" )

  Extracts the skeleton nodes (as a 2D array) from a .mat file
  from the given key.
"""
function read_nodes_from_mat( filename, key=skel_file_node_key )
  
  n = Array{Int,2}();

  MAT.matopen( filename ) do f
    n = MAT.read( f, key );
  end

  convert(Array{Int,2},n);
end

"""

    read_soma_coords( cell_info_filename )

  Create a mapping from cell_id to (currently warped and zscaled)
  soma coordinate from a .mat file containing a cell_info struct.
"""
function read_soma_coords( cell_info_filename, key=soma_coord_key )

  id_to_coord = Dict()

  MAT.matopen( cell_info_filename ) do f
    cell_info = MAT.read( f, "cell_info" );

    cell_ids = convert(Vector{Int},vec(cell_info["cell_id"]));
    soma_coords = cell_info[ key ];

    id_to_coord = [ cell_ids[i] => soma_coords[i] 
                         for i=1:length(cell_ids)];
  end

  id_to_coord;
end


function scale_z_to_orig( coords )
  z_factor = voxel_dimensions[2] / voxel_dimensions[3];
  scale_z(coords, z_factor);
end


function scale_z_to_y_scale( coords )
  z_factor = voxel_dimensions[3] / voxel_dimensions[2];
  scale_z(coords, z_factor);
end


function scale_z( coords, factor )
  multiplier = ones(size(coords,2),1);
  multiplier[end] = factor;
  return coords .* multiplier';
end


function extract_id( filename )
  int_str = match(cell_id_regexp, filename).captures[1]  
  parse(Int, int_str)
end


function find_node_centroid( nodes )
  mean(nodes,1)
end

function find_arbor_vector( centroid, soma_coord )
  (centroid' - soma_coord)[2:end]
end

function compute_asym_index( arbor_vector, nodes, soma_coord )

  nodes .-= soma_coord';
  #limiting to yz projection
  nodes = nodes[:,2:end]

  dot_products = nodes * arbor_vector;
  forward_nodes = count( x -> x > 0, dot_products );
  #backward_nodes = count( x -> x < 0, dot_products );

  forward_nodes / size(nodes,1);
end

function omni_spatial_to_volume_coord( spatial_coord, mip_change=0 )

  return (spatial_coord ./ voxel_dimensions) * (2^(-mip_change));

end

function find_skel_file( cell_id, skeleton_directory=default_skeleton_directory )
  searchdir( skeleton_directory, string(cell_id) )[1];
end


"""
Used to find average 2an direction
"""
function find_average_directions( cell_ids,
             skeleton_directory=default_skeleton_directory,
             cell_info_filename=default_cell_info_filename)

  num_cells = length(cell_ids);
  id_to_coords = read_soma_coords( cell_info_filename );

  unit_avs = Array{Float64,2}(num_cells, 2);

  for (i,cell_id) in enumerate(cell_ids)

    println(cell_id)

    skel_file = find_skel_file( cell_id, skeleton_directory );

    #importing soma coordinate
    #see main for some notes on this process
    warped_zscaled_soma = id_to_coords[ cell_id ];
    if typeof(warped_zscaled_soma) == Array{Union{},2} error("cell $(cell_id) not a ganglion cell?") end
    warped_zscaled_soma = warped_zscaled_soma[[3,1,2]];

    #indexing from zyx -> xyz
    nodes = read_nodes_from_mat("$skeleton_directory$skel_file")[:,[3,2,1]];
    nodes = scale_z_to_y_scale(nodes);

    skeleton_centroid = find_node_centroid( nodes );

    arbor_vector = find_arbor_vector( skeleton_centroid, warped_zscaled_soma );

    arbor_vector ./= norm(arbor_vector);

    unit_avs[i,:] = arbor_vector;

  end

  println(mean(unit_avs,1))
  mean(unit_avs,1)

end


function main( classification_file = default_classification_file,
                   skeleton_directory=default_skeleton_directory, 
	           cell_info_filename=default_cell_info_filename,
	           output_filename=default_output_filename )
  
  #all_skel_files = searchdir( skeleton_directory, ".mat" )

  id_to_type = read_classification_file( default_classification_file, "," );
  id_to_coords = read_soma_coords( cell_info_filename );

  id_to_avec = Dict{ Int, Vector{Float64} }();
  id_to_asym = Dict{ Int, Float64 }();
  id_to_centroid = Dict{ Int, Vector{Float64} }();
  id_to_soma = Dict{ Int, Vector{Float64} }();
  id_to_asym2an = Dict{ Int, Float64 }();
  id_to_typeasym = Dict{ Int, Float64}();

  num_per_type = Dict{ AbstractString, Int }();
  for cell_type in values(id_to_type)
    num_per_type[cell_type] = get(num_per_type,cell_type,0) + 1;
  end
  type_to_direction = [ cell_type => Float64[0,0] for cell_type in keys(num_per_type) ];

  #for f in all_skel_files #old
  # first loop to define average direction
  for cell_id in keys(id_to_type)

    f = find_skel_file( cell_id, skeleton_directory );

    println(cell_id)

    #indexing shifts from yzx to xyz (see note in constants section)
    # it also makes the changes the type to Vector{Float64}
    warped_zscaled_soma = id_to_coords[ cell_id ];
    if typeof(warped_zscaled_soma) == Array{Union{},2} continue end #non-ganglion cells
    warped_zscaled_soma = warped_zscaled_soma[[3,1,2]];
    
    #indexing from zyx->xyz (see note in constants section)
    nodes = read_nodes_from_mat("$skeleton_directory$f")[:,[3,2,1]];
    nodes = scale_z_to_y_scale(nodes);


    skeleton_centroid = find_node_centroid( nodes );

    arbor_vector = find_arbor_vector( skeleton_centroid, warped_zscaled_soma );

    asym_index = compute_asym_index( arbor_vector, nodes, warped_zscaled_soma );
    asym_2an   = compute_asym_index( dir_2an',     nodes, warped_zscaled_soma );

    id_to_avec[ cell_id ] = vec(scale_z_to_orig(arbor_vector));
    id_to_asym[ cell_id ] = asym_index;
    id_to_centroid[ cell_id ] = vec(scale_z_to_orig(skeleton_centroid));
    id_to_soma[ cell_id ] = vec(scale_z_to_orig(warped_zscaled_soma));
    id_to_asym2an[ cell_id ] = asym_2an;

    #adding arbor vector to type average so far
    cell_type = id_to_type[ cell_id ];
    type_to_direction[ cell_type ] += (arbor_vector ./ num_per_type[ cell_type ] );
  end

  println(type_to_direction);

  #second loop to compute index on the type-defined directions
  for cell_id in keys(id_to_type)

    #reload soma
    #indexing shifts from yzx to xyz (see note in constants section)
    # it also makes the changes the type to Vector{Float64}
    warped_zscaled_soma = id_to_coords[ cell_id ];
    if typeof(warped_zscaled_soma) == Array{Union{},2} continue end #non-ganglion cells
    warped_zscaled_soma = warped_zscaled_soma[[3,1,2]];

    #reload nodes
    f = find_skel_file( cell_id, skeleton_directory );
    #indexing from zyx->xyz (see note in constants section)
    nodes = read_nodes_from_mat("$skeleton_directory$f")[:,[3,2,1]];
    nodes = scale_z_to_y_scale(nodes);

    cell_type = id_to_type[ cell_id ];
    type_direction = type_to_direction[ cell_type ];

    type_asym_index = compute_asym_index( type_direction, nodes, warped_zscaled_soma );
    id_to_typeasym[ cell_id ] = type_asym_index;
  end

  io_u.write_map_file( output_filename, 
  	              id_to_avec, id_to_centroid, id_to_soma,
                      id_to_asym, id_to_asym2an, id_to_typeasym );

end

main()

end #module end
