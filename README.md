# e2198-gc-analysis
Connectomics analysis of Ganglion cells traced in e2198.

Please note two submodules are included in this repository. They constitute code bases
used by several subsequent analyses. To download all relavent data and code from the submodules,
use git command `git submodule update --init`.

When appropriate, additional documentations are available as the README.md file in specific subfolders.


## EM data, cell reconstructions / `e2198_reconstruction`
This submodule contains many of the initial steps, and more general facilities that are needed by almost all subsequent analyses. Especially code to create the initial cell_info structure is availble here. (see README in the submodule for additional details)

Importantly, the final cell type assignments as a result of the hierarchical clustering are in the cell_info_typedef_gc.m file contained in this repository.

Routines that carry special interests in this submodule include:
* Computational flattening and downsampling of the EM volume


## Calcium imaging data and analyses / `e2198_Ca_imaging`

Raw calcium imaging data, and all analyses pertaining to this data are in this submodule.

Some earlier hierarchical clustering efforts are also in this folder.

Relations amongst the EM data orientation, the visual stimuli directions, and the final polar coordinates are documented in this submodule.

Routines that carry special interests in this submodule include:
* Contact analysis


## Segmentation data
TO BE UPLOADED.


## Data of special interests, by directory name:

### Data
Computed features for each cell. All features used for hierarchical clustering are included here.

The cell_info structure array in 'cell_info.mat' is an aggregate of features from all
 other files, and for ganglion cells also includes the assigned cell type. (For more
 details on the cell_info structure, see the `e2198_reconstruction` submodule).

Note mapping of cell IDs and their assigned cell types are also available in `e2198_reconstruction/cell_info_typedef_gc.m`. 

### skeletons_GC
Skeletons of GC dendritic arbors.

p - Point cloud representation of cell segmentation (N_p x 3) [66 x 92 x 66 nm^3]

n - Coordinates of skeleton nodes (N_n x 3) [66 x 92 x 66 nm^3]

e - Connectivity of skeleton nodes represented by its row index in n (ex. [1, 2] : n(1,:) --> n(2,:)) (N_e x 2) 

bn - Coordinates of branch nodes (N_bn x 3) [66 x 92 x 66 nm^3]

rad - Distance to the boundary from each skeleton nodes (N_n x 1) [arbitrary unit]  
*These values cannot be converted to real unit due to anisotropy of voxel.

root - Coordinates of root nodes (N_rt x 3) [66 x 92 x 66 nm^3]

dest - Coordinates of destination nodes (N_dt x 3) [66 x 92 x 66 nm^3]


## Code for specific analyses, by directory name:

### analysis
MATLAB codes and data used for analysis. Run MATLAB codes within this folder to avoid dependency issues.

### arbor_segregation
Code for quantifying segregation of distribution of cell features. Run it within the folder "analysis" to avoid dependency issues.

### density conservation
Code used to test "density conservation" hypothesis. Run it within the folder "analysis" to avoid dependency issues.

### arbor_asymmetry
Code for quantifying the asymmetry of dendritic arbors.

### n764
Convolutional network parameters. This is the network for segmenting somata for the soma size analyeses.
This was also the network used to generate the original segmentation used in Eyewire.org

### skeletonization	
Code for automated skeletonization.

### soma_measurement
Code to compute soma size.


## Code for plotting figures

### plotting
Code used to create figures in the manuscript. For MATLAB codes, run it within the folder "analysis" to avoid dependency issues.
