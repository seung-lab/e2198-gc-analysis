# e2198-gc-analysis
Connectomics analysis of Ganglion cells traced in e2198.

Please note two submodules are included in this repository. They constitute code bases
used by several subsequent analyses. To download all relavent data and code from the submodules,
use git command `git submodule update --init`.

When appropriate, additional documentations are available as the README.md file in specific subfolders.


## EM data, cell reconstructions / `e2198_reconstruction`
This submodule contains many of the initial steps, and more general facilities that are needed by almost all subsequent analyses. Especially code to create the initial cell_info structure is availble here. (see for additional details)

Importantly, the final cell type assignments as a result of the hierarchical clustering are in the cell_info_typedef_gc.m file contained in this repository.

Routines that carry special interests in this submodule include:
* Computational flattening and downsampling of the EM volume


## Calcium imaging data and analyses / `e2198_Ca_imaging`

Raw calcium imaging data, and all analyses pertaining to this data are in this submodule.

Some earlier hierarchical clustering efforts are also in this folder.

Routines that carry special interests in this submodule include:
* Contact analysis


## Segmentation data
TO BE UPLOADED.


## Data of special interests, by directory name:

### cell_features
Computed features for each cell. All features used for hierarchical clustering are here.
Note mapping of cell IDs and their assigned cell types are available in `e2198_reconstruction/cell_info_typedef_gc.m`.

### ganglion_skeletons
Skeletized RGC dendritic arbors.


## Code for specific analyses, by directory name:

### arbor_asymmetry
Code for quantifying the asymmetry of dendritic arbors. Earlier versions of the The same 

### n764
Convolutional network parameters. This is the network for segmenting somata for the soma size analyeses.
This was also the network used to generate the original segmentation used in Eyewire.org

### plotting
R code used to create some figures in the manuscript.

### skeletonization	
Code for automated skeletonization.

### soma_measurement
Code to compute soma size.


