MeshOpt is copyright (C) 2014-2015 J. Moore, and is distributed under the terms
of the GNU General Public License, Version 2 or later.

MeshOpt generates arbitrarily high order curvilinear meshes from GMSH input and
a reference geometry using the medhods developed by Abel Gargallo-Peiro and 
Xevi Roca in:

"Optimization of a regularized distortion measure to generate curved high-order unstructured tetrahedral meshes", International Journal for Numerical Methods in Engineering


#### 1.0 GENERAL

Currently, the only supported geometry format is *.stp. The geometry used to 
generate the GMSH mesh must be composed of closed surfaces. 

A more comprehensive guide on how to use MeshOpt can be found in the wiki.


#### 2.0 INSTALLATION
MeshOpt is distributed as a binary for 64-bit Linux. No compilation is required.

#### 3.0 USAGE
Several demo cases are included in the directory Demo. For each case, the 
following are required to be in the case directory:

a) A valid GMSH file. Make sure to run 'gmsh -check yourmesh.msh' to ensure 
   there are no errors.

2) A .stp or .step file of your geometry, containing no naked edges 
   (i.e. must be composed of closed surfaces).

3) A configuration file named MeshOpt.config. This file specifies the mesh 
   order to generate, boundary layer thickness, minimum acceptible mesh quality,
   etc. Examples are provided in the Demo cases. 

MeshOpt is run by executing the binary Meshopt in the case directory. MeshOpt 
does not accept any command line arguments; everything is specified in the 
configuration file.

#### 4.0 BUGS/SUPPORT
Please feel free to contact me at johnpmooreiv@gmail.com with any bug reports 
or questions you may have.

