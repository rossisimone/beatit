## Dependencies
To install and run this code, you need:
- VTK (https://vtk.org/)
- Libmesh (https://libmesh.github.io/)
- Boost (https://www.boost.org/) 
- Eigen (https://eigen.tuxfamily.org/)
- PETSC (https://petsc.org/release/)

Instructions to install Beatit can be found at https://github.com/rossisimone/beatit/


## Preprocessing steps

### Segmentation of images
Making sure that there are no overlapping structures. 

### Input file format, boundary sets and landmark points
The method requires the definition of common boundary sets (endocardium, epicardium, mitral valve ring, pulmonary veins) and two landmark points, one for the left atrial appendage and one for the fossa ovalis. The landmark points should be added to the input file data.beat and is explained in the section "Input variables" below.

To assign the sidesets, you can follow the schematic illustrated in the diagram 

<img src="https://github.com/rossisimone/beatit/tree/master/examples/example_AFib_fibers/assigning_boundary_sets.png" width="128"/>

In case you need to assign your sidesets and does not have an easy way to do so, you can use the https://github.com/rossisimone/beatit/tree/master/examples/example_remap_region/brute_force.cpp by providing exodus surface files with the boundary set information and a volumetric mesh in exodus format. 



## Input variables 
The variable `common_left_trunk` should be set to `true` if the anatomy has one common left pulmonary vein. Say the landmark points for the LAA and FO are (-0.3, -2.94, -0.59) with radius 0.02 (chosen such that it encolses one mesh node) and (3.14, 3.61, 3.27) with radius 0.02 (chosen such that it encolses one mesh node), then this is how they should be set in the data.beat input file:
```
nodeset_x = '-0.3, 3.14'
nodeset_y = '-2.94, -3.61'
nodeset_z = '-0.59, 3.27'
nodeset_r = '0.02, 0.02'
```
In case you want your Bachmann's Bundle to cross the septum in a higher position, you can decrease the value set for `t_bachmann`.

### Mesh example
You can find a Mesh example here: https://osf.io/2q3up/?view_only=19abb61193ba4384a504392ca6769703
Note, the original mesh is the 03patient from http://doi.org/doi:10.18742/RDM01-289

The input file  `mesh`, `la_mesh` variables are set to be the input mesh names:
```
mesh=input_mesh.e
la_mesh=input_mesh.e
```
 --- under construction: include input files  ----

## Thresholds 

### Video Tutorial

 --- under construction ----

### Correspondence of threshold symbols in the paper and in the code

| Symbol paper |	Name of varible in the code |
|--------------|------------------------------|
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{(ENDO_EPI)}}}#gh-light-mode-only">
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{(ENDO_EPI)}}}#gh-dark-mode-only"> | hard coded | 
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{AF}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{AF}}}#gh-dark-mode-only"> | t_anterior_floor  |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{PF}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{PF}}}#gh-dark-mode-only"> |	t_floor |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{RA}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{RA}}}#gh-dark-mode-only"> |	t_antra_rsvp or t_antra_ripv|
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{AP}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{AP}}}#gh-dark-mode-only"> |	t_anterior_posterior |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{RIPV}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{RIPV}}}#gh-dark-mode-only"> |	t_septum_ripv_2 |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{LR}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{LR}}}#gh-dark-mode-only"> |	t_left_right |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{LA}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{LA}}}#gh-dark-mode-only"> |	t_antra_left |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{BB}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{BB}}}#gh-dark-mode-only"> |	t_bachmann |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{LLR}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{LLR}}}#gh-dark-mode-only"> |	t_left_lateral_ridge |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{AL}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{AL}}}#gh-dark-mode-only"> |	t_anterior_lateral |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{AR}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{AR}}}#gh-dark-mode-only"> |	t_anterior_top_1 |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{EL}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{EL}}}#gh-dark-mode-only"> |	t_lateral_u5 |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{LAL}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{LAL}}}#gh-dark-mode-only"> |	t_anterior_bottom_1 |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{LE}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{LE}}}#gh-dark-mode-only"> |	t_anterior_top_2 |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{FO}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{FO}}}#gh-dark-mode-only"> |	t_fo |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{S}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{S}}}#gh-dark-mode-only"> |	t_septum_ripv_1 |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{L}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{L}}}#gh-dark-mode-only"> |	t_lateral |
| <img src="https://render.githubusercontent.com/render/math?math={\\alpha_{\text{LAA}}#gh-light-mode-only"> 
<img src="https://render.githubusercontent.com/render/math?math={\color{white}\alpha_{\text{LAA}}}#gh-dark-mode-only"> |	t_laa |



## Output file format
The output file will be an Exodu0s II file called `input_mesh_afib_fibers.e'.  
