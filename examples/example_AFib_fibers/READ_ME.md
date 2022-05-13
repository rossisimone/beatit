## Dependencies
To install and run this code, you need:
- VTK (https://vtk.org/)
- Libmesh (https://libmesh.github.io/)


## Install BeatIt
We use CMake to configure the BeatIt library.
```
git clone https://github.com/rossisimone/beatit.git

```

<under construction>

## Preprocessing steps

### Segmentation of images
Making sure that there are no overlapping structures. 

### Input file format, boundary sets and landmark points
The method requires the definition of common boundary sets (endocardium, epicardium, mitral valve ring, pulmonary veins) and two landmark points, one for the left atrial appendage and one for the fossa ovalis.


In case you need to assign your sidesets and soes not have an easy way to do so, you can use the brute_force.cpp? number of stl files + one volumetric exodus file

<Under construction>

### Output file format
The output file will be an exodus file format called `<input-file-name>_afib_fibers.e'.  

## Thresholds 

### Video Tutorial

<under construction>

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




