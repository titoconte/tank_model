#Tank model
This project goal is develop a numeric modelling of a rotational fluid inside a cylinder to simulate the physical model.The scrips are in Fortran.

## Scripts

The main script is called model_polito.f. It's based in model_polito.m, a MatLab script witch are slower to process the data as the model demend. The other fortran scrips are used as model's functions.

## Input File
The param.txt is a file to read the model inputs. In orther:
* g=gravity (m/s²)
* dt = time step (seconds)
* ntime = number of time steps
* nr=number of ratio cells
* nth = number of theta cells
* difz = vertical diffusion
* difT= turbullent diffusion
* rfric = friction coeeficient 
* fsave = plot frequency (1/s)
* rpm = rotation per minute 
* rho = water density (kg/m³)
* Raio = tank ratio (meters)

The dep.dep file is the bathymetry file 
