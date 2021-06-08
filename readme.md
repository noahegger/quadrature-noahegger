# Program Goal

   The program has two major components: First, the user defines a box-shaped nuclear reactor with a detector placed a variable distance away from one of the x faces. The distance of the detector from the reactor is defined by the user (user defines minimum and maximum distance of the detector as well as the step interval distance). The program calculates the flux at the detector location, which starts at the user defined minimum location and moves one step at a time to the maximum location, while calculating flux at each step. 

   The second part of the program calculates the flux from the same box but with a hollow spherical region within the box, centered on the center of the box. The radius of the spherical region starts at a minimum value and steps to a maximum value taking flux measurements at each step. The detector is stationary at a set distance from the reactor while the hollow region radius increases. These minimum, maximum, and step radius values are provided by the user, as well as the detector distance. 

# Directions for Usage

   Navigate to the relevant directory through terminal. Make sure all `.f90` files are in the same directory. Type "make" into the terminal as this will compile the files and create an executable called `nuclear-reactor`. In terminal type "./nuclear-reactor" to run the program. You will be prompted to enter values for the depth, width, and height of the box, as well as the y-coordinate of the detector, the minimum and maximum distance of the detector from the reactor in addition to the step size, the number of intervals over which the integrations will be approximated (for Newton-Coates quadrature), and the number of sampling points (for the Monte Carlo quadrature). The program will write the flux values (for each detector distance using each integration method) into a file called `results_basic.dat`. You can visualize them by running the jupyter notebook `.ipynb` file.

   The program will then request another set of inputs from the user for the hollow reactor section of the program. Enter the position of the detector (stationary for this part of the program), the minumum and maximum radius of the spherical hollow region, and the step size of this radius. The fluxes will be measured for the box with a hollow region growing from the minimum radius to the maximum radius by the step size of the radius given by the user. The flux values will be written intp a file called `results_advanced.dat`. You can visualize the data by running the jupyter notebook `.ipynb` file.
    
# Contents
   
   This program contains the files `neutron_flux.f90`, `read_write.f90`, `types.f90`, `main.f90`, `quadrature.f90`, as well as `plots.ipynb`, this `readme.md`, and the `makefile`. 

`neutron_flux.f90` computes the flux values for the first and second parts of the program (box reactor vs box reactor with spherical hollow shell)
`quadrature.f90` uses numerical integration techniques of quadrature to evaluate integrations to compute flux values. Contains Boole's Rule and Monte Carlo quadrature.
`read_write.f90` takes inputs from user and writes data to results files. also initializes calculations and steps through detector distance and radius.
`types.f90` contains types of inputs to be used (real double precision, pi ~ 3.14159..., etc)
`main.f90` calls to read/write module to execute the program. 
`plots.ipynb` can be opened with jupyter notebook, contains graphs of flux over distances and radius measures. 
`makefile` allows user to type "make" into terminal to compile the program.
`readme.md` this file. Contains instructions for the use of this program. 
   `
