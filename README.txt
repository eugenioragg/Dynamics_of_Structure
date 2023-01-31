Sdofvib
-------

Sdofvib is a Python package for single-degree-of-freedom vibration analysis.


### INSTALLATION

Download the zip file sdofvib.zip, extract it and place it in your path.


### CONTENT

Driver    : sdofvib.py
Functions : timestep.py
			lineplot.py


### USAGE

Run the package from the driver file sdofvib.
Fill out the input data under INPUT DATA and run.
The driver then calls timestep(), which returns the response in terms of displacement, velocity and acceleration, as well as the discrete times and forces.

For post-processing, the package contains the plot function lineplot().
This is a general plot function for 2D line plots, eg. time histories.