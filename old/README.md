# OrbitSim
#### A basic gravity/orbit simulator written in python.

## Notes
There are a few things to note:
- The code is very messy, unrefined, and scattered all over the place
- Simulating more than 5 bodies at once will cause the efficiency to drop rapidly
- In order to animate the simlution, ffmpeg will need to be installed on the machine running the code (animating may also take a while)
- The code is very inneficient, and relies on the most basic form of numerical integration

## Running the Sim
Currently the notebook is a mess, some of the cells have errors, left over from experimentation. The runsim.py file is a better depiction of how a simulation can be run.

Among the libraries that would need to be installed is:
- numpy
- scipy
- matplotlib

## Future
In the future I plan to rewrite the simulator in a faster language (C++ most likely), and implementing different algorithms to improve efficiency in general, as well as 
allow a large amount of bodies to be simulated. I also want to investigate either a graphics library or even game engine to make less painful animations, and perhaps have
the simulation be rendered realtime.


