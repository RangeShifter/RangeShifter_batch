# RangeShifter Batch Mode <img src="doc/RS_logo.png" align="right" height = 100/>
C++ code for the RangeShifter v2 batch mode application

[RangeShifter](https://rangeshifter.github.io/)
is an eco-evolutionary modelling platform that is becoming 
increasingly used worldwide for both theoretical and applied purposes
[(Bocedi et al. 2014)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12162).

RangeShifter is a spatially-explicit, individual-based simulation platform that 
allows modelling species’ population and range dynamics, such as expansion and shifting, and 
patch connectivity by linking complex local population dynamics and dispersal 
behaviour, while also taking into account inter-individual variability and 
evolutionary processes. RangeShifter is highly flexible in terms of the spatial 
resolution and extent, and regarding the complexity of the considered ecological 
processes. Due to its modular structure, the level of detail in genetics, demographic and 
dispersal processes can be easily adapted to different research questions and 
available data.

This repo contains the source code for the Batch Mode interface of RangeShifter.
In Batch Mode, RangeShifter can be run from the command line (e.g., `./rangeshifter.exe`) within a project directory containing a set of input files.
This allows the user to run large batches of simulations with different parameters, which would need to be specified individually in the GUI version.
The Batch Mode also enables running RangeShifter on machines with a non-interactive interface, for example a high-performance cluster.

## Usage

RangeShifter can be built in Windows or Linux. To build the programm you need to set the operating system in the file Version.h.

For Windows: set #define LINUX_CLUSTER 0
For Linux: set #define LINUX_CLUSTER 1

Then, one can build RangeShifter e.g. with the GNU compiler:

```bash
g++ -o rs.exe ./src/*.cpp ./src/RScore/*.cpp
```

For instructions on how to setup the project directory and input files, please refer to section 3.3 of the [User Manual](https://raw.githubusercontent.com/RangeShifter/RangeShifter-software-and-documentation/master/RangeShifter_v2.0_UserManual.pdf), and to the [documentation repository](https://github.com/RangeShifter/RangeShifter-software-and-documentation) for examples.

## Contributing

Have you spotted a typo, discovered a bug you could fix or would like to propose a new feature? 
We welcome contributions! Please refer to our [contributing guidelines](https://github.com/RangeShifter/RangeShifter_batch_dev/contributing.md) for how to proceed.

## See also

- Compiled software and documentation
- [RScore](https://github.com/RangeShifter/RScore), source for RangeShifter's core code
- [RangeShiftR-pkg](https://github.com/RangeShifter/RangeShiftR-pkg), source for the R interface

## References

 - Bocedi G, Palmer SCF, Pe’er G, Heikkinen RK, Matsinos YG, Watts K, Travis JMJ (2014). 
 *RangeShifter: A Platform for Modelling Spatial Eco-Evolutionary Dynamics and 
 Species’ Responses to Environmental Changes.* Methods in Ecology and Evolution 5: 388–96. 
 - Bocedi G, Palmer SCF, Malchow AK, Zurell D, Watts K, Travis JMJ (2021) *RangeShifter 2.0: An extended and enhanced platform for modelling spatial eco-evolutionary dynamics and species’ responses to environmental changes.* Ecography 44:1453-1462.
