# RangeShifter Batch Mode <img src="doc/RS_logo.png" align="right" height = 100/>
C++ code for the RangeShifter v2 batch mode application

<img title="" src="https://github.com/RangeShifter/RangeShifter_batch_dev/blob/main/doc/rs_batch_logo.png" alt="" align="right" height="150">

[RangeShifter](https://rangeshifter.github.io/) is an eco-evolutionary modelling platform that is becoming 
increasingly used worldwide for both theoretical and applied purposes [(Bocedi et al. 2014)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12162).

RangeShifter is a spatially-explicit, individual-based simulation platform that 
allows modelling species’ population and range dynamics, such as expansion and shifting, and patch connectivity by linking complex local population dynamics and dispersal behaviour, while also taking into account inter-individual variability and 
evolutionary processes. RangeShifter is highly flexible in terms of the spatial 
resolution and extent, and regarding the complexity of the considered ecological 
processes. Due to its modular structure, the level of detail in genetics, demographic and dispersal processes can be easily adapted to different research questions and 
available data.

This repo contains the source code for the Batch Mode interface of RangeShifter.
In Batch Mode, RangeShifter can be run from the command line (e.g., `./rangeshifter.exe`) within a project directory containing a set of input files.
This allows the user to run large batches of simulations with different parameters, which would need to be specified individually in the GUI version.
The Batch Mode also enables running RangeShifter on machines with a non-interactive interface, for example a high-performance cluster.

## Building RangeShifter

The compiled software can be found in the [Software and Documentation](https://github.com/RangeShifter/RangeShifter-software-and-documentation) repo. 

Building RangeShifter from the source code requires CMake. If you haven't done so yet, you will need to [download and install it](https://cmake.org/download/).

RangeShifter can then be configured and built (out-of-source) from `CMakeLists.txt`, with the usual CMake commands:

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

If you use Visual Studio as your IDE, CMake should be recognised automatically when `RangeShifter_batch_dev` is opened as a new folder. 
Visual Studio will take care of the configuration, and you only need to select target RangeShifter.exe before pressing the build button.

Alternatively, RangeShifter can also be built directly with the GNU C++ compiler, in which case some #define macros must be passed to it:

```bash
g++ -o RangeShifter.exe ./src/*.cpp ./src/RScore/*.cpp -DRSDEBUG -DRSWIN64 -DLINUX_CLUSTER
```

## Running RangeShifter

For instructions on how to setup the project directory and input files, please refer to section 3.3 of the [User Manual](https://raw.githubusercontent.com/RangeShifter/RangeShifter-software-and-documentation/master/RangeShifter_v2.0_UserManual.pdf), and to the [documentation repository](https://github.com/RangeShifter/RangeShifter-software-and-documentation) for examples.

## Contributing

See [CONTRIBUTING](https://github.com/RangeShifter/RangeShifter_batch_dev/blob/main/CONTRIBUTING.md)

## See also

- [Compiled software and documentation](https://github.com/RangeShifter/RangeShifter-software-and-documentation)
- [RScore](https://github.com/RangeShifter/RScore), source for RangeShifter's core code
- [RangeShiftR-pkg](https://github.com/RangeShifter/RangeShiftR-pkg), source for the R interface

## Maintainer

- [@TheoPannetier](https://github.com/TheoPannetier)

## References

 - Bocedi G, Palmer SCF, Pe’er G, Heikkinen RK, Matsinos YG, Watts K, Travis JMJ (2014). 
 *RangeShifter: A Platform for Modelling Spatial Eco-Evolutionary Dynamics and 
 Species’ Responses to Environmental Changes.* Methods in Ecology and Evolution 5: 388–96. 
 - Bocedi G, Palmer SCF, Malchow AK, Zurell D, Watts K, Travis JMJ (2021) *RangeShifter 2.0: An extended and enhanced platform for modelling spatial eco-evolutionary dynamics and species’ responses to environmental changes.* Ecography 44:1453-1462.
