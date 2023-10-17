# Contributing to RangeShifter

## Repo structure
RangeShifter is distributed with three user interfaces, each living in their own repo:
1. the RangeShifter GUI (clickable Windows interface)*
2. [RangeShifter Batch Mode](https://github.com/RangeShifter/RangeShifter_batch) (command line interface)
3. the [RangeShiftR package](https://github.com/RangeShifter/RangeShiftR-pkg) (R interface)

All three share the same source code for the core simulation (i.e., the actual model), which lives in its own repo ([RScore](https://github.com/RangeShifter/RScore)).
Each of the interfaces keeps a copy of RScore, and is kept in sync with it via a git subtree in its RScore/ folder.

⚠️ If you wish to propose a change to the core code of the simulation, please do so *in the [RScore](https://github.com/RangeShifter/RScore) repo*, rather than in the RScore folder in one of the interfaces.

*The RangeShifter GUI is currently being rewritten, and is not open source yet.
