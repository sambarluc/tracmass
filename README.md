Tracmass - MITgcm
=================

This is a fork of the TRACMASS Lagrangian trajectory code.

The aim of this fork is to adapt TRACMASS to work with MITgcm binary output, and its shaved cells in the vertical. The code has also been stripped of all calendar functionality, since it easily lead to conflicts with MITgcm own calendar (the code obviously still keeps track of time).

The fork has been tested on cartesian rectangular and curvilinear MITgcm grids, following several thousands particles for several weeks. No particle hit the boundary, provided a sufficiently high number of intermediate iterations between GCM outputs was used.

The time analytic scheme was tested, modifying the code in a few places to deal with the full 3D flow, but the code tends to produce underflow errors. I did not investigate the cause since the standard iteration method gives good results. I suspect the issue is a trivial one (which does not mean that the solution is trivial).

The fork also contains a python module which reads the TRACMASS output and presents it in an easy to access format based on xarray data objects. Ideally, in the future I will write some pyhton interface between MITgcm and TRACMASS, in particular to set the simulation time in a more streamlined way.

Comments, suggestions, bug reports, etc. are welcome. To do so please open an issue.
