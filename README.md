
Echelle Layout Program 2
===============================================

Marc Schneider, 2021 - 2022
marc.schneider@kit.edu

When designing integrated wavelength division multiplexing systems using
photonic chips, different methods for wavelength filtering are available.
One of them are planar concave gratings or Échelle gratings. For such
Echelle gratings, differend construction methods exist, like Rowland
Circle method and Two Stigmatic Point (TSP) method.

This MATLAB app implements my version of the TSP method to calculate
arbitrary Echelle gratings with Bragg grating reflectors. From this, a
simulation project for COMSOL Multiphysics can be generated to simulate
the Echelle grating in 2D, as well as a GDSII layout for using the
grating in a chip design. Additionally a sample file for Synopsys
OptoDesigner is generated.

The project is published under the MIT license (see LICENSE), with the
following exception:

- The function 'gdsii_boundarytext_Bevel2.m' is a modified version of
  a function from the gdsii-toolbox of Ulf Griesmann
  (https://github.com/ulfgri/gdsii-toolbox). The toolbox is in the
  public domain, so I don't want to put my additions under a more
  restrictive license. Therefore this function is also in the public domain.
- The included logo pictures are not covered by the MIT license.

This software was developed at the Karlsruhe Institute of Technology (KIT),
Germany. This software is an experimental system. KIT assumes no
responsibility whatsoever for its use by other parties, and makes no
guarantees, expressed or implied, about its quality, reliability, or any
other characteristics.



To run the program to its full extend, you still need
- GDSII Toolbox v1.41 from Ulf Griesmann (https://github.com/ulfgri/gdsii-toolbox)
- COMSOL Multiphysics (installed with MATLAB support)

To speed up calculations, it's preferred to install the
- Parallel Computing Toolbox (Matlab)


The file "EffectiveRefractiveIndices_Wavelength1480-1630nm_SiliconThickness205-225nm+245-255nm.xlsx"
contains the effective refractive indices of the layer stack used for the
Echelle grating. The given file contains calculated effective refractive
indices for a layer stack of silicon dioxide - silicon - silicon dioxide
for wavelengths between 1480nm and 1630nm and for silicon thicknesses
between 205nm and 225nm as well as between 245nm and 255nm in 1nm steps.
The data is used for the Material Chooser. Extend the file to your own needs.

The main program is "MarcEchelle.mlapp" and runs in Matlabs App Designer.
Or at least it should, if all required parts are included in Matlabs search path.


Help
----
If you find a bug in the software, please send a message to 
marc.schneider@kit.edu and I will try to fix it.
