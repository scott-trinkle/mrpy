'''
I have created complete wrapper functions for the "tckgen" and "tckmap"
commands from MRTrix3.

These have the exact same capability as the command-line MRTrix3 functions,
themselves, except with options presented as keyword arguments to Python
functions. For information on all possible options, see the documentation:

tckgen:
https://mrtrix.readthedocs.io/en/latest/reference/commands/tckgen.html#tckgen

tckmap:
https://mrtrix.readthedocs.io/en/latest/reference/commands/tckmap.html

or, run "tckgen" or "tckmap" in the terminal
'''

# First import the package
import mrpy as mr

# Let's do a whole-brain tractography run from the HARDI dataset.
mr.tckgen(source='odfs.nii.gz',  # input odf file
          tracks='wb_tract.tck',  # output .tck file
          select=50000,  # default is 500000, making smaller for time
          # seeding from the whole brain mask
          seed_rejection='brain_mask.nii.gz',
          force=True)  # overwrite existing files

# And that's it! All other values are set to sensible defaults, but we can
# change anything we want, just use the keyword from the tckgen documentation
# and input whatever values you want. The exception to that is
# angle/curvature. MRTrix3 uses "angle" as a keyword, I prefer to use
# "curvature". The MRPy "tckgen" function will accept both - if you input
# "angle", it will use it, if you input "curvature", it will convert it to an
# angle and use that. If you input both, it will use the input for "angle"

# Now let's create two track density images (TDI) from our .tck file
mr.tckmap(tracks='wb_tract.tck',  # the tracks we just generated
          output='wb_tract.nii.gz',  # the output tdi image
          vox=0.05)  # we want the output spacing to be 50 um (0.05 mm)

# This will give us a voxel image where the count in each voxel is a weighted
# sum of streamlines passing through that voxel. We can also create a separate
# tdi with "directionally encoded color" (DEC), where each voxel is assigned
# an RGB color corresponding to the primary orientation of the streamlines.

# We can do that with the "dec=True" option
mr.tckmap(tracks='wb_tract.tck',
          output='wb_tract_dec.nii.gz',  # note the filename change
          vox=0.05,  # same output spacing
          dec=True)
