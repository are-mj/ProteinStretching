# ProteinStreching
Matlab files for 
- Identifying unfolding/refolding events using optical tweezers
- Calculating probability density vs. transition force histogram arrays
- Fitting model parameters to the Dudko and Bell models

The repository contains files for reading experiment output files created by 
the miniTweezers instrument from TweezersLab (http://tweezerslab.unipr.it)
and extracting unfolding and refolding events when stretching individual protein molecules.
It also contains tools for processing the data to probability densities and model parameters

Matlab Toolboxes used:

- Signal Processing Toolbox
Used by analyse_file and other functions using the experiment text files

- Optimization and Statistics and Machine Learning Toolboxes
Used by fit_unfold_parameters and fit_refold_parameters

All Matlab files have been tested on version R2023a, but may well work with several older versions.