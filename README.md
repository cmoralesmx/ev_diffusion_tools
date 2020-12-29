# EV diffusion tools

The python scripts in this repository depend on the following libraries

- holoviews
- bokeh
- jupyter-lab

It is recommended to create a conda environment where these dependencies are installed.
To do so, the following instructions should be executed. [source](https://holoviews.org/install.html)
Note: the name of the environment is optional, the rest is mandatory

```python
# modern versions of conda
conda create --name phd python=3.6
conda activate phd
conda install -c pyviz holoviews bokeh
conda install -c conda-forge jupyterlab
# nodejs is needed for the following jupyter-lab extension and for producing some EV plots
conda install nodejs
jupyter labextension install @pyviz/jupyterlab_pyviz
conda install shapely xmltodict datashader selenium pillow lxml
# the following is needed if exporting to SVG?
#npm install -g phantomjs-prebuilt
```
The conda environment produced was last checked on the 15 of December, 2020

This repository contain scripts to execute most of the functionality needed
for the EV model. Some are jupyter notebooks, they depend on elements only
available on this platform i.e.: interactive widgets

## Pre simulation
- `00_Xsec-contours.ipynb` contains the image processing code to identify the
contours of the cross-sections of the oviduct. The five key regions of the
oviduct are prepared for processing in single independent notebook cells.
This way, each can be processed independently from the others and there is no
guessing about the pipeline to execute in the future.
Each of these cross sections will produce a set of polygons stored in 'pickle'
format for the next processing step

- `01_setup_ampulla_image.ipynb` contains code for producing the FlameGPU
compatible data from the polygons produced from the previous file. The section
titled "simple tests" contains code to produce simple test scenarios useful for
debugging the FlameGPU model. After this, some detail is provided about the EV
behaviour in our model. This is followed by the section where each of the cross
sections of the oviduct can produce the FlameGPU file format of the corresponding
initial states file (0.xml)


## Pre analysis
- `02_define_ROIs.ipynb` This file is contains the functionality for producing
the Regions of Interest (ROIs) to use for analysing the data in the next step. It
requires the same polygon files produced by the contour detection notebook
(00_Xsec-contours.ipynb`). Unfortunately, this process is highly manual. The
ROIs must be defined by creating rectangles/polygons in the provided areai
(V1, V2, manually defined ROI). The produced ROIs must be saved manually too.
Functionallity is provided for checking the size of the polygons. This check
can only show the size of the active polygon, which was selected/modified.
Similarly, code is provided for saving the data in the correct data format and
for visualising without modifying for ROIs created earlier.

## Analysis
For the data analysis stage, the pipeline changed a lot. It started as a custom
implementation in Python using pandas, evolving into up to four alternate
versions in parallel until eventually only two survived. These are represented
by the notebook which file names are suffixed by v2 and v4. Each corresponding
to a different version of the ROIs.
However, as the work progressed, and the queries became more complex, I had to
stop using these tools and instead rely on TIBCO Statistica and R.
To do so, two Python CLI scripts were created. The first one, `locator.py` does
identify which EVs are within the ROIs and produces the corresponding output.
The second one, `statistics.py` processes those files to generate CSV formated
output. These can be later used in numerical computations packages for doing
the statistical analysis. Both scrpts are self documented, still their
documentation is reproduced further down this file for convenience.

- `03_analysis_last_state_only_v2.ipynb` DEPRECATED
- `03_analysis_last_state_only_v4.ipynb` DEPRECATED
- `03_analysis_v2_batch.py` DEPRECATED
- `03_analysis_v2b_batch.py` DEPRECATED
- `locator.py` Identifies the EVs positioned within the ROIs. Uses multi
programming to reduce the time required for checking the N agents against
the M regions of interest
- `statistics.py` Exports .CSV files from the files produced by `locator.py`
These files can be loaded in Statistica, R or any numerical package for
further computations.

### Demos, tests, Examples

- `0d_analysis_displac.ipynb` An early version of the analysis of EV
displacement. It can be used to check the displacement corresponds to a random
walk.

- `Diffusion_AnalyticalSolution.ipynb` contains an early exploration of this
process from a theoretical perspective
- `analysis_last_state.py` A CLI version of the the deprecated custom
analytical tools
- `brownianMotion.ipynb` An early implementation of the simulation of brownian
motion targeting the widgets in jupyter notebook for visualisation purposes
- `collision_detection.ipynb` Created for debugging the collision detection
- `collision_detection2.ipynb` Created for debugging the collition detection
- `count_disabled.py` Reads an XML save file and counts the EVs in a disabled
state.
- `cross_section_starting_points.py` DEPRECATED An early implementation of the
algorithm computing the valid starting points for a cross section
- `data.npy` was this created by a third-party tool?
- `evDif_v36_continueSim.py` Processes an XML save state file produced by the
FlameGPU-based EV model to allow continuing the simulation taking a given
save state as starting point. This is needed to clear runtime values < 1e^-4
which are recorded as 0.0000 by the FlameGPU file format. These values must be
cleared prior to restarting the simulation otherwise the state of the
simulation becomes erroneous.
- `samples_hv_barplots.ipynb` Plotting samples, holoviews version
- `samples_pandas.ipynb` Plotting samples, pandas version
- `samples_stats.ipynb` Stattistical samples using pandas
