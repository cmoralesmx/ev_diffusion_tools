import numpy as np
import pandas as pd
import datashader as ds
import datashader.transfer_functions as tf

import bz2
import copy
import pickle
import time

from lxml import etree
from pathlib import Path
from shapely.geometry import MultiPolygon, Polygon, Point, LinearRing
from shapely.geometry.polygon import orient
from shapely.prepared import prep

import analysis_last_state as als
from ev_model import persistency, environment, elements
from ev_model.utilities import geometry as evmetry

targets = {}
# need to add the 
targets['ampulla1'] = {'buffer':334}
targets['ampulla2'] = {'buffer':334}
targets['infundibulum'] = {'buffer':77}
targets['ia-junction'] = {'buffer':80}
targets['utj'] = {'buffer':34}
targets['isthmus'] = {'buffer':32}
targets['edges'] = [v for v in np.arange(0.04, 0.17, 0.01)]
targets['sizes'] = [v for v in np.arange(0.04, 0.17, 0.01)]

#sections = ['utj', 'isthmus', 'ia-junction', 'ampulla1', 'ampulla2']
sections = ['']

# deffinitions for v2 of the analysis
classes2 = {'ordering':[('narrow_end','#30a2da'), ('wide_end','#fc4f30'), ('narrow_lumen','#e5ae38'),('wide_lumen','green')]}
classes2['utj'] = {'narrow_end':[1,2,3, 7,8, 13,21], 'wide_end':[4,5,6, 9,14,19], 'narrow_lumen':[15, 16, 17,18], 'wide_lumen':[10,11,12, 20]}
classes2['isthmus'] = {'narrow_end':[1,2,3,4,5,6], 'wide_end':[7,8, 9], 'narrow_lumen':[10, 11, 12, 13, 14, 15], 'wide_lumen':[16,17,18,19,20,21]}
classes2['ia-junction'] = {'narrow_end':[1,2,3,4,5,6,7,8], 'wide_end':[9,10,11,12,13,14,15,16], 'narrow_lumen':[17,18,19,20,21,22,23,24], 'wide_lumen':[25,26,27,28,29,30,31,32]}
classes2['ampulla1'] = {'narrow_end':[1,2,3,4,5,6], 'wide_end':[7,8,9,10, 11, 12], 'narrow_lumen':[13,14,15,16,17,18], 'wide_lumen':[19, 20, 21, 22, 23, 24]}
classes2['ampulla2'] = {'narrow_end':[1,2,3,4,5,6], 'wide_end':[7,8,9,10, 11, 12], 'narrow_lumen':[13,14,15,16,17,18], 'wide_lumen':[19, 20, 21, 22, 23, 24]}

# select which verion of the analysis to use
classes = copy.deepcopy(classes2)

rois_available = als.display_rois_available()
print('--'*10)
experiments_available = als.display_experiments_available()

# v2 compatible
def display_rois_without_selection(user_rois, hvPolys_selected, section_name, rois_per_category, ordering, export_svg=False, png_size=600, legend_position='top_right', labels=True):
    #upolys = hv.Polygons(user_polys if user_polys else [])
    coords_and_category = []
    i = 0
    for roi_id in range(len(user_rois)):
        for class_id in range(len(ordering)):
            #print('class_id:',class_id)
            if roi_id + 1 in rois_per_category[ordering[class_id][0]]:
                # store the coords, category name
                #print(user_rois[roi_id])
                coords_and_category.append( (user_rois[roi_id], ordering[class_id][0], i) )
                #print('id:', i, 'category:', ordering[class_id][0])
                i += 1

    upolys = hv.Polygons([{('x', 'y'):coords, 'category': cl} for (coords, cl, idx) in coords_and_category], vdims='category')
    
    if labels:
        # http://holoviews.org/getting_started/Gridded_Datasets.html
        centroids = [Polygon([pair for pair in zip(p['x'],p['y'])]).centroid.coords[:] for p in upolys.data]

        # https://holoviews.org/reference/elements/bokeh/Labels.html
        labels = hv.Labels([(centroid[0][0],centroid[0][1], f'{coords_and_category[i][2]+1}') for i, centroid in enumerate(centroids)])

    if len(hvPolys_selected) > 1:
        hvPolys_selected_list = hv.Overlay([p for p in reversed(hvPolys_selected)])
    else:
        hvPolys_selected_list = hv.Polygons(hvPolys_selected, label='cross_section').opts(color='darkgrey')

    # produce a HoloViews layout with all the elements
    if labels:
        all_elements = hvPolys_selected_list * upolys.opts(fill_alpha=0.5, cmap=[colour for category, colour in ordering], color='category') * labels
    else:
        all_elements = hvPolys_selected_list * upolys.opts(fill_alpha=0.5, cmap=[colour for category, colour in ordering], color='category')
    # save the plot as svg and png
    hv.save(all_elements.opts(opts.Polygons(show_legend=True, legend_position=legend_position)), f'./resources/analysis/output/{section_name}_cross_section_and_ROIs.png', fmt='png', size=png_size)
    
    if export_svg:
        render =  hv.render(all_elements, backend='bokeh')
        render.output_backend = "svg"
        export_svgs(render, filename=f'./resources/analysis/output/{section_name}_cross_section_and_ROIs.svg')
        
    return all_elements


section = 'ampulla1'
distances_selected = [160,320] #[0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320]
replicates, target_iteration, base_path = 6, 2340000, 'D:\\iterations\\v33\\amp20diestrous\\' #'/home/cmoralesmx/ev_iters/v32/amp/apop_2h01,sec_09,dt01'
#replicates, target_iteration, base_path = 6, 5760000, '/home/cmoralesmx/ev_iters/v32/amp/apop-2h_sec-09_dt-01'
#replicates, target_iteration, base_path = 5, 540000, 'D:\\iterations\\v32\\amp\\apop_2h,sec_09,dt01' #'/home/cmoralesmx/ev_iters/v32/amp/apop_2h01,sec_09,dt01'

rois_to_load = 0
experiment_to_load = None #6

analysis_setup = als.prepare_analysis(section, targets, distances_selected, base_path, replicates, target_iteration)
user_rois_loaded, user_rois, analysis_setup['prep_polys'] = als.load_rois_from_file(rois_available[rois_to_load])

# EVs
analysis_setup = als.load_evs(analysis_setup)
analysis_setup['evs_per_replicate'] = [len(li) for li in analysis_setup['evs_d_useful']]
print(f'Evs per replicate (in full cross section) mean:{np.mean(analysis_setup["evs_per_replicate"])}, sd:{np.std(analysis_setup["evs_per_replicate"]):.2f}')

base_polygon = als.produce_base_polygons(environment.load_polygons('ampulla' if section in ['ampulla1', 'ampulla2'] else section))
#base_hvpolygon = als.produce_hvPolygon(base_polygon)
#all_elements = display_rois_without_selection(user_rois if 'user_rois' in locals() else None, [base_hvpolygon], section, classes[section], classes['ordering'], legend_position='top_right')


# WITHOUT distances
# if the counts do not exist, compute them and save to file
[analysis_setup['evs_in_roi_replicate_objects'], analysis_setup['evs_in_roi_size_replicate_counts'], 
analysis_setup['evs_in_roi_size_replicate']] = als.identify_evs_per_roi(
    analysis_setup['prep_polys'], analysis_setup['evs_d_useful'], targets['sizes'], 
    distance_polygons=None) # distance polygons here are ordered (N-0) yet the result distance indexes must be (0-N)

# export the identified evs per ROI for later reuse
def pickle_data_to_compressed_file(data, name):
    from os import path
    p = path.join('resources','analysis','output',f'{name}.pickle.bz2')
    with bz2.BZ2File(p, 'wb') as compressed_output_file:
        pickle.dump(data, compressed_output_file)
        print(f'data in {name} saved to {p}')

file_names = ['evs_in_roi_replicate_objects', 'evs_in_roi_size_replicate_counts',
    'evs_in_roi_size_replicate']
for fn in file_names:
    pickle_data_to_compressed_file(analysis_setup[fn], fn)