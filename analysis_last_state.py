# -*- coding: utf-8 -*-
"""\
USE: python <PROGNAME> OPTIONS
examples:
python <PROGNAME> -s isthmus -S 0.2 10
python <PROGNAME> -s ampulla -S 0.1 15 -P resources/analysis_simplified

Computes 1um cannals from the boundaries of the cross section towards the
centroid of the cross sections from the five oviduct sections available.

Valid segments: [isthmus, utj, ia-junction, infundibulum, ampulla]
    complexity:  lowest                                    highest
OPTIONS:
    -h : print this help message and exit
    -s SEGMENT : Compute internal shapes at distance from the external shape in the segment
    -p SEGMENT : Produce prepared polygons for topological queries.
    -S TOLERANCE: Enables shape simplification using the tolerance specified [0.0 - 1.0]
    -P PATH: Alternative path where resources are expected
    d : distance to compute [1-10]
"""
import bz2
import copy
import getopt
import holoviews as hv
import numpy as np
import pandas as pd
import pickle
import re
import sys
import time

from bokeh.io import export_svgs
from datetime import datetime
from holoviews import opts
from holoviews import Path as hvPath
from holoviews import Polygons as hvPolygons
from holoviews import streams
from pathlib import Path
from shapely.geometry import MultiPolygon, Point, Polygon, LineString, LinearRing
from shapely.geometry.polygon import orient
from shapely.prepared import prep
#hv.extension('bokeh')
#from holoviews import streams
#import datashader.transfer_functions as tf
#from holoviews.operation.datashader import datashade
#from random import random
from ev_model import persistency, environment, elements
from ev_model.utilities import geometry as evmetry

#%%
def load_evs_from_xml(directory, file_name):
    """
    Produces a dictionary containing:
    'dictionaries': a list of dictionaries represeinting the EVs
    'objects' a list of objects representing the EVs
    """
    evs = persistency.read_xml(directory, file_name, ev_only=True)

    return evs

def load_experiment_data(base_path, replicates=None, target_iteration='2880000'):
    """
    """
    if replicates:
        evs = []
        print('Processing', replicates, 'replicates')
        for i in range(1, replicates+1):
            evs_d = load_experiment_replicate(base_path + '_r' + str(i) + '/', target_iteration)
            evs.append(evs_d)
    else:
        print('Processing a single repeat')
        if base_path[-1:] != '/':
            base_path += '/'
        evs_d = load_experiment_replicate(base_path, target_iteration)
    print('Done loading data')
    return evs

def load_experiment_replicate(base_path, target_iteration, debug=False):
    if type(target_iteration) is not str:
        target_iteration = str(target_iteration)
    # target pickle
    tp = base_path + target_iteration + '.xml' + '.pickle'
    tp_file = Path(tp)

    if tp_file.exists() and not tp_file.is_dir():
        if debug:
            print('Pickled dictionary files exist at')
            print(tp)
        
        # load the pickled dictionaries
        with bz2.BZ2File(tp, 'rb') as pickled_file:
            evs = pickle.load(pickled_file)
            if debug:
                print('Dictionary loaded from:', tp)
        
    else:
        # compute and save
        print('No pickled dictionary files exist yet. Computing now, wait...')
        evs = persistency.read_xml(base_path, target_iteration + '.xml', ev_only=True)
        
        #pickle_file = bz2.BZ2File(tp, 'wb')
        with bz2.BZ2File(tp, 'wb') as pickle_file:
            pickle.dump(evs, pickle_file)
            print('Dictionary saved to:',tp)
    if debug:
        print('Total evs loaded as dictionaries:',len(evs))
    return evs

def filter_oob_from_repeats(shape, all_evs_in_repeats, debug=False):
    if type(all_evs_in_repeats[0]) is not list:
        all_evs_in_repeats = list(all_evs_in_repeats)
    
    total_i = 0
    total_oob = 0
    r = 0
    # iterate over the replicates
    for evs_in_repeat in all_evs_in_repeats:
        oob = []
        # iterate over he evs
        n_evs = len(evs_in_repeat)
        divisor = n_evs // 4
        if debug:
            print('Processing', n_evs,'EVs in replicate',r)
        for i in range(n_evs):
            if i % divisor == 0 and debug:
                print('Checked',i,'/',n_evs)
            ev_p = Point(evs_in_repeat[i]['x'], evs_in_repeat[i]['y'])
            if not shape.contains(ev_p):
                oob.append(i)
                total_oob += 1
        if debug:
            print('OOB in this replicate',len(oob))
        total_i += i

        # remove the identified EVs in current repeat
        if debug:
            print('previous size of replicate', len(all_evs_in_repeats[r]))
        for i in reversed(oob):
            del(all_evs_in_repeats[r][i])
        if debug:
            print('new size of replicate', len(all_evs_in_repeats[r]))
        r += 1
    
    print('Total EVs processed', total_i, 'oob', total_oob)
    return all_evs_in_repeats

# analysis setup
def prepare_analysis(section, targets, distances_selected, base_path, replicates, target_iteration):
    analysis = dict()
    analysis['shrinked_shplyPolys'] = load_shrinked_polygons(section, targets[section]['buffer'])
    print(len(analysis['shrinked_shplyPolys']))
    analysis['optimized_shplyPolys'] = compute_optimized_polygons(analysis['shrinked_shplyPolys'])
    print(len(analysis['optimized_shplyPolys']))
    analysis['hvPolys'] = produce_hvPolygons(analysis['shrinked_shplyPolys'])
    analysis['distances_selected'] = distances_selected
    analysis['base_path'] = base_path
    analysis['replicates'] = replicates
    analysis['target_iteration'] = str(target_iteration)

    analysis['opt_polys_selected'] = []
    analysis['polys_selected'] = []
    analysis['hvPolys_selected'] = []
    analysis['from_to'] = []

    prev_size = len(analysis['optimized_shplyPolys'])-1
    for i in range(len(distances_selected)-2, -1, -1):
        lo = distances_selected[i]
        hi = prev_size
        print(f'EVs within [{lo}-{hi})um from the epithelial tissue')
        prev_size = distances_selected[i]
        analysis['from_to'].append([lo, hi])
        
        analysis['opt_polys_selected'].append(analysis['optimized_shplyPolys'][distances_selected[i]])

        analysis['polys_selected'].append(analysis['shrinked_shplyPolys'][distances_selected[i]])
        analysis['hvPolys_selected'].append(analysis['hvPolys'][distances_selected[i]])
    return analysis

def load_evs(analysis_setup):
    analysis_setup['evs_d'] = load_experiment_data(analysis_setup['base_path'], analysis_setup['replicates'], target_iteration=analysis_setup['target_iteration'])
    # identify and filter oob
    analysis_setup['evs_d_useful'] = filter_oob_from_repeats(analysis_setup['optimized_shplyPolys'][0], analysis_setup['evs_d'])
    return analysis_setup

def display_polygons_available():
    p = './resources/analysis/'
    # find the files starting as 'user_polygons'
    pt = Path(p)

    polygons_available = []
    for f in pt.iterdir():
        if f.name[:13] == 'user_polygons':
            polygons_available.append(f)

    for pair in enumerate(polygons_available):
        print(pair[0], '-', pair[1].name)
    
    return polygons_available

def select_and_load_polygons_specified(polygons_to_load, polygons_available):
    with open(polygons_available[polygons_to_load], 'rb') as polygons_file:
        user_polys_loaded = pickle.load(polygons_file)
        user_polys = copy.copy(user_polys_loaded)
    print('Polygons loaded from:', polygons_available[polygons_to_load].name)
    
    return user_polys_loaded, user_polys

def display_polygons_loaded(user_polys, hvPolys_selected):
    upolys = hv.Polygons(user_polys if user_polys else [])

    # http://holoviews.org/getting_started/Gridded_Datasets.html
    centroids = [Polygon([pair for pair in zip(p['x'],p['y'])]).centroid.coords[:] for p in upolys.data]

    # https://holoviews.org/reference/elements/bokeh/Labels.html
    labels = hv.Labels([(cent[0][0],cent[0][1],i+1) for i, cent in enumerate(centroids)])

    poly_streams = streams.PolyDraw(source=upolys, drag=True)

    n_selected = len(hvPolys_selected) 
    if n_selected == 2:
        hvPolys_selected_list = hvPolys_selected[1] * hvPolys_selected[0]
    elif n_selected == 3:
        hvPolys_selected_list = hvPolys_selected[2] * hvPolys_selected[1] * hvPolys_selected[0]
    elif n_selected == 4:
        hvPolys_selected_list = hvPolys_selected[3] * hvPolys_selected[2] * hvPolys_selected[1] * hvPolys_selected[0]
    elif n_selected == 5:
        hvPolys_selected_list = hvPolys_selected[4] * hvPolys_selected[3] * hvPolys_selected[2] * hvPolys_selected[1] * hvPolys_selected[0]
    elif n_selected == 6:
        hvPolys_selected_list = hvPolys_selected[5] * hvPolys_selected[4] * hvPolys_selected[3] * hvPolys_selected[2] * hvPolys_selected[1] * hvPolys_selected[0]
    elif n_selected == 7:
        hvPolys_selected_list = hvPolys_selected[6] * hvPolys_selected[5] * hvPolys_selected[4] * hvPolys_selected[3] * hvPolys_selected[2] * hvPolys_selected[1] * hvPolys_selected[0]
    elif n_selected == 8:
        hvPolys_selected_list = hvPolys_selected[7] * hvPolys_selected[6] * hvPolys_selected[5] * hvPolys_selected[4] * hvPolys_selected[3] * hvPolys_selected[2] * hvPolys_selected[1] * hvPolys_selected[0]
    elif n_selected == 9:
        hvPolys_selected_list = hvPolys_selected[8] * hvPolys_selected[7] * hvPolys_selected[6] * hvPolys_selected[5] * hvPolys_selected[4] * hvPolys_selected[3] * hvPolys_selected[2] * hvPolys_selected[1] * hvPolys_selected[0]
    elif n_selected == 10:
        hvPolys_selected_list = hvPolys_selected[9] * hvPolys_selected[8] * hvPolys_selected[7] * hvPolys_selected[6] * hvPolys_selected[5] * hvPolys_selected[4] * hvPolys_selected[3] * hvPolys_selected[2] * hvPolys_selected[1] * hvPolys_selected[0]
        
    return (hvPolys_selected_list * upolys.opts(fill_alpha=0.5, active_tools=['poly_draw']) * labels), poly_streams

# stats and plots
def save_polygons_for_reuse(poly_streams, user_polys_loaded, section):
    # fetch the user-created polygons and prepare their coordinates for future use
    user_polys = []
    for i in range((len(poly_streams.data['xs']))):
        if len(poly_streams.data['xs'][i]) > 3:
            user_polys.append([pair for pair in zip(poly_streams.data['xs'][i], poly_streams.data['ys'][i])])

    # Create prepared versions of the user-provided polygons for querying
    n_polys = len(user_polys)
    if n_polys == 1:
        user_polys = list(user_polys)
    # polys = []
    prep_polys = []
    for coords in user_polys:
        p = Polygon(coords)
        #polys.append(p)
        prep_polys.append(prep(p))

    # export the user-selected polygons if needed
    if user_polys_loaded == user_polys:
        print('The polygons did not change')
    else:
        print('The polygons changed, saving to a new file now...')
        dt = datetime.now()
        user_polys_file = f"./resources/analysis/user_polygons_{section}_{dt.strftime('%Y-%m-%d_%H-%M-%S')}.pickle"
        with open(user_polys_file, 'wb') as polys_file:
            pickle.dump(user_polys, polys_file)
        user_polys_loaded = copy.copy(user_polys)

        print(len(user_polys),'user-created polygons are ready for querying. Their coordinates were saved to:', user_polys_file)
    
    return user_polys, prep_polys, user_polys_loaded

def identify_evs_per_polygon(prep_polys, evs_in_replicates, sizes, distance_polygons=None):
    """
    Identifies which of the evs are lying within the user-provided polygons
    Returns:
    a list of the EVs in each polygon (per replicate)
    a list of frequencies of the ev-size per polygon per replicate
    """
    # then, the prepared polygons should be used to produce the statistical information needed
    n_sizes = len(sizes)
    n_polys = len(prep_polys)
    n_replicates = len(evs_in_replicates)

    if 'distance_polygons' in locals():
        print('received', len(distance_polygons), 'distance polygons and',len(prep_polys),'prepared polygons')

    # list used for displaying the EVs using holoviews
    # has the format: evs_per_replicate_in_polygon[replicate][polygon] = [evs in polygon]
    evs_per_replicate_in_polygon = [[list() for p in range(len(prep_polys))] for i in range(n_replicates)]

    # list used for generating the statistical data required
    evs_in_polygon_per_size = np.zeros([len(prep_polys), n_sizes, n_replicates]) # polygons, sizes, replicates
    # polygons, distance polygons, sizes, replicates
    evs_at_distance_in_polygons_per_size = np.zeros([len(prep_polys), len(distance_polygons), n_sizes, n_replicates])

    # this same loop can produce the counts needed for stats purposes
    print('There are',n_replicates,'replicates,', n_polys, 'polygons')
    for r in range(n_replicates):
        print('Processing replicate',r,':', end='')
        # check the evs per replicate against this polygon
        for p in range(n_polys):
            print(' poly', p, 'checking', len(evs_in_replicates[r]),'evs,', end='')
            for ev in evs_in_replicates[r]:
                point = Point(ev['x'], ev['y'])
                if prep_polys[p].contains(point):
                    evs_per_replicate_in_polygon[r][p].append(copy.deepcopy(ev))

                    # if distance-polygons are provided, check them here
                    # bear in mind, the distance polygons are inverted from farthest to nearest to the boundaries
                    # this is to prevent the EVs from being detected at more than 1 distance
                    if distance_polygons:
                        for dp in range(len(distance_polygons)):
                            if distance_polygons[dp].contains(point):
                                distance_polygon = dp
                                break

                    # identify the size range for the ev
                    for s in range(n_sizes):
                        if ev['radius_um'] < sizes[s]:
                            evs_in_polygon_per_size[p, s, r] += 1

                            if distance_polygons:
                                # evs_in_distance_per_polygon_per_size
                                evs_at_distance_in_polygons_per_size[p, distance_polygon, s, r] += 1
                            break
        print()
    if distance_polygons:
        return evs_per_replicate_in_polygon, evs_in_polygon_per_size, evs_at_distance_in_polygons_per_size
    else:
        return evs_per_replicate_in_polygon, evs_in_polygon_per_size

def compute_basic_stats_per_polygon(evs_in_polygon_per_size):
    counts = {}
    max_freq = 0
    for p in range(len(evs_in_polygon_per_size)):
        for s in range(len(evs_in_polygon_per_size[0])):
            numbers = []
            for r in range(len(evs_in_polygon_per_size[0][0])):    
                numbers.append(evs_in_polygon_per_size[p, s, r])
            if p not in counts:
                counts[p] = dict()

            if max_freq < max(numbers):
                max_freq = max(numbers)
            counts[p][s] = {'counts':numbers, 'mean':np.mean(numbers), 'std':np.std(numbers)}
    return counts, max_freq

def compute_stats_per_distance_per_polygon(evs_at_distance_in_polygon_per_size):
    counts = {}
    max_freq = 0
    for p in range(len(evs_at_distance_in_polygon_per_size)): # user polygons
        for d in range(len(evs_at_distance_in_polygon_per_size[0])):
            for s in range(len(evs_at_distance_in_polygon_per_size[0][0])):
                numbers = []
                for r in range(len(evs_at_distance_in_polygon_per_size[0][0][0])):    
                    numbers.append(evs_at_distance_in_polygon_per_size[p, d, s, r])
                if p not in counts:
                    counts[p] = dict()
                if d not in counts[p]:
                    counts[p][d] = dict()

                if max_freq < max(numbers):
                    max_freq = max(numbers)
                counts[p][d][s] = {'counts':numbers, 'mean':np.mean(numbers), 'std':np.std(numbers)}
    return counts, max_freq

def produce_size_distribution_histograms_per_polygon(counts, sizes, edges, max_freq, columns=2, section_name=None):
    freqs = {}
    errors = {}
    histograms = []
    edges = np.array(edges)
    #print(edges)
    # shapes
    for i in range(len(counts)):
        freqs[i] = []
        errors[i] = []
        means = []
        stds = []
        # sizes
        for s in range(len(counts[i])):
            freqs[i].append(counts[i][s]['mean'])
            errors[i].append((sizes[s], counts[i][s]['mean'], counts[i][s]['std']))
            # for creating the table to display values to the user
            means.append(counts[i][s]['mean'])
            stds.append(counts[i][s]['std'])
        freqs[i] = np.array(freqs[i])

        # http://holoviews.org/user_guide/Tabular_Datasets.html
        size_ranges = [f'[{edges[i]:.2f},{edges[i+1]:.2f})' for i in range(len(edges)-1)]
        
        vtable = hv.Table({'Radius range':size_ranges, 'Mean':means[1:], 'SD':stds[1:]}, ['Radius range', 'Mean', 'SD'])

        # http://holoviews.org/Reference_Manual/holoviews.element.html
        histo = hv.Histogram((edges, freqs[i])).opts(
                title=f"Ev size distribution in ROI #{i+1}", 
                xlabel='EV radius in um', 
                ylabel='Frequency',
                xticks= [(edges[i+1],size_ranges[i]) for i in range(len(edges)-1)],
                xrotation=45,
                ylim=(0, max_freq * 1.02), xlim=(0.04,0.17),
                fill_alpha=0.5
            )
        hist_w_errors = (hv.ErrorBars(errors[i]) * histo)
        hv.save(hist_w_errors, f'./resources/analysis/output/{section_name}_output_roi_{i+1}.png', fmt='png')
        
        # exporting directly from bokeh works but depends on the following dependencies plus prior to launching the jupyter-lab executing in the terminal: export OPENSSL_CONF=/etc/ssl/
        #!conda install -c bokeh selenium -y
        #!conda install selenium pillow -y
        #!npm install -g phantomjs-prebuilt
        #
        render =  hv.render(hist_w_errors, backend='bokeh')
        render.output_backend = "svg"
        export_svgs(render, filename=f'./resources/analysis/output/{section_name}_output_roi_{i+1}.svg')

        histograms.append(hist_w_errors + vtable.opts(
            height=400, title=f'Values producing the histogram for ROI #{i+1}')
            )
    
    return hv.Layout(histograms).cols(columns)

def produce_size_distribution_histograms_per_polygon_at_distance(counts, sizes, edges, max_freq, from_to, columns=2, section_name=None):
    freqs = {}
    errors = {}
    histograms = []
    edges = np.array(edges)
    # user polygons
    for p in range(len(counts)):
        # distance
        for d in range(len(counts[p])):
            inverted_d = (len(counts[p]) - 1) -d
            freqs = []
            errors = []
            means = []
            stds = []
            # sizes
            for s in range(len(counts[p][inverted_d])):
                freqs.append(counts[p][inverted_d][s]['mean'])
                errors.append((sizes[s], counts[p][inverted_d][s]['mean'], counts[p][inverted_d][s]['std']))
                # for creating the table to display values to the user
                means.append(counts[p][inverted_d][s]['mean'])
                stds.append(counts[p][inverted_d][s]['std'])
            freqs = np.array(freqs)

            # http://holoviews.org/user_guide/Tabular_Datasets.html
            size_ranges = [f'[{edges[e]:.2f},{edges[e+1]:.2f})' for e in range(len(edges)-1)]
            
            vtable = hv.Table({'Radius range':size_ranges, 'Mean':means[1:], 'SD':stds[1:]}, ['Radius range', 'Mean', 'SD'])

            # http://holoviews.org/Reference_Manual/holoviews.element.html
            histo = hv.Histogram((edges, freqs)).opts(
                    title=f"Ev size distribution, ROI #{p+1} within {from_to[inverted_d]}um", 
                    xlabel='EV radius in um', 
                    ylabel='Frequency',
                    xticks= [(edges[e+1], size_ranges[e]) for e in range(len(edges)-1)],
                    xrotation=45,
                    ylim=(0, max_freq * 1.02), xlim=(0.04,0.17),
                    fill_alpha=0.5
                )
            
            hist_w_errors = (hv.ErrorBars(errors) * histo)
            hv.save(hist_w_errors, f'./resources/analysis/output/{section_name}_output_roi_{p+1}_within_{from_to[inverted_d]}um.png', fmt='png')
        
            # exporting directly from bokeh works but depends on the following dependencies 
            # plus executing in the terminal: export OPENSSL_CONF=/etc/ssl/ prior to launching the jupyter-lab
            #!conda install -c bokeh selenium -y
            #!conda install selenium pillow -y
            #!npm install -g phantomjs-prebuilt
            #
            render =  hv.render(hist_w_errors, backend='bokeh')
            render.output_backend = "svg"
            export_svgs(render, filename=f'./resources/analysis/output/{section_name}_output_roi_{p+1}_within_{from_to[inverted_d]}um.svg')

            histograms.append(hist_w_errors + vtable.opts(
                height=400, title=f'Histogram values, ROI #{p+1} within {from_to[inverted_d]}um')
                )
    return hv.Layout(histograms).cols(columns)

def load_walls_from_xml(directory, file_name, parameters = ['id', 'x', 'y', 'p1x', 'p1y', 'p2x', 'p2y', 'cell_direction']):
    """
    Produces a list of dictionaries and a dataframe with the EVs in a state file specified
    """
    [_, _, secretory, ciliary, _, _, _, _,
        grid_all_walls, _] = persistency.read_xml(directory, file_name)
    all_walls = secretory + ciliary
    all_walls_df = persistency.agent_data_to_data_frame(all_walls, parameters)
    return all_walls, all_walls_df, grid_all_walls

# v2 functions start here

def produce_base_polygons(polygons_loaded):
    """
    Accepts an array of Shapely polygons.
    Identifies the external polygons and the polygons describing holes in those polygons
    Returns an array of Shapely polygons with corresponding holes (where appropiate)
    """
    n = len(polygons_loaded)
    print('Total shapes',n)
    internals = []
    shply_polys = []
    for e in range(n):
        shply_holes = []
        current_internals = []
        for i in range(n):
            if i != e and polygons_loaded[e].contains(polygons_loaded[i]):
                internals.append(i)
                current_internals.append(i)
                shply_holes.append(polygons_loaded[i].exterior.coords)

        if e not in internals:
            #print('External shape',e,'has',len(hv_holes),'holes', end='')
            if len(current_internals) > 0:
                #print(' shapes',current_internals)
                shply_poly = Polygon(polygons_loaded[e].exterior.coords, holes=shply_holes)
            else:
                #print()
                shply_poly = Polygon(polygons_loaded[e].exterior.coords)
            
            # we store the external polygons to a list
            shply_polys.append(shply_poly)
    
    # create the global polygon
    return MultiPolygon(shply_polys)

def produce_hvPolygon(shply_polygons):
    """
    Accepts a Shapely polygon or MultiPolygons (with or without holes).
    Every polygon is assumed to be an external polygon which may contain holes.
    Returns a single holoviews Polygons composed of the n polygons received
    """
    hv_polys = []
    if type(shply_polygons) is Polygon:
        # produce the holoview Polygons directly
        ext_coords = np.array(shply_polygons.exterior.coords)
        hv_polys.append(hvPolygons([{'x':ext_coords[:,0], 'y':ext_coords[:,1]}]))
        
    elif type(shply_polygons) is MultiPolygon:
        for poly in shply_polygons:
            ext_coords = np.array(poly.exterior.coords)
            ni = len(poly.interiors)
            if ni > 0:
                hv_holes = []
                # gather the coordinates for the holes
                for i in range(ni):
                    hv_holes.append(poly.interiors[i].coords[:])
                hv_polys.append(hvPolygons([{'x':ext_coords[:,0], 'y':ext_coords[:,1], 'holes':[hv_holes]}]))
            else:
                # no holes
                hv_polys.append(hvPolygons([{'x':ext_coords[:,0], 'y':ext_coords[:,1]}]))
    else:
        raise TypeError('The provided element is not supported. Expecting a Polygon or MultiPolygon')
    return hvPolygons(hv_polys)

def produce_hvPolygons(shply_polygons):
    """
    Produces the required Shapely polygons for a list of polygons provided
    If the provided element is single Polygon or MultiPolygon, calls the 
    """
    if type(shply_polygons) is list:
        polygons = []
        for poly in shply_polygons:
            polygons.append(produce_hvPolygon(poly))
        return polygons
    elif type(shply_polygons) is Polygon or type(shply_polygons) is MultiPolygon:
        return produce_hvPolygon(shply_polygons)

def compute_shrinked_polygons(source_poly, limit):
    """
    Shrinked polygons can be used for querying the presence of EVs at specific
    distances from the EV boundaries. However, for querying more than a few EVs,
    it is best to use the optimized versions instead of the regular polygons
    """
    polygons_at_d = []
    
    # the original MultiPolygon is extended by 1 um to capture any EV in initial state
    poly_at_d = source_poly.buffer(1)
    polygons_at_d.append(poly_at_d)

    # compute the shrinked down shapely polygons
    for d in range(1, limit):
        poly_at_d = source_poly.buffer(-d)
        polygons_at_d.append(poly_at_d)
    print('Computed',len(polygons_at_d),'shrinked polygons')
    return polygons_at_d

def compute_optimized_polygons(polygons_at_d):
    """
    Uses the polygons provided to generate representation optimized for 
    querying
    """
    prep_polygons_at_d = []
    for poly_at_d in polygons_at_d:
        prep_polygons_at_d.append(prep(poly_at_d))
    return prep_polygons_at_d

def save_shrinked_polygons(polygons_at_d, section_name):
    p = './resources/analysis/shrinked_polygons_' + section_name + '_' + str(len(polygons_at_d)) + '.pickle'
    with bz2.BZ2File(p, 'wb') as pickled_polys:
        pickle.dump(polygons_at_d, pickled_polys)
        print('Polygons for', section_name,'saved as', p)

def load_shrinked_polygons(section_name, limit):
    p = './resources/analysis/shrinked_polygons_' + section_name + '_' + str(limit) + '.pickle'
    with bz2.BZ2File(p, 'rb') as pickled_polys:
        polygons_at_d = pickle.load(pickled_polys)
        print('Polygons loaded from', p)
    return polygons_at_d

# v2 functions end here

# V1 functions start here
# shapely
def produce_polygon_with_holes(polygon, holes):
    if len(polygon.exterior.coords) > 3:
        if type(holes) == type(list()):
            if len(holes) > 1:
                return Polygon(polygon.exterior.coords[:], [h.coords[:] for h in holes])
            else:
                raise ValueError('Less than 1 hole in the list of holes')
        else:
            raise ValueError('Holes is not a list')
    else:
        raise ValueError('Polygon has less coordinates than the min required (3).')

#%%
def produce_path_original(original):
    """
    Produces a HoloViews Path from a Shapely Polygon without holes
    """
    ocoords = np.array(original.exterior.coords[:])
    original = hvPath({'x':ocoords[:,0], 'y':ocoords[:,1]})
    return original

def produce_path_subtracted(subtracted):
    """
    Produces a HoloViews Path from a Shapely Polygon with holes
    """
    ecoords = np.array(subtracted.exterior.coords[:])

    if len(subtracted.interiors) > 1:
        icoords=[interior.coords[:] for interior in subtracted.interiors]
        hv_path = hvPath(icoords) * hvPath({'x':ecoords[:,0], 'y':ecoords[:,1]})
    elif len(subtracted.interiors) == 1:
        icoords = np.array(subtracted.interiors[0].coords[:])
        hv_path = hvPath({'x':icoords[:,0], 'y':icoords[:,1]}) * hvPath({'x':ecoords[:,0], 'y':ecoords[:,1]})
    else:
        hv_path = hvPath({'x':ecoords[:,0], 'y':ecoords[:,1]})

    return hv_path

#%%
# this function calls a non-existent function
def produce_paths(section, polygons, distances = [1,2,3,4,5,6]):
    data = dict()
    data[section] = {
            'subtracted_polygons': dict(),  # polygons
            'original_paths': None,  # polygons
            'subtracted_paths': dict(),
            'within_distance': dict()
        }

    t_section = 0
    print('   produce_paths: Working on section:',section,'containing',len(polygons),'polygons. Distances:',distances)
    for distance in distances:
        print(f'    produce_paths: Processing distance {distance}')
        start = time.time()
        for polygon in polygons:
            print('      produce_paths: nodes in exterior coords of polygon:', len(polygon.exterior.coords))
            data[section]['original_paths'] = produce_path_original(polygon)
            subtracted_polygon = shape_subtract(polygon, distance, 'right')
            data[section]['subtracted_polygons'][distance] = subtracted_polygon
            data[section]['subtracted_paths'][distance] =  produce_path_subtracted(subtracted_polygon)
        t_sec = (time.time() - start)
        t_section += t_sec
        if t_sec/60 < 1.:
            print(f' Done in {t_sec:.2f} seconds')
        else:
            print(f' Done in {t_sec/60:.2f} minutes (actually {t_sec:2f} seconds)')
    print(f'All distances processed for the section in {t_section/60:.2f} minutes')
    return data

#%% Shapely methods
def holes_in_shape_at_d_1(linestring, simplifying = None):
    """
    Computes the holes in this shape at a distance of 1 from the edges.
    """
    holes = []
    print(f"      {'+' if linestring.is_ring else '-'}ring, {'+' if linestring.is_simple else '-'}simple, {'+' if linestring.is_valid else '-'}valid")
    hole = linestring.parallel_offset(1,'left', resolution=0, mitre_limit=1.0)
    if type(hole) == LineString:
        if simplifying:
            hole = hole.simplify(simplifying)
            if type(hole) is not LineString:
                print('Is not a LineString anymore!')
        if len(hole.coords) > 3:
            print(f"        {len(hole.coords)} coords, {'+' if hole.is_ring else '-'}ring, {'+' if hole.is_simple else '-'}simple, {'+' if hole.is_valid else '-'}valid")
            holes.append(hole)
    else: # multi holes
        #print('MultiLineString?',type(hole))
        if len(hole) > 0:
            for innerhole in hole:
                if simplifying:
                    innerhole = innerhole.simplify(simplifying)
                if len(innerhole.coords) > 3:
                    print(f"          {len(innerhole.coords)} coords, {'+' if innerhole.is_ring else '-'}ring, {'+' if innerhole.is_simple else '-'}simple, {'+' if innerhole.is_valid else '-'}valid")
                    holes.append(innerhole)
    return holes

def produce_interiors_from_exteriors(exteriors, simplifying=None):
    # the polygons here would already be simplified
    interiors = [list() for i in range(len(exteriors))]

    inner_holes = 0
    for idx in range(len(exteriors)):
        interiors[idx] = produce_list_of_LinearRings(holes_in_shape_at_d_1(
            exteriors[idx], simplifying))
        print(f'      {idx + 1:3d}/{len(exteriors)} -> {len(interiors[idx])} shapes')
        inner_holes += len(interiors[idx])

    print('Computed', inner_holes,'interior shapes for this distance')
    return interiors

def persist_shapes_to_disk(d, section_name, exteriors, interiors, base_path='resources/analysis'):
    interior_pickle = f'{base_path}/data_v2_{section_name}_{d}_interiors.pickle.bz2'
    exterior_pickle = f'{base_path}/data_v2_{section_name}_{d}_exteriors.pickle.bz2'
    #print('exteriors:', len(exteriors), 'interiors', len(interiors))
    with bz2.BZ2File(interior_pickle, 'w') as target1:
        pickle.dump(interiors, target1)
        print('Data saved as:', interior_pickle)
    with bz2.BZ2File(exterior_pickle, 'w') as target2:
        pickle.dump(exteriors, target2)
        print('Data saved as:', exterior_pickle)

def load_persisted_data(section_name, d, base_path='resources/analysis'):
    interior_pickle = f'{base_path}/data_v2_{section_name}_{d}_interiors.pickle.bz2'
    exterior_pickle = f'{base_path}/data_v2_{section_name}_{d}_exteriors.pickle.bz2'
    print(f'    loading data for d={d} from {interior_pickle} and {exterior_pickle}')

    with bz2.BZ2File(exterior_pickle, 'r') as source:
        exteriors = pickle.load(source)
    with bz2.BZ2File(interior_pickle, 'r') as source:
        interiors = pickle.load(source)
    # verify the data loaded is a list of LinearRings
    print('Checking elements loaded as exteriors')
    check_list_of_LinearRings(exteriors)
    print('Checking elements loaded as interiors')
    check_list_of_LinearRings(interiors)
    return exteriors, interiors

def check_list_of_LinearRings(elements):
    for element in elements:
        if type(element) == type(list()):
            print('It is a list with', len(element), 'elements')
            for elem in element:
                print(type(elem))
                if type(elem) != LinearRing:
                    raise ValueError('Is not a LinearRing but a', type(elem))
        else:
            print(type(element))
            if type(element) is not LinearRing:
                raise ValueError('Is not a list but a', type(element))

def get_LinearRing_gt3(element, t):
    """
    Creates a LinearRing using the 'element' provided.
    This function will only return the LinearRing created if
    it is composed of more than 3 coordinates
    """
    if t is LineString:
        lr = LinearRing(element)
    elif t is LinearRing:
        lr = element

    if len(lr.coords) > 3:
        if lr.is_ccw:
            return lr
        else:
            return orient(Polygon(lr)).exterior

def produce_list_of_LinearRings(elements):
    new_line_strings = []
    # flatten the nested lists if needed
    for element in elements:
        t = type(element)
        if t is type(list()):
            for ls in element:
                tls = type(ls)
                if tls is LineString or tls is LinearRing:
                    lr = get_LinearRing_gt3(ls, tls)
                    if lr:
                        new_line_strings.append(lr)
                else:
                    print('Interior in List not a LineString or LinearRing but a', tls)
        elif t is LinearRing or t is LineString:
            lr = get_LinearRing_gt3(element, t)
            if lr:
                new_line_strings.append(lr)
        else:
            print('interior is not a list or LineString but a', t)
    return new_line_strings

def load_all_holes_in_polygons(section_name, target_distance, base_path='resources/analysis'):
    last_d = identify_last_distance_available(section_name, base_path)

    shapes_at_d = {} # key: 'section_name', values: dictionaries of {distance: holes}

    if last_d >= target_distance:
        # load last results available or target_distance and store in the general collection
        for d in range(1, target_distance + 1 if last_d > target_distance else last_d + 1):
            shapes_at_d[d] = dict()
            shapes_at_d[d]['exterior'], shapes_at_d[d]['interior'] = load_persisted_data(section_name, d, base_path)
            print('  Loaded shapes for distance: ', d)
        return shapes_at_d
    else:
        print('Required not available, needs to be computed first')
        print('Last distance available', last_d)

    return None


#%%
def compute_all_holes_in_polygons(section_name, target_distance, base_path='resources/analysis', simplifying=None):
    """
    1. Check what distances have been computed in the past
    2. If files for previous distances exist:
        load them to the global collection of shapes
    3. Identify what distances must be computed
    4. Compute the missing distances required persisting to disk at each level

    The sources used for computing the shapes is related to the distance being
    computed. For d=1, the original polygons are used. For larger values of d,
    the interior shapes computed for the preceding distance are used. Note: not
    all polygons will exist at every distance, some will eventually disapear as
    their size will be progresively reduced through the process.

    structures in use
    shapes_at_d = {'utj':{
    1 : {'exterior':[0,1,2,3],
          'interior':[[0], <- produced from exterior 0
                      [0], <- produced from exterior 1
                      [0, 1] <- produced from exterior 2
                      [0]]}, <- produced from exterior 3
    2: {'exterior':[0, 1, 2, 3, 4], <- the interior shapes from the last distance
                                        become the exterior's of the next distance
        'interior':[[0], [], [0, 1], [0]]} <- thus, these would be produced from
    }, 'ampulla: {1:{'exterior'=[], 'interior':[]}} ..
                                  the new exteriors, not from the previous exteriors
    For exteriors with no interior shapes, an empty list should be stored.
    The link between interiors-exteriors is given by the indexes from the containing
    lists. Therefore, the shapes in interior[4] would match to those in exterior[4]

    The lists of exteriors and interiors for a given distance will be persisted
    to disk immediately after their computation is completed. These files would
    be independent from those produced for other distances.
    }
    """
    last_d = identify_last_distance_available(section_name, base_path)

    shapes_at_d = {} # key: 'section_name', values: dictionaries of {distance: holes}

    if last_d > 0:
        # load last results available or target_distance and store in the general collection
        d = target_distance if target_distance < last_d else last_d
        shapes_at_d[d] = dict()
        shapes_at_d[d]['exterior'], shapes_at_d[d]['interior'] = load_persisted_data(section_name, d, base_path)
        print('  Loaded shapes for distance: ', d)

    # NOW we CAN compute the holes for the distances between last_d and target_distance
    r = [c for c in range(last_d + 1, target_distance + 1)]
    for d in r:
        if d == 1:
            if simplifying:
                exteriors = [polygon.simplify(simplifying).exterior for polygon in environment.load_polygons(section_name)]
            else:
                exteriors = [polygon.exterior for polygon in environment.load_polygons(section_name)]

            print('    Start processing from d=0...\n      The initial polygon(s) contain(s)', len(exteriors), 'exterior LineStrings.')
        else:
            # obtain the new exteriors from the interiors from the previous distance
            exteriors = produce_list_of_LinearRings(shapes_at_d[d - 1]['interior'])
            print('    Processing distance', d, 'converting', len(exteriors),'(d-1 interior shapes) into exterior shapes for this distance')
        start = time.time()
        interiors = produce_interiors_from_exteriors(exteriors, simplifying)

        if len(interiors) > 1:
            persist_shapes_to_disk(d, section_name, exteriors, interiors, base_path)

            # add to the general collection
            shapes_at_d[d] = {}
            shapes_at_d[d]['exterior'] = exteriors
            shapes_at_d[d]['interior'] = interiors

            elapsed = time.time() - start
            if elapsed < 60:
                print(f'    Elapsed {elapsed:.2f} seconds')
            else:
                m = elapsed // 60
                print(f"    Elapsed {m} {'minutes' if m > 1 else 'minute'} and {elapsed - (m * 60): .2f} seconds")
        else:
            print('0 interiors produced, aborting persistency to disk')
            print('Aborting the rest of the process')
            break

    return shapes_at_d

#%%
def identify_last_distance_available(section_name, base_path='resources/analysis/'):
    target_dir = Path(base_path)
    if target_dir.exists() and target_dir.is_dir():
        print(f'Target path {base_path} exists and is a directory, checking it\'s content')
        pattern = f'^data_v2_{section_name}_([0-9]*)_interiors.pickle.bz2$'
        max_d = 0
        for f in target_dir.iterdir():
            res = re.search(pattern, f.name)
            if res:
                d = int(res.group(1))
                if max_d < d:
                    max_d = d
        if max_d > 0:
            print('  This section has been processed before, largest value is', max_d)
            return max_d
        else:
            print('  No previous results available')
            return 0
    else:
        raise ValueError('target DOES NOT exist or IS NOT a directory')

def produce_prepared_polygons(section_name, d, base_path='resources/analysis'):
    from shapely.prepared import prep

    last_d = identify_last_distance_available(section_name, base_path)

    if last_d == 0:
        print('ERROR, exterior and interior shapes are not available')
    elif d > last_d:
        print('ERROR, exterior and interior shapes are not available at the desired distance')
    else:
        # load previous results and store in the general collection
        exteriors, interiors = load_persisted_data(section_name, d, base_path)
        print('  Loaded shapes for distance: ', d)

        polygons = []
        prepared_polygons = []

        for eid in range(len(exteriors)):
            poly = Polygon(exteriors[eid], [interior for interior in interiors])

            polygons.append(poly)
            prepared_polygons.append(prep(poly))

        return polygons, prepared_polygons

def show_help(progname):
    progname = progname.split('/')[-1]
    help = __doc__.replace('<PROGNAME>', progname)
    print(help, file=sys.stderr)

if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], 'hS:s:P:')
    opts = dict(opts)
    distance = args
    base_path = 'resources/analysis'
    simplifying = None

    if '-h' in opts:
        show_help(sys.argv[0])
        sys.exit()

    if '-S' in opts:
        simplifying = float(opts['-S'])
        print('Polygon simplification enabled, tolerance:', simplifying)

    if '-P' in opts:
        base_path = opts['P']
        print('Using non-default base path:', base_path)

    if '-s' in opts:
        if len(args) == 1:
            d = args[0]
            section = opts['-s']
            print('Processing section', section, 'distance',d)
            compute_all_holes_in_polygons(section, int(d), base_path, simplifying)
        else:
            print('ERROR No distance specified')
            print(help, file=sys.stderr)
            sys.exit
    elif '-p' in opts:
        if len(args) == 1:
            d = args[0]
            section = opts['-o']
            print('Producing prepared polygons for topological querys for section', section, 'distance', d)
            produce_prepared_polygons(section, distance, base_path)
        else:
            print('ERROR, no distance specified')
            print(help, file=sys.stderr)
            sys.exit
    else:
        show_help(sys.argv[0])
        sys.exit()

