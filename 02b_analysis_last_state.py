# -*- coding: utf-8 -*-
"""\
USE: python <PROGNAME> (options) 
OPTIONS:
    -h : print this help message and exit
    -s SEGMENT : Oviduct section to compute [utj, isthmus, ampulla, ia-junction, infundibulum]
    d : max distance to compute [1-10]
"""

import numpy as np
#import pandas as pd
#import datashader as ds
from holoviews import Path as hvPath
import time
#hv.extension('bokeh')
#from holoviews import streams
#import datashader.transfer_functions as tf
#from holoviews.operation.datashader import datashade
#from random import random
import pickle
import bz2
import time
from shapely.geometry import Polygon, LineString
from ev_model import persistency, environment, elements
from ev_model.utilities import geometry as evmetry
from pathlib import Path
import re
import getopt
import sys

#%%
def load_system_state(directory, file_name):
    itno, _, secretory, ciliary, evs, *grids = persistency.read_xml(directory, file_name)
    parameters = ['id', 'x', 'y', 'radius_um', 'age']
    evs_df = persistency.agent_data_to_data_frame(evs, parameters)
    return evs, evs_df


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

#%% Plotly methods
def holes_in_shape_at_d_1(linestring):
    """
    Computes the holes in this shape at a distance of 1 from the edges.
    """
    holes = []
    print(f"{'+' if linestring.is_ring else '-'}ring, {'+' if linestring.is_simple else '-'}simple, {'+' if linestring.is_valid else '-'}valid")
    hole = linestring.parallel_offset(1,'right', resolution=0, mitre_limit=1.0)
    if type(hole) == LineString:
        if len(hole.coords) > 3:
            print(f"{len(hole.coords)} coords, {'+' if hole.is_ring else '-'}ring, {'+' if hole.is_simple else '-'}simple, {'+' if hole.is_valid else '-'}valid")
            holes.append(hole)
    else: # multi holes
        #print('MultiLineString?',type(hole))
        if len(hole) > 0:
            for innerhole in hole:
                if len(innerhole.coords) > 3:
                    print(f"{len(innerhole.coords)} coords, {'+' if innerhole.is_ring else '-'}ring, {'+' if innerhole.is_simple else '-'}simple, {'+' if innerhole.is_valid else '-'}valid")
                    holes.append(innerhole)
    return holes

def produce_interiors_from_exteriors(exteriors):
    interiors = [list() for i in range(len(exteriors))]
    
    inner_holes = 0
    for idx in range(len(exteriors)):
        interiors[idx] = produce_list_of_line_strings(holes_in_shape_at_d_1(
            exteriors[idx]))
        print(f'LineString {idx + 1:3d}/{len(exteriors)} -> {len(interiors[idx])} shapes')
        inner_holes += len(interiors[idx])
    
    print('Computed', inner_holes,'interior shapes for this distance')
    return interiors

def persist_shapes_to_disk(d, interiors, exteriors, interior_fname, exterior_fname):
    interior_pickle = f'resources/analysis/{interior_fname}{d}.pickle.bz2'
    exterior_pickle = f'resources/analysis/{exterior_fname}{d}.pickle.bz2'
    with bz2.BZ2File(interior_pickle, 'w') as target:
        pickle.dump(interiors, target)
        print('Data saved as:', interior_pickle)
    with bz2.BZ2File(exterior_pickle, 'w') as target:
        pickle.dump(exteriors, target)
        print('Data saved as:', exterior_pickle)

def produce_list_of_line_strings(elements):
    new_line_strings = []
    # flatten the nested lists if needed
    for element in elements:
        t = type(element)
        if t is type(list()):
            for ls in element:
                tls = type(ls)
                if tls is LineString:
                    if len(ls.coords) > 3:
                        new_line_strings.append(ls)
                else: 
                    print('Interior not a LineString but a', tls)
        elif t is LineString:
            new_line_strings.append(element)
        else:
            print('interior is not a list or LineString but a', t)
    return new_line_strings

#%%
def compute_all_holes_in_polygons(section, target_distance):
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
    shapes_at_d = {
    1 : {'exterior':[0,1,2,3],
          'interior':[[0], <- produced from exterior 0
                      [0], <- produced from exterior 1
                      [0, 1] <- produced from exterior 2
                      [0]]}, <- produced from exterior 3
    2: {'exterior':[0, 1, 2, 3, 4], <- the interior shapes from the last distance
                                        become the exterior's of the next distance
        'interior':[[0], [], [0, 1], [0]]} <- thus, these would be produced from 
                                  the new exteriors, not from the previous exteriors
    For exteriors with no interior shapes, an empty list should be stored
    The link between interiors-exteriors is given by the indexes from the containing
    lists. Therefore, the shapes in interior[4] would match to those in exterior[4]
    
    The lists of exteriors and interiors for a given distance will be persisted 
    to disk immediately after their computation is completed. These files would
    be independent from those produced for other distances.
    }
    """
    interior_fname = f'data_v2_{section}_interior_d_'
    exterior_fname = f'data_v2_{section}_exterior_d_'
    
    last_d = identify_last_distance_available(interior_fname, target_distance)
    
    shapes_at_d = {} # key: polygon ID, values: dictionary of {distance: holes}

    if last_d > 0:
        # load previous results and store in the general collection
        for d in range(1, last_d + 1):
            shapes_at_d[d] = {'exterior':list(), 'interior':list()}
            
            interior_pickle = 'resources/analysis/' + interior_fname + str(d) + '.pickle.bz2'
            exterior_pickle = 'resources/analysis/' + exterior_fname + str(d) + '.pickle.bz2'
            print('loading data for d=',d,'from', interior_pickle, 'and', exterior_pickle)
            with bz2.BZ2File(exterior_pickle, 'r') as source:
                shapes_at_d[d]['exterior'] = pickle.load(source)
            with bz2.BZ2File(interior_pickle, 'r') as source:
                shapes_at_d[d]['interior'] = pickle.load(source)
            # verify the data loaded is a list of LineString
            for interior in shapes_at_d[d]['interior']:
                if type(interior) == type(list()):
                    for ls in interior:
                        if type(ls) != LineString:
                            raise ValueError('Interior not a LineString but a', type(ls))
                else:
                    raise ValueError('interior not a list but a', type(interior))
            print('  Loaded shapes for distance: ', d)
    
    # NOW we CAN compute the holes for the distances between last_d and target_distance
    for d in range(last_d + 1, target_distance + 1):
        if last_d == 0:
            exteriors = [polygon.exterior for polygon in environment.load_polygons(section)]
            print('This section contains', len(exteriors), 'exterior LineStrings.\nStart processing from d=0...')
        else:
            # obtain the new exteriors from the interiors from the previous distance
            exteriors = produce_list_of_line_strings(shapes_at_d[d - 1]['interior'])
            print('Processing distance', d, 'converting', len(exteriors),'LineStrings (d-1 interior shapes) into exterior shapes for this distance')
        start = time.time()
        interiors = produce_interiors_from_exteriors(exteriors)

        persist_shapes_to_disk(d, interiors, exteriors, interior_fname, exterior_fname)
        elapsed = time.time() - start
        if elapsed < 60:
            print(f'Elapsed {elapsed} seconds')
        else:
            m = elapsed // 60
            print(f'Elapsed {m} minutes and {elapsed - (m * 60)}seconds')

        # add to the general collection
        shapes_at_d[d] = {}
        shapes_at_d[d]['exterior'] = exteriors
        shapes_at_d[d]['interior'] = interiors
        
    return shapes_at_d

#%%
def identify_last_distance_available(fname, target_distance, dir_name='resources/analysis/'):
    target_dir = Path(dir_name)
    if target_dir.exists() and target_dir.is_dir():
        print('target exists and is a directory')
        pattern = '^'+fname+'([0-9]*).pickle.bz2$'
        max_d = 0
        for f in target_dir.iterdir():
            res = re.search(pattern, f.name)
            if res:
                d = int(res.group(1))
                if max_d < d:
                    max_d = d
                    fname = f.name
        if max_d > 0:
            print('This section has been processed before, largest value is', max_d)
            return max_d
        else:
            print('No previous results available')
            return 0
    else:
        raise ValueError('target DOES NOT exist or IS NOT a directory')
        
#%%
#section = 'ia-junction'
#shapes_at_d = compute_all_holes_in_polygons(section, 3)
#identify_last_distance_available(section, 6, './resources/analysis')

#produce_polygon_with_holes()

if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], 'hs:')
    opts = dict(opts)
    distance = args
    
    if '-h' in opts:
        progname = sys.argv[0]
        progname = progname.split('/')[-1]
        help = __doc__.replace('<PROGNAME>', progname, 1)
        print(help, file=sys.stderr)
        sys.exit()
    
    if '-s' in opts:
        if len(args) == 1:
            d = args[0]
            section = opts['-s']
            print('Processing section', section, 'distance',d)
            compute_all_holes_in_polygons(section, int(d))
        else:
            print('ERROR No distance specified')
            print(help, file=sys.stderr)
            sys.exit