import bz2
import copy
import pickle
from multiprocessing import Process, Manager, freeze_support

import numpy as np
import pandas as pd

import analysis_last_state as als
from ev_model import environment


def compute_areas_and_evs_per_cross_section(section, analysis_setup,
                                            user_rois):
    """
    Compute the area of the outer most polygon of the specified section.

    Args:
        section: the section to compute the area from, can be one of ['isthmus',
        'ampulla', 'utj', 'ampulla1', 'ampulla2', 'ia-junction']
        analysis_setup: a dictionary produced by X?
        user_rois: the Regions of Interest to intersect and

    Returns:
        [
        the total area of cross-section in um^2,
        the total area of cross-section in mm^2,
        a list of areas per ROI,
        a list of counts of EVs per ROI,
        ]
    """
    from shapely.geometry import Polygon

    intersecting_polys = []
    base_polygon = als.produce_base_polygons(
        environment.load_polygons(
            'ampulla' if section in
            ['ampulla', 'ampulla1', 'ampulla2'] else section))
    area_mm2 = base_polygon.area / 1000000.
    print(f'Area of base polygon: {base_polygon.area:.3f} um^2 '
          f'or {area_mm2} mm^2')
    roi_areas_um2 = []
    evs_per_roi_counts = []
    print('total ROIs:', len(user_rois))

    if 'evs_in_roi_replicate_objects' in analysis_setup:
        for user_rois_idx in range(len(user_rois)):
            # Compute the polygons resulting from the intersection of
            # the ROIs and the base polygon
            user_poly = Polygon(user_rois[user_rois_idx])
            intersected_polygon = base_polygon.intersection(user_poly)

            # Compute the number concentration per replicate in this ROI
            evs_per_replicate = []
            for rep_idx in range(
                    len(analysis_setup['evs_in_roi_replicate_objects'][0])):
                evs_per_replicate.append(
                    len(analysis_setup['evs_in_roi_replicate_objects']
                        [user_rois_idx][rep_idx]))
            # collect the total EVs per replicate in ROI
            evs_per_roi_counts.append(evs_per_replicate)
            # collect the area per ROI
            roi_areas_um2.append(intersected_polygon.area)
            # collect the intersected polygons
            intersecting_polys.append(intersected_polygon)
    else:
        raise RuntimeError(
            'evs_per_replicate_in_polygon does not exist in analysis_setup.'
            'Did you load the experiment counts?')
    return base_polygon.area, area_mm2, roi_areas_um2, evs_per_roi_counts


def not_valid_ev_locator(replicate_id, evs_in_replicate, poly_coords):
    """
    Identify EVs out of bounds of the polygons representing the environment.

    Under normal circumstances, when new EVs are introduced in the system, they
    are created outside the polygon and are displaced to the interior mimmicking
    the biological release from the cell

    Althought, this is fine for the model, they could be considered not valid
    for some processes i.e.: for debugging collision detection

    Args:
        replicate_id: the index of the replicate to check
        evs_in_replicate: a list of dictionaries, one for each EV in this
                          replicate
        poly_coords: A list of coordinates. They must produce a valid Polygon
    Returns:
        void, stores the index of the EV OOB in the multiprocessign queue
    """
    from shapely.geometry import Point
    from shapely.prepared import prep
    not_valid = []
    ppol = prep(poly_coords)
    for i in range(len(evs_in_replicate)):
        if evs_in_replicate[i]['defaultAt'] > 0 and evs_in_replicate[i][
                'disabledAt'] < 1:
            p = Point(evs_in_replicate[i]['x'], evs_in_replicate[i]['y'])

            if not ppol.contains(p):
                not_valid.append(i)
        else:
            not_valid.append(i)
    _q.put((replicate_id, not_valid))


def ev_in_roi_locator(rep_id, roi_id, evs_in_replicate, roi_coords):
    """
    Finds what EVs are inside which ROI per replicate

    Args:
        rep_id: the index of the replicate to work with
        roi_id: the index of the Region of Interest to work with
        evs_in_replicate: a list of dictionaries, one per EV in replicate
        roi_coords: a list of coordinates for each vertice of a polygon
    Returns:
        void, uses the multiprocessign queue to store the values identified
    """
    from shapely.geometry import Polygon, Point
    from shapely.prepared import prep

    pol = Polygon(roi_coords)
    ppol = prep(pol)
    evs_here = []
    for ev in evs_in_replicate:
        p = Point(ev['x'], ev['y'])
        if ppol.contains(p):
            evs_here.append(ev)
    _q.put((roi_id, rep_id, evs_here))


def identify_valid_evs_multiprocess(poly_coords, evs_in_replicates):
    """
    Identify EVs out of bounds. Multiprocessing version

    More details are available on `not_valid_ev_locator`

    Args:
        poly_coords: A list of lists of vertices per polygon
    Returns:
        evs_in_replicates: A list of list of indices of the EVs per replicate
    """
    n_replicates = len(evs_in_replicates)
    print('Will check', n_replicates, 'replicates. \nProcessing, wait...')

    processes = []
    # iterate the list of replicates, each is handled by 1 processor
    for rep_id in range(n_replicates):
        p = Process(target=not_valid_ev_locator,
                    args=(rep_id, evs_in_replicates[rep_id], poly_coords))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()
    total_i = 0
    total_non_valid = 0
    for p in processes:
        res = _q.get()
        total_non_valid += len(res[1])
        nrep = len(evs_in_replicates[res[0]])
        print(f'EV in replicate: {nrep}, non-valid: {len(res[1])} '
              f'({(100 / nrep) * len(res[1]) :.3f}%)')
        for i in reversed(res[1]):
            del (evs_in_replicates[res[0]][i])
        total_i += len(evs_in_replicates[res[0]])

    print(f'Total EVs processed: {total_i}, non-valid: {total_non_valid} '
          f'({(100 / total_i) * total_non_valid:.3f}%)')
    return evs_in_replicates


def identify_evs_per_roi_multiprocess(polys, evs_in_replicates):
    """
    Find the EVs located inside the ROIs. Multiprocessing version

    Args:
        polys: a list of shapely polygons
        evs_in_replicates: a list of lists of dictionaries, one per EV
    Returns:
       evs_in_roi_replicate_objects: a list of lists of EVs per ROI per Replicate
       evs_in_roi_replicate_counts: a list of the total EVs per ROI per Replicate
       evs_in_roi_replicate_radius_age: a list of tuples (ROI id, replicate id,
                                        EV radius, age)
    """
    n_polys = len(polys)
    n_replicates = len(evs_in_replicates)
    evs_in_roi_replicate_objects = [[list() for p in range(n_replicates)]
                                    for i in range(n_polys)]
    print('n_polys:', n_polys, 'n_replicates:', n_replicates)
    evs_in_roi_replicate_counts = np.zeros([n_polys, n_replicates])
    print(evs_in_roi_replicate_counts.shape)
    evs_in_roi_replicate_radius_age = list()

    print('Will check', n_polys, 'rois in', n_replicates,
          'replicates. \nProcessing, wait...')
    for rep_id in range(n_replicates):
        print('Rep', rep_id, end='... ')
        processes = []
        # iterate the list of ROIs, each ROI is handled by 1 processor
        for roi_id in range(n_polys):
            p = Process(target=ev_in_roi_locator,
                        args=(rep_id, roi_id, evs_in_replicates[rep_id],
                              polys[roi_id]))
            processes.append(p)
            p.start()

        for p in processes:
            p.join()

        for p in processes:
            # res [0:ROI id, 1:Replicate id]
            res = _q.get()

            for ev in res[2]:
                evs_in_roi_replicate_objects[res[0]][res[1]].append(ev)

                evs_in_roi_replicate_counts[res[0], res[1]] += 1

                evs_in_roi_replicate_radius_age.append(
                    (res[0], res[1], ev['radius_um'], ev['age']))
        print('COMPLETED')
    # Display a list of Replicate - ROI - # of elements
    # which can be used for an initial assesment of the system or task
    print('Rep ROI Element')
    total_in_rois = 0
    for roi in range(len(evs_in_roi_replicate_objects)):
        for rep in range(len(evs_in_roi_replicate_objects[roi])):
            total_in_rois += len(evs_in_roi_replicate_objects[roi][rep])
            print(f'{rep:2d}   {roi:2d}   '
                  f'{len(evs_in_roi_replicate_objects[roi][rep]):4d}')
    print('Total EVs in ROIs:', total_in_rois)
    return [
        evs_in_roi_replicate_objects, evs_in_roi_replicate_counts,
        evs_in_roi_replicate_radius_age
    ]


def pickle_data_to_compressed_file(data, name, version, iteration):
    """Saves the data as pickle objects"""

    rao = './resources/analysis/output/'
    vni = f'{version}_{name}_{iteration}'
    with bz2.BZ2File(f'{rao}{vni}.pickle.bz2', 'wb') as compressed_output_file:
        pickle.dump(data, compressed_output_file)
        print(f'data in {name} saved to {rao}{vni}.pickle.bz2')


_targets = {}
# We define the maximum distance for each section that will still produce
# valid polygons when shrinking the original polygon to produce inner regions
# These inner polygons can be used to identify the EVs located at specific
# distance from the edges of the oviduct
_targets['infundibulum'] = {'buffer': 77}
_targets['ia-junction'] = {'buffer': 80}
_targets['utj'] = {'buffer': 34}
_targets['isthmus'] = {'buffer': 32}
_targets['ampulla'] = {'buffer': 334}

sections = ['isthmus', 'ampulla']

def main(section, base_path, iteration, replicates, load_rois, min_distance,
         max_distance, version, streaming, force_overwrite):
    """
    Identify what EVs are within the regions of interest (ROI) per replicate.
    This computation is intensive. The counts produced are stored for later use

    Args:
        section: the name of the oviduct section to work with
        base_path: a file path to take as the working directory for this script
        iteration: the iteration number to work with
        replicates: the number of replicates to work with
        load_rois: the file containing the polygons to load
        min_distance: ??
        max_distance: ??
        version: string to identify what version of the experiment to work with
        streaming: ??
        force_overwrite: Flag to enable rewriting the output files if they exist
    """
    distances_selected = [min_distance, max_distance]

    analysis_setup = als.prepare_analysis(section, _targets, distances_selected,
                                          base_path, replicates, iteration)
    [_, user_rois,
     analysis_setup['prep_polys']] = als.load_rois_from_file(load_rois)

    analysis_setup['evs_d'] = als.load_experiment_data(
        analysis_setup['base_path'],
        analysis_setup['replicates'],
        target_iteration=analysis_setup['target_iteration'],
        streaming=streaming,
        force_overwrite=force_overwrite)
    base_polygon = als.produce_base_polygons(
        environment.load_polygons(section))
    analysis_setup['evs_d_useful'] = identify_valid_evs_multiprocess(
        base_polygon, analysis_setup['evs_d'])

    analysis_setup['evs_per_replicate'] = [
        len(li) for li in analysis_setup['evs_d_useful']
    ]
    print(f'Evs per replicate '
          f'mean:{np.mean(analysis_setup["evs_per_replicate"])}, '
          f'sd:{np.std(analysis_setup["evs_per_replicate"]):.2f}')
    # base_polygon

    # MULTIPROCESS
    [
        analysis_setup['evs_in_roi_replicate_objects'],
        analysis_setup['evs_in_roi_replicate_counts'],
        analysis_setup['evs_in_roi_replicate_radius_age']
    ] = identify_evs_per_roi_multiprocess(user_rois,
                                          analysis_setup['evs_d_useful'])

    stats = dict()

    print('Computing stats for', section)
    [
        stats['total_area_um2'], stats['total_area_mm2'], stats['roi_areas'],
        stats['total_evs_per_roi']
    ] = compute_areas_and_evs_per_cross_section(section, analysis_setup,
                                                user_rois)

    # total evs per cross section
    stats['total_evs_per_section'] = analysis_setup['evs_per_replicate']
    stats['evs_in_roi_replicate_counts'] = analysis_setup[
        'evs_in_roi_replicate_counts']
    stats['evs_in_roi_replicate_radius_age'] = analysis_setup[
        'evs_in_roi_replicate_radius_age']

    stats['evs_per_replicate'] = [
        (rep, ev['radius_um'], ev['age'])
        for rep in range(len(analysis_setup['evs_d_useful']))
        for ev in analysis_setup['evs_d_useful'][rep]
    ]

    m_index = pd.MultiIndex.from_tuples([(section, t[0])
                                         for t in stats['evs_per_replicate']],
                                        names=['section', 'replicate'])
    stats['evs_per_replicate_df'] = pd.DataFrame(
        [t[1:] for t in stats['evs_per_replicate']],
        index=m_index,
        columns=['radius_um', 'age'])

    m_index = pd.MultiIndex.from_tuples(
        [(section, ) + t[:-2]
         for t in analysis_setup['evs_in_roi_replicate_radius_age']],
        names=['section', 'roi', 'replicate'])
    stats['data_frame'] = pd.DataFrame(
        [t[-2:] for t in analysis_setup['evs_in_roi_replicate_radius_age']],
        index=m_index,
        columns=['radius_um', 'age'])

    print(f'Cross section total evs:'
          f'{np.mean(stats["total_evs_per_section"]):.2f}')

    # SAVE the stats for this section

    sf = f'./resources/analysis/output/stats_{version}_{iteration}_'\
         f'{load_rois.name[5:14]}.pickle.bz2'
    with bz2.BZ2File(sf, 'wb') as stats_file:
        pickle.dump(stats, stats_file)
        print('stats saved to ', sf)


_m = Manager()
_q = _m.Queue()

if __name__ == '__main__':
    """
    Process XML files produced by the FlameGPU-based EV model to find the agents
    located inside the regions of interest

    Args:
        path: string location in the file system where to look for the input files
        section: string to identify section of the oviduct to work with
        version: string to identify the version of the experiment to work with
        iteration: the iteration number to work with
        replicates: the number of replicates to work with
        polygons: the polygons to load
        force: force writing the output if the files exist
    Returns:
        void: the output is written to screen and files are saved to drive
    """
    import argparse
    freeze_support()

    parser = argparse.ArgumentParser('Analyse EV experimental data')
    parser.add_argument('path',
                        metavar='base_path',
                        help='Directory where the target files are located')
    parser.add_argument(
        'section',
        choices=['isthmus', 'ampulla'],
        help='The oviduct section to analyse. '
        'The name is used to select the ROIs to load/use during processing')
    parser.add_argument(
        'version',
        metavar='version',
        help='Aids identification of the data produced. '
        'The files generated by this script will have this value in their name'
    )
    parser.add_argument('-i',
                        '--iteration',
                        nargs='?',
                        type=int,
                        help='Iteration number to read. Default 14.4k')
    parser.add_argument('-r',
                        '--replicates',
                        nargs='?',
                        type=int,
                        help='The number of replicates to read')
    parser.add_argument('-p',
                        '--polygons',
                        nargs='?',
                        type=int,
                        help='Polygons to load for the ROIs')
    parser.add_argument(
        '-f',
        '--force',
        action='store_true',
        help='Force creating the pickle files overwriting existing files')
    args = parser.parse_args()

    print('--' * 10)
    rois_available = als.display_rois_available()
    print('--' * 10)

    i = args.iteration if args.iteration else 14400000
    r = args.replicates if args.replicates else 1
    if args.section == 'isthmus':
        rois = 1
        min_d, max_d = 14, 30
    else:
        rois = 2
        min_d, max_d = 160, 320

    load_rois = rois_available[args.polygons]

    print('      path:', args.path)
    print('   section:', args.section)
    print('   version:', args.version)
    print(' iteration:', i)
    print('replicates:', r)
    print('   ROIS id:', args.polygons if args.polygons else '---', ', file:',
          load_rois)

    main(args.section,
         args.path,
         i,
         r,
         load_rois,
         min_d,
         max_d,
         args.version,
         streaming=True,
         force_overwrite=args.force)
