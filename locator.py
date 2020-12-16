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

    output:
    total area of cross-section in um^2
    total area of cross-section in mm^2
    list of areas per ROI
    list of counts of EVs per ROI
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
    q.put((replicate_id, not_valid))


def ev_in_roi_locator(rep_id, roi_id, evs_in_replicate, roi_coords):
    """
    Finds what EVs are inside a given ROI in a specific replicate
    """
    from shapely.geometry import Polygon, Point
    from shapely.prepared import prep

    # n_sizes = len(size_ranges)
    pol = Polygon(roi_coords)
    ppol = prep(pol)
    evs_here = []
    # print(f'rep {rep_id}, roi {roi_id}')
    for ev in evs_in_replicate:
        p = Point(ev['x'], ev['y'])
        if ppol.contains(p):
            evs_here.append(ev)
    q.put((roi_id, rep_id, evs_here))


def identify_valid_evs_multiprocess(poly_coords, evs_in_replicates):
    n_replicates = len(evs_in_replicates)
    print('Will check', n_replicates, 'replicates. \nProcessing, wait...')

    processes = []
    for rep_id in range(n_replicates):
        # iterate the list of replicates, each is handled by 1 processor
        p = Process(target=not_valid_ev_locator,
                    args=(rep_id, evs_in_replicates[rep_id], poly_coords))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()
    total_i = 0
    total_non_valid = 0
    for p in processes:
        res = q.get()
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
    n_polys = len(polys)
    n_replicates = len(evs_in_replicates)
    evs_in_roi_replicate_objects = [[list() for p in range(n_replicates)]
                                    for i in range(n_polys)]
    print('n_polys:', n_polys, 'n_replicates:', n_replicates)
    evs_in_roi_replicate_counts = np.zeros([n_polys, n_replicates])
    print(evs_in_roi_replicate_counts.shape)
    evs_in_roi_replicate_radius_age = list()

    # n_loops = n_polys//n_processors + n_polys % n_processors
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
            res = q.get()

            for ev in res[2]:
                evs_in_roi_replicate_objects[res[0]][res[1]].append(ev)

                evs_in_roi_replicate_counts[res[0], res[1]] += 1

                evs_in_roi_replicate_radius_age.append(
                    (res[0], res[1], ev['radius_um'], ev['age']))
        print('COMPLETED')
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


targets = {}

targets['infundibulum'] = {'buffer': 77}
targets['ia-junction'] = {'buffer': 80}
targets['utj'] = {'buffer': 34}
targets['isthmus'] = {'buffer': 32}
targets['ampulla'] = {'buffer': 334}

sections = ['isthmus', 'ampulla']

# definitions to use for the analysis
classes2 = {
    'ordering': [('narrow_end', '#30a2da'), ('wide_end', '#fc4f30'),
                 ('narrow_lumen', '#e5ae38'), ('wide_lumen', 'green')]
}
classes2['utj'] = {
    'narrow_end': [1, 2, 3, 7, 8, 13, 21],
    'wide_end': [4, 5, 6, 9, 14, 19],
    'narrow_lumen': [15, 16, 17, 18],
    'wide_lumen': [10, 11, 12, 20]
}
classes2['isthmus'] = {
    'narrow_end': [1, 2, 3, 4, 5, 6],
    'wide_end': [7, 8, 9],
    'narrow_lumen': [10, 11, 12, 13, 14, 15],
    'wide_lumen': [16, 17, 18, 19, 20, 21]
}
classes2['ia-junction'] = {
    'narrow_end': [1, 2, 3, 4, 5, 6, 7, 8],
    'wide_end': [9, 10, 11, 12, 13, 14, 15, 16],
    'narrow_lumen': [17, 18, 19, 20, 21, 22, 23, 24],
    'wide_lumen': [25, 26, 27, 28, 29, 30, 31, 32]
}
classes2['ampulla'] = {
    'narrow_end': [1, 2, 3, 4, 5, 6],
    'wide_end': [7, 8, 9, 10, 11, 12],
    'narrow_lumen': [13, 14, 15, 16, 17, 18],
    'wide_lumen': [19, 20, 21, 22, 23, 24]
}

# deffinitions for v3 of the analysis
# color palete from https://learnui.design/tools/data-color-picker.html#palette
classes3 = {
    'ordering': [('narrow_end', '#003f5c'), ('narrow_lumen', '#444e86'),
                 ('lumen_centre', '#955196'), ('near_epithelium', '#dd5182'),
                 ('far_epithelium', '#ff6e54'), ('wide_end', '#ffa600')]
}
classes3['isthmus'] = {
    'narrow_end': [1, 2, 4, 5, 6, 7],
    'wide_end': [],
    'narrow_lumen': [9, 10, 11, 12, 13, 14],
    'lumen_centre': [17, 18],
    'near_epithelium': [1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14],
    'far_epithelium': [3, 8, 15, 16, 19, 20, 21]
}
classes3['ampulla'] = {
    'narrow_end': [1, 2, 3, 4, 5, 6],
    'wide_end': [7, 8, 9, 10, 11, 12],
    'narrow_lumen': [13, 14, 15, 16, 17, 18],
    'lumen_centre': [19, 20, 21, 22, 23, 24],
    'near_epithelium': [32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
    'far_epithelium': [25, 26, 27, 28, 29, 30, 31]
}

# select which verion of the analysis to use
classes = copy.deepcopy(classes2)


def pickle_data_to_compressed_file(data, name, version, iteration):
    rao = './resources/analysis/output/'
    vni = f'{version}_{name}_{iteration}'
    with bz2.BZ2File(f'{rao}{vni}.pickle.bz2', 'wb') as compressed_output_file:
        pickle.dump(data, compressed_output_file)
        print(f'data in {name} saved to {rao}{vni}.pickle.bz2')


def main(section, base_path, iteration, replicates, load_rois, min_distance,
         max_distance, version, streaming, force_overwrite):
    distances_selected = [min_distance, max_distance]

    analysis_setup = als.prepare_analysis(section, targets, distances_selected,
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


m = Manager()
q = m.Queue()

if __name__ == '__main__':
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
