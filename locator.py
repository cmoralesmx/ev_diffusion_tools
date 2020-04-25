import copy
from multiprocessing import Process, Manager, Queue, freeze_support
from pathlib import Path

from lxml import etree

import analysis_last_state as als


def f(rep_id, roi_id, evs_in_replicate, roi_coords, size_ranges):
    """
    Finds what EVs are inside a given ROI in a specific replicate
    """
    from shapely.prepared import prep

    n_sizes = len(size_ranges)
    pol = Polygon(roi_coords)
    ppol = prep(pol)
    evs_here = []
    for ev in evs_in_replicate:
        p = Point(ev['x'], ev['y'])
        size_id = None
        if ppol.contains(p):
            # identify the size range for the ev
            for s in range(n_sizes):
                if ev['radius_um'] < size_ranges[s]:
                    size_id = s
                    break
            evs_here.append((size_id, ev))
    q.put((roi_id, rep_id, evs_here))

def identify_evs_per_roi_multiprocess(polys, evs_in_replicates, sizes, n_processors):
    n_sizes = len(sizes)
    n_polys = len(polys)
    n_replicates = len(evs_in_replicates)
    evs_in_roi_replicate_objects = [[list() for p in range(n_replicates)] for i in range(n_polys)]
    evs_in_roi_size_replicate_counts = np.zeros([n_polys, n_sizes, n_replicates])

    n_loops = n_polys//n_processors + n_polys % n_processors
    print('Will perform',n_loops,'loops per replicate. \nProcessing:')
    for rep_id in range(n_replicates):
        processes = []
        # iterate the list of ROIs, each ROI is handled by 1 processor
        for roi_id in range(n_polys):
            p = Process(target=f, args=(rep_id, roi_id, evs_in_replicates[rep_id], polys[roi_id], sizes))
            processes.append(p)
            p.start()
            
        for p in processes:
            p.join()

        for p in processes:
            res = q.get()
            res[0] # roi_id
            res[1] # rep_id
            res[2] # [size_id, ev]
            for tup in res[2]:
                evs_in_roi_replicate_objects[res[0]][res[1]].append(tup[1])
                evs_in_roi_size_replicate_counts[res[0]][tup[0]][res[1]] += 1
                evs_in_roi_size_replicate.append((res[0],tup[0], res[1], tup[1]['radius_um'], tup[1]['age']))

    return evs_in_roi_replicate_objects, evs_in_roi_size_replicate_counts, evs_in_roi_size_replicate



def load():
    import numpy as np

    #def display_data_available():
    rois_available = als.display_rois_available()
    print('--'*10)
    experiments_available = als.display_experiments_available()
    analysis_setup = als.prepare_analysis(section, targets, distances_selected, base_path, replicates, target_iteration)
    user_rois_loaded, user_rois, analysis_setup['prep_polys'] = als.load_rois_from_file(rois_available[rois_to_load])

    analysis_setup = als.load_evs(analysis_setup)
    analysis_setup['evs_per_replicate'] = [len(li) for li in analysis_setup['evs_d_useful']]
    print(f'Evs per replicate (in full cross section) mean:{np.mean(analysis_setup["evs_per_replicate"])}, sd:{np.std(analysis_setup["evs_per_replicate"]):.2f}')

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

    sections = ['utj', 'isthmus', 'ia-junction', 'ampulla1', 'ampulla2']

    # deffinitions for v1 of the analysis
    classes1 = {'ordering': [('wide','#30a2da'), ('narrow','#fc4f30'), ('isolated','#e5ae38'), ('bottleneck','green')]}
    classes1['utj'] = {'narrow':[5,6,7,9], 'bottleneck':[8], 'wide':[1,2,3, 10, 11], 'isolated':[4]}
    classes1['isthmus'] = {'narrow':[3,4,5,7,8,9], 'bottleneck':[], 'wide':[1,2,6, 10], 'isolated':[]}
    classes1['ia-junction'] = {'narrow':[5,6,7,13,14,15,16,17], 'bottleneck':[4,8,18], 'wide':[1,2,3,9,10,11,12], 'isolated':[19]}
    classes1['ampulla1'] = {'narrow':[6,8,10,14,17,21,23,25,26,27], 'bottleneck':[7,9,11,12,13,15,16,19,20], 'wide':[1,2,3,4,5,18,24,28], 'isolated':[22]}
    classes1['ampulla2'] = {'narrow':[6,8,10,14,17,21,23,25,26,27], 'bottleneck':[7,9,11,12,13,15,16,19,20], 'wide':[1,2,3,4,5,18,24,28], 'isolated':[22]}

    # deffinitions for v2 of the analysis
    classes2 = {'ordering':[('narrow_end','#30a2da'), ('wide_end','#fc4f30'), ('narrow_lumen','#e5ae38'),('wide_lumen','green')]}
    classes2['utj'] = {'narrow_end':[1,2,3, 7,8, 13,21], 'wide_end':[4,5,6, 9,14,19], 'narrow_lumen':[15, 16, 17,18], 'wide_lumen':[10,11,12, 20]}
    classes2['isthmus'] = {'narrow_end':[1,2,3,4,5,6], 'wide_end':[7,8, 9], 'narrow_lumen':[10, 11, 12, 13, 14, 15], 'wide_lumen':[16,17,18,19,20,21]}
    classes2['ia-junction'] = {'narrow_end':[1,2,3,4,5,6,7,8], 'wide_end':[9,10,11,12,13,14,15,16], 'narrow_lumen':[17,18,19,20,21,22,23,24], 'wide_lumen':[25,26,27,28,29,30,31,32]}
    classes2['ampulla1'] = {'narrow_end':[1,2,3,4,5,6], 'wide_end':[7,8,9,10, 11, 12], 'narrow_lumen':[13,14,15,16,17,18], 'wide_lumen':[19, 20, 21, 22, 23, 24]}
    classes2['ampulla2'] = {'narrow_end':[1,2,3,4,5,6], 'wide_end':[7,8,9,10, 11, 12], 'narrow_lumen':[13,14,15,16,17,18], 'wide_lumen':[19, 20, 21, 22, 23, 24]}

    # select which verion of the analysis to use
    classes = copy.deepcopy(classes2)

def wtf():
    section = 'ampulla1'
    distances_selected = [160,320]
    
    replicates, target_iteration, base_path = 6, 2340000, 'D:\\iterations\\v33\\amp20estrous\\'

    rois_to_load = 0
    experiment_to_load = 7

    print(section, distances_selected, replicates, target_iteration, base_path)

    # base_polygon

    # WITHOUT distances MULTIPROCESS
    #[analysis_setup['evs_in_roi_replicate_objects'],
    #analysis_setup['evs_in_roi_size_replicate_counts'], 
    #analysis_setup['evs_in_roi_size_replicate']] = identify_evs_per_roi_multiprocess(user_rois, analysis_setup['evs_d_useful'], targets['sizes'], 4)

m = Manager()
#.Queue()

#if __name__ == '__main__':
if __name__ ==  '__main__':
    freeze_support()

    wtf()