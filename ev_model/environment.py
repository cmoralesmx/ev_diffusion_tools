import numpy as np
import math
import random
import pickle
from math import radians
import itertools
import time

from ev_model import elements, environment, persistency
from ev_model.utilities import geometry as evmetry

sections = {
    'utj': 'resources/img_polygons/3-1to3-5_L2-3_bl-utj.pickle',
    'isthmus': 'resources/img_polygons/3-1to3-5_L2-3_tl-isth.pickle',
    'ia-junction': 'resources/img_polygons/3-1to3-5_L2-3_c.pickle',
    'ampulla': 'resources/img_polygons/3-1to3-5_L2-3_tr-amp_2.pickle',
    'infundibulum': 'resources/img_polygons/3-1to3-5_L2-3_br-inf.pickle'
}

# identify the coordinates for the segments, the midpoints,
# the directions, the lengths 
def process_edge(segments, p1, p2, avg_cell_diameter, section_id=0, 
        inside_another=False, debug_log=False):
    """
    Given that the normals 'point out' from a CCW described shape, we restrict 
    their use to shapes laying inside other shapes.
    Note to self: The binormal points to the inside of the same shape

    inputs:
    segments - list of segments
    p1, p2 - segment endpoints
    avg_cell_diameter - the average cell diameter to use
    section_id -
    internal - if the shape is contained by another we use the normal, otherwise,
                we use the binormal
    """
    segments_text = []
    distance_between_points = evmetry.two_points_distance(p1, p2)
    if debug_log:
        segments_text.append(f'processing edge p1:({p1[0]:-06.3f}, {p1[1]:-06.3f}), p2({p2[0]:-06.3f}, {p2[1]:-06.3f}) dBwPts = {distance_between_points:-06.3f}\n')

    # is this section wider/larger than the avg cell diameter?
    if distance_between_points >= avg_cell_diameter:
        source_p = p1
        # how many cells should we fit?
        n = distance_between_points / avg_cell_diameter
        cells_to_fit = int(np.ceil(n)) if n % 1 > 0.5 else int(n)
        target_cell_diameter = distance_between_points / cells_to_fit

        if debug_log:
            segments_text.append(f'fits {cells_to_fit} cells, target diameter: {target_cell_diameter}\n')
        p1p2 = p2 - p1

        if cells_to_fit >= 2:
            # for more than two cells, we divide the segment into n subsegments
            for st in range(cells_to_fit - 1):
                # compute the coordinates of the new end point for the subsegment
                dest_p = evmetry.project_point_on_direction(source=source_p,
                                        direction = p1p2, dist=target_cell_diameter)
                # compute the mid point between the source and the new end point
                mp = evmetry.middle_point(source_p, dest_p)
                # compute the normal for the sub segment, 
                # it must match the normal for the parent segment
                seg = evmetry.Segment(source_p[0], source_p[1], dest_p[0],
                                        dest_p[1], mp[0], mp[1], inside_another)
                segments.append(seg)
                if debug_log:
                    segments_text.append(f'  >>> segment {st} from p1:({source_p[0]}, {source_p[1]}) to p2:({dest_p[0]}, {dest_p[1]})\n')
                source_p = dest_p
            # for the last segment
            mp = evmetry.middle_point(source_p, p2)
            seg = evmetry.Segment(source_p[0], source_p[1], p2[0], p2[1],
                                    mp[0], mp[1], inside_another)
            segments.append(seg)
            if debug_log:
                segments_text.append(f'  >>> segment {st} from p1:({source_p[0]}, {source_p[1]}) to p2:({p2[0]}, {p2[1]})\n')
            source_p = p2
        else:
            # fit one cell
            dest_p = evmetry.project_point_on_direction(source=source_p,
                        direction = p1p2, dist=1./distance_between_points
                        * avg_cell_diameter)
            mp = evmetry.middle_point(source_p, dest_p)

            seg = evmetry.Segment(source_p[0], source_p[1], dest_p[0],
                                    dest_p[1], mp[0], mp[1], inside_another)
            segments.append(seg)
            if debug_log:
                segments_text.append(f'  >> segment A from p1:({source_p[0]}, {source_p[1]}) to p2:({dest_p[0]}, {dest_p[1]})\n')
            source_p = dest_p

            mp = evmetry.middle_point(source_p, p2)
            seg = evmetry.Segment(source_p[0], source_p[1], p2[0], p2[1],
                                    mp[0], mp[1], inside_another)
            segments.append(seg)
            # check_normal(seg, internal)
            if debug_log:
                segments_text.append(f'  >> segment B from p1:({source_p[0]}, {source_p[1]}) to p2:({dest_p[0]}, {dest_p[1]})\n')
    else:
        # only one cell fits in this section
        mp = evmetry.middle_point(p1, p2)
        seg = evmetry.Segment(p1[0], p1[1], p2[0], p2[1], mp[0], mp[1], inside_another)
        segments.append(seg)
        # check_normal(seg, internal)
        if debug_log:
            segments_text.append(f'  > segment from p1:({p1[0]}, {p1[1]}) to p2:({p2[0]}, {p2[1]})\n')

    if debug_log:
        with open('segments.log', 'a') as outfile:
            for line in segments_text:
                outfile.write(line)
    return segments

def identify_segments_in_polygons(polygons, avg_cell_diameter):
    segments = []
    n = len(polygons)
    for i in range(n):
        # we need to check if this polygon is inside another
        cur_pol = polygons[i]
        inside_another = True
        e_coords = np.array(cur_pol.exterior.coords)

        for idx in range(len(e_coords) - 1):
            p1 = e_coords[idx]
            p2 = e_coords[idx + 1]
            segments = process_edge(segments, p1, p2, avg_cell_diameter, i, inside_another=inside_another)
    return segments

def compute_source_points_locations_in_cell(p1, p2, largest_ev_diameter_um=0.4, debug_log=False):
    """
    Each cell can have many source points if its size is larger than the largest sized EV.
    This function computes these source points per cell.
    """
    direction = p2 - p1
    sps = []
    d = evmetry.two_points_distance(p1, p2)
    max_possible_source_points_in_cell = int(d / (1.1 * largest_ev_diameter_um))
    distance_between_source_points =  d / max_possible_source_points_in_cell
    if debug_log:
        print(f'direction: {direction}, d: {d:6.3f}, max_possible_source_points_in_cell: {(max_possible_source_points_in_cell-1 if max_possible_source_points_in_cell > 1 else max_possible_source_points_in_cell) }, distance_between_source_points: {distance_between_source_points}')
    if max_possible_source_points_in_cell > 1:
        for i in range(1, max_possible_source_points_in_cell):
            sps.append(evmetry.project_point_on_direction(source=p1, direction=direction, dist=i * distance_between_source_points))
    else:
        sps.append(evmetry.middle_point(p1, p2))
    return sps

def identify_cell_range(x, y, grid, grid_size):
    """
    Computes a logical grid of n by n cells. Each cell has a size or grid_size
    in each of its dimmensions.
    Upon receiving a pair of values for a co-ordinate, it estimates the cell
    that would store this co-ordinate and the eight nearest neighbouring cells.
    If the current cell is already stored in the grid, it would return the set
    of cells
    """
    current_cell = (x // grid_size, y // grid_size)
    min_x = int(current_cell[0] - 1 if current_cell[0] > 1 else current_cell[0])
    min_y = int(current_cell[1] - 1 if current_cell[1] > 1 else current_cell[1])
    # print(f'min_x:{min_x}, min_y:{min_y}')
    if current_cell in grid:
        return current_cell, [c for c in [(x,y)
                for x in range(min_x, min_x + 3)
                for y in range(min_y, min_y + 3)] if c in grid]
    else:
        return current_cell, [current_cell]

def filter_source_points(grid, largest_ev_diameter_um=0.4, debug_log=False):
    """
    Checks the source points in neighbouring cells. 
    If they are too close to each other discard one of them randomly
    """
    min_space = largest_ev_diameter_um * 2.5
    dropped_count = 0
    emptied = []
    for k in grid.keys():
        minx = int(k[0] if k[0] == 1 else k[0] - 1)
        miny = int(k[1] if k[1] == 1 else k[1] - 1)
        if debug_log:
            print('Will start checking at',minx, miny)
        neighborhoods = [(x,y) for x in range(minx, minx + 3) for y in range(miny, miny + 3)]
        for cellid in [c for c in neighborhoods if c in grid]:
            if len(grid[cellid]) > 2 :
                if debug_log:
                    print(f'would check the sp in the {len(grid[cellid])} cells at {cellid}')
                for comb in itertools.combinations(grid[cellid], 2):
                    if debug_log:
                        print(f'{comb[0]._id} locations: {comb[0].source_points_locations}')
                        print(f'{comb[1]._id} locations: {comb[1].source_points_locations}')
                    max_loops = len(comb[0].source_points_locations) + len(comb[1].source_points_locations)
                    for _ in range(max_loops):
                        r1 = [i for i in range(len(comb[0].source_points_locations))]
                        r2 = [j for j in range(len(comb[1].source_points_locations))]
                        if debug_log:
                            print(r1)
                            print(r2)
                        for pair in [(e1, e2) for e1 in r1 for e2 in r2]:
                            if debug_log:
                                print(f'e1:{pair[0]}, e2:{pair[1]}')
                            if (np.dot(comb[0].unit_normal, comb[1].unit_normal) 
                                    <= 0 and evmetry.two_points_distance(
                                            comb[0].source_points_locations[pair[0]], 
                                            comb[1].source_points_locations[pair[1]])
                                            < min_space):
                                if debug_log:
                                    print(f'Will drop one of these indices: {pair[0]} or {pair[1]}')
                                dropped_count += 1
                                if random.random() < 0.5:
                                    del(comb[0].source_points_locations[pair[0]])
                                    if len(comb[0].source_points_locations) == 0:
                                        emptied.append( (cellid, comb[0]) )
                                    if debug_log:
                                        print(f'Chosen: {pair[0]}')
                                else:
                                    del(comb[1].source_points_locations[pair[1]])
                                    if len(comb[1].source_points_locations) == 0:
                                        emptied.append( (cellid, comb[1]) )
                                    if debug_log:
                                        print(f'Chosen: {pair[1]}')
                                break
                    if debug_log:
                        print('done with this pair')
                # checked = True
    print("Removed",dropped_count," potentially conflictive source points." \
            "As a result,", len(emptied),"secretory cells have no source points. These will become ciliary cells.")
    return grid, emptied

def transform_empty_secretories_in_ciliaries(grid, emptied):
    new_ciliaries = []
    for grid_id, empty_secretory in emptied:
        for secretory_id in range(len(grid[grid_id])):
            if grid[grid_id][secretory_id]._id == empty_secretory._id:
                # make this a new ciliary and remove from secretories
                new_c = elements.CiliaryCell(empty_secretory._id,
                    empty_secretory.p1[0], empty_secretory.p1[1],
                    empty_secretory.p2[0], empty_secretory.p2[1],
                    empty_secretory.mp[0], empty_secretory.mp[1],
                    empty_secretory.direction_length, empty_secretory.internal)
                new_ciliaries.append(new_c)
                del(grid[grid_id][secretory_id])
                break
    return grid, new_ciliaries

def updated_list_of_secretory_cells(grid):
    return list(grid.values())

def produce_cells(segments, probabilistic = False, threshold_secretory = 0.8,
        largest_ev_radius_um=0.2, debug_log=False, grid_size=5):
    """
    This function uses the provided segments to identify/assign the type of cell it will become.

    Edge case scenarios:
    - Two mutually facing secretory cells too close (d < 2xlargest_ev_radius_um) to each other will send EVs out of bounds.
        Mitigation:
            1 Meassure the distance between segments.
            2 keep track of what segment originated the new Secretory cell
            3 Before a new Secretory cell is introduced, check the distance to the other segments already originating a Secretory cell
            4 If within the minimum distance, avoid creating a Secretory cell
    """
    # use the identified segments to assign the cells
    cells = []
    secretory = []
    ciliary = []
    grid = {}

    cell_id = 0

    largest_ev_diameter_um = largest_ev_radius_um * 2
    shortest_secretory_separation = largest_ev_diameter_um * 2.5

    for segment in segments:
        passes_secretory_checks = True
        # if the distance between the segment's points is shorter than the
        # radius of the largest EV we must create a ciliary cell
        if segment.distance_p1p2 < largest_ev_diameter_um:
            passes_secretory_checks = False
        else:
            # if two segments facing opposing directions are closer than
            # 2xlargest_ev_diameter_um, only one of them can become a Secretory cell
            current_cell, cell_range = identify_cell_range(segment.mp[0],
                segment.mp[1], grid, grid_size)

            if current_cell not in grid:
                # the cell can safely become a secretory cell
                pass
            else:
                # check secretory cells in close range
                for cell_id in cell_range:
                    for secretory_cell in grid[cell_id]:
                        d = evmetry.two_points_distance(segment.location,
                            secretory_cell.location)
                        dp = np.dot(segment.unit_normal, secretory_cell.unit_normal)
                        if dp < 1. and d < shortest_secretory_separation:
                            passes_secretory_checks = False
                            break

            if probabilistic and random.random() < threshold_secretory:
                    passes_secretory_checks = False

        if passes_secretory_checks:
            cell = elements.SecretoryCell(cell_id, segment.p1[0], segment.p1[1],
                segment.p2[0], segment.p2[1], segment.mp[0], segment.mp[1],
                segment.direction_length, segment.internal)
            cell.source_points_locations = list(
                    compute_source_points_locations_in_cell(cell.p1, cell.p2,
                largest_ev_diameter_um))

            if debug_log:
                print(f'''Cell {cell._id} had {len(cell.source_points_locations)} 
                    segment.urce point locations. Is it internal? {segment.internal}''')
            secretory.append(cell)
            if current_cell not in grid:
                grid[current_cell] = list()
            grid[current_cell].append(cell)
        else:
            cell = elements.CiliaryCell(cell_id, segment.p1[0], segment.p1[1],
                segment.p2[0], segment.p2[1], segment.mp[0], segment.mp[1],
                segment.direction_length, segment.internal)
            ciliary.append(cell)

        cells.append(cell)
        cell_id += 1
        cell = None

    return cells, secretory, ciliary, grid

def at_safe_distance_to_near_secretory_cells(secretory, segment, shortest_secretory_separation):
    # checks secretory cells already allocated in the region
    # if two segments facing opposing directions are closer than
    # 2xlargest_ev_diameter_um, only one of them can become a Secretory cell
    _, cell_range = secretory.nearest_cells(segment.mp[0], segment.mp[1])
    for cell_id in cell_range:
        for secretory_cell in secretory[cell_id].values():
            d = evmetry.two_points_distance(segment.location, secretory_cell.location)
            dp = np.dot(segment.unit_normal, secretory_cell.unit_normal)
            if dp < 1. and d < shortest_secretory_separation:
                return False
    return True
    
def produce_cell_grids(segments, probabilistic = False, threshold_secretory = 0.8,
        largest_ev_radius_um=0.2, debug_log=False, grid_size=5):
    """
    Tests a segment against the prevously allocated secretory cells to determine
    if a new Secretory cell can be produced. Otherwise, the new cell will be a 
    Ciliary cell.

    Edge case scenarios:
    - Two mutually facing secretory cells very close (d < 2xlargest_ev_radius_um)
      to each other will send EVs out of bounds.
      Mitigation:
            1 Meassure the distance between segments.
            2 keep track of what segment originated the new Secretory cell
            3 Before a new Secretory cell is introduced, check the distance to the
                other segments already originating a Secretory cell
            4 If within the minimum distance, avoid creating a Secretory cell
    """
    # use the identified segments to assign the cells
    all_cells = elements.GridStorage(grid_size)
    secretory = elements.GridStorage(grid_size)
    ciliary = elements.GridStorage(grid_size)

    cell_id = 0

    largest_ev_diameter_um = largest_ev_radius_um * 2
    shortest_secretory_separation = largest_ev_diameter_um * 2.5

    for segment in segments:
        passes_secretory_checks = False
        # if the distance between the segment's points is shorter than the
        # radius of the largest EV we must create a ciliary cell
        if segment.distance_p1p2 > largest_ev_diameter_um:
            if at_safe_distance_to_near_secretory_cells(secretory, segment,
                    shortest_secretory_separation):
                if probabilistic and random.random() > threshold_secretory:
                    passes_secretory_checks = True

        if passes_secretory_checks:
            cell = elements.SecretoryCell(cell_id, segment.p1[0], segment.p1[1],
                segment.p2[0], segment.p2[1], segment.mp[0], segment.mp[1],
                segment.direction_length, segment.internal)
            cell.source_points_locations = list(
                    compute_source_points_locations_in_cell(cell.p1, cell.p2,
                largest_ev_diameter_um))

            if debug_log:
                print(f'''Cell {cell._id} had {len(cell.source_points_locations)} 
                    source point locations. Is it internal? {segment.internal}''')
            secretory.store(cell)
        else:
            cell = elements.CiliaryCell(cell_id, segment.p1[0], segment.p1[1],
                segment.p2[0], segment.p2[1], segment.mp[0], segment.mp[1],
                segment.direction_length, segment.internal)
            ciliary.store(cell)

        all_cells.store(cell)
        cell_id += 1
        cell = None

    return all_cells, secretory, ciliary

def generate_source_points_from_cells(cells, debug_log=False):
    """
    Collects the source points as a unique list
    """
    source_points = {'location':[], 'direction':[]}

    cells_with_n_sps = {}

    for cell_idx in range(len(cells)):
        cell = cells[cell_idx]
        n_sps = len(cell.source_points_locations)
        if(debug_log):
            print(f'surce points in this cell: {n_sps}')
        for sp in cell.source_points_locations:
            source_points['location'].append(sp)
            source_points['direction'].append(cell.unit_normal)

        if n_sps not in cells_with_n_sps:
            cells_with_n_sps[n_sps] = []
        cells_with_n_sps[n_sps].append(cell_idx)

    max_spts_per_cell = np.max([k for k in cells_with_n_sps.keys()])
    print(f'MAX NUMBER OF SOURCE POINTS PER CELL is {max_spts_per_cell}, {len(cells_with_n_sps[max_spts_per_cell])} cells have such many cells')
    if 0 in cells_with_n_sps:
        print('Cells with 0 source points:', len(cells_with_n_sps[0]),'indices:',
                ', '.join([str(n) for n in cells_with_n_sps[0]]))
    else:
        print('No cells with 0 source points')
    return source_points, max_spts_per_cell

def generate_evs_from_source_points(max_evs, secretory, new_evs_interval_seconds, debug_log=True):
    """
    Generates the N requested EVs.
    We create a list of ids for the secretory cells which gets shuffled.
    Then, we randomly choose a cell to produce a new EV by iterating from this list.
    The actual source point is selected randomly
    """
    source_point_indexes = [i for i in range(len(secretory))]
    random.shuffle(source_point_indexes)
    evs = []
    j = 0

    for i in range(max_evs):
        # generate a new EV in order to allocate the values non-dependent on the current/previous position
        ev = elements.EV(_id=i)
        sc_idx = source_point_indexes[j]
        sc = secretory[sc_idx]

        # we randomly select from the source points available in this cell
        sp_idx = random.randrange(len(sc.source_points_locations))
        starting_point = sc.source_points_locations[sp_idx]
        displacement_distance = ev.radius_um * ((0.4 * random.random()) + 0.1)

        ev.x_1=starting_point[0]
        ev.y_1=starting_point[1]

        # To prevent false positive collisions, we displace the current starting_point by the length of the EV radius
        # and on the same direction of the unit normal of the base source point        
        new_point = evmetry.project_point_on_direction(starting_point, sc.unit_normal, displacement_distance)
        if debug_log:
            print(f'    new ev_radius:{ev.radius_um}, displacement_distance:{displacement_distance}], new_point:[x:{new_point[0]}, y:{new_point[1]}]')

        # assign the location for the EV
        ev.x = float(new_point[0])
        ev.y = float(new_point[1])
        # The base direction for the EV matches the unit normal of the source point
        ev.base_direction_x = sc.unit_normal[0]
        ev.base_direction_y = sc.unit_normal[1]

        sc.time_to_next_secretion_attempt = new_evs_interval_seconds
        sc.last_source_point_secreting = sp_idx
        evs.append(ev)

        # prevent evs with crazy velocities
        if np.abs(ev.vx) >= 10 or np.abs(ev.vy) >= 10:
            print(f'warning: {i} dirx:{starting_point.dirx} diry:{starting_point.diry}')

        j += 1

    return evs

def generate_evs_in_circle(n, circle_radius, random_direction=False):
    """
    inputs:
    n - number of EVs to generate
    circle_radius - desired radius of the circle to use
 
    output:
    evs - a list of EV objects
    """
    evs = []
    for i in range(n):
        # compute a random initial position
        a = random.random() * 2 * np.pi
        r = circle_radius * np.sqrt(random.random())
        x = r * np.cos(a)
        y = r * np.sin(a)
        magnitude = np.sqrt(x*x + y*y)

        if random_direction:
            # compute a random angle
            a = random.random() * 2 * np.pi
            direction = np.array([np.cos(a), np.sin(a)])
        else:
            direction = np.array([x/magnitude, y/magnitude])

        ev = elements.EV(_id=i, x_1=x, y_1=y)
        ev = evmetry.project_point_on_direction(ev, direction[0], direction[1], ev.radius_um)
        # assign initial velocity
        ev.vx = 0 * direction[0]
        ev.vy = 0 * direction[1]
        evs.append(ev)
    return evs

# secretory, ciliary, collection_ciliaries, collection_secretories
def export_cell_distributions(base_filename, threshold_secretory, cells=None,
                                secretory=None, ciliary=None, grid=None):
    s = int(threshold_secretory * 100)
    c = str(100 - s)
    s = str(s)
    filename = f'resources/distributions/{base_filename}_{s}-{c}_s-c_'
    if secretory:
        flnm = f'{filename}secretory.pickle'
        with open(flnm, 'bw') as pickle_target:
            pickle.dump(secretory, pickle_target)
            print(f'Exported file: {flnm}')
    if ciliary:
        flnm = f'{filename}ciliary.pickle'
        with open(flnm, 'bw') as pickle_target:
            pickle.dump(ciliary, pickle_target)
            print(f'Exported file: {flnm}')
    if cells:
        flnm = f'{filename}allCells.pickle'
        with open(flnm, 'bw') as pickle_target:
            pickle.dump(cells, pickle_target)
            print(f'Exported file: {flnm}')
            
def export_cell_distribution_grids(base_filename, threshold_secretory, cells=None,
        secretory=None, ciliary=None):
    s = int(threshold_secretory * 100)
    c = str(100 - s)
    s = str(s)
    filename = f'resources/distributions/{base_filename}_{s}-{c}_s-c_'
    if secretory:
        flnm = f'{filename}secretory_grid.pickle'
        with open(flnm, 'bw') as pickle_target:
            pickle.dump(secretory, pickle_target)
            print(f'Exported file: {flnm}')
    if ciliary:
        flnm = f'{filename}ciliary_grid.pickle'
        with open(flnm, 'bw') as pickle_target:
            pickle.dump(ciliary, pickle_target)
            print(f'Exported file: {flnm}')
    if cells:
        flnm = f'{filename}allCells_grid.pickle'
        with open(flnm, 'bw') as pickle_target:
            pickle.dump(cells, pickle_target)
            print(f'Exported file: {flnm}')

def load_cell_distributions(base_filename, pct_sec, pct_cil):
    s = str(pct_sec)
    c = str(pct_cil)
    filename = f'resources/distributions/{base_filename}_{s}-{c}_s-c_'
    with open(f'{filename}allCells.pickle', 'rb') as pickle_source:
        loaded_cells = pickle.load(pickle_source)

    with open(f'{filename}secretory.pickle', 'rb') as pickle_source:
        loaded_secretory = pickle.load(pickle_source)

    with open(f'{filename}ciliary.pickle', 'rb') as pickle_source:
        loaded_ciliary = pickle.load(pickle_source)

    return loaded_cells, loaded_secretory, loaded_ciliary
    
def load_cell_distribution_grids(base_filename, pct_sec, pct_cil):
    s = str(pct_sec)
    c = str(pct_cil)
    filename = f'resources/distributions/{base_filename}_{s}-{c}_s-c_'
    with open(f'{filename}allCells_grid.pickle', 'rb') as pickle_source:
        loaded_cells = pickle.load(pickle_source)

    with open(f'{filename}secretory_grid.pickle', 'rb') as pickle_source:
        loaded_secretory = pickle.load(pickle_source)

    with open(f'{filename}ciliary_grid.pickle', 'rb') as pickle_source:
        loaded_ciliary = pickle.load(pickle_source)
    return loaded_cells, loaded_secretory, loaded_ciliary

def load_polygons(section):
    global sections
    if section not in sections:
        raise ValueError('section ',section,'does not exist in sections',sections)
    else:
        source = sections[section]
        with open(source, 'rb') as poly_source:
            polygons = pickle.load(poly_source)
        return polygons

def setup_experiment(version, working_section, degrees_of_freedom = 2,
                    dt=0.005, avg_cell_diameter=4.5, threshold_secretory=0.5,
                    target_spts_per_cell=11, new_evs_interval_seconds=20,
                    new_evs_threshold=0.97, n=0, debug_log=False,
                    plot_distributions=False, return_distributions=False):

    global sections

    start = time.time()
    s = time.time()
    base_filename = sections[working_section][len('resources/img_polygons/'):- len('.pickle')]
    print(f'Base filename: {base_filename}')
    with open(sections[working_section], 'rb') as pickled_source:
        polygons = pickle.load(pickled_source)

    segments = environment.identify_segments_in_polygons(polygons, avg_cell_diameter)
    print(f'Computed {len(segments)} segments in all polygons. Elapsed {time.time()-s:.2f}')

    s = time.time()
    # generate the hypothetical cell distributions
    cells, secretory, ciliary, grid = environment.produce_cells(segments, True,
                    threshold_secretory, grid_size=np.ceil(avg_cell_diameter))
    print(f'OLD METHOD: Generated {len(cells)} cells. The distribution is: {len(secretory)}'
        ' secretory and {len(ciliary)} ciliary. Elapsed {time.time()-s:.2f}')
    
    s = time.time()
    cells_g, secretory_g, ciliary_g = environment.produce_cell_grids(segments,
        True, threshold_secretory, grid_size=np.ceil(avg_cell_diameter))
    print(f'Generated {cells_g.total_elements} cells. The distribution is: '
        '{secretory_g.total_elements} secretory and {ciliary_g.total_elements}'
        'ciliary. Elapsed {time.time()-s:.2f}')

    # we must filter those source points too close to each other
    grid, emptied = environment.filter_source_points(grid)
    grid, new_ciliaries = environment.transform_empty_secretories_in_ciliaries(
        grid, emptied)
    ciliary.extend(new_ciliaries)
    secretory = environment.updated_list_of_secretory_cells(grid)
    cells = secretory + ciliary
    print(f'Updated lists: {len(cells)} cells, {len(secretory)} secretory, {len(ciliary)} ciliary')
    print(f'Filtered source points. Elapsed {time.time()-s:.2f}')

    # save the distributions for future reuse
    environment.export_cell_distributions(base_filename, threshold_secretory,
        cells, secretory, ciliary)
    environment.export_cell_distribution_grids(base_filename, threshold_secretory, 
        cells_g, secretory_g, ciliary_g)

    s = time.time()
    # the max EVs to generate is capped at len(secretory)
    if n >= len(secretory)/4:
        n = int(len(secretory)/4)
        print(f"WARNING: max EV introduced in initial state caped at 1/4 len(secretory), or: {n}")

    new_evs = environment.generate_evs_from_source_points(n, secretory,
                            new_evs_interval_seconds, debug_log=debug_log)
    if n > 1:
        print(f'EVs generated: {len(new_evs)}. Elapsed {time.time() - s}')
    else:
        print('No EVs were generated')

    # produce the initial states file 0.xml
    persistency.save_0xml_file(new_evs, n, dt, degrees_of_freedom, base_filename,
            target_spts_per_cell, version=version, boundaries=True,
            ev_ev_collisions=True, cells=cells, new_evs=True,
            new_evs_interval_seconds=new_evs_interval_seconds,
            new_evs_threshold=new_evs_threshold)
    end = time.time()
    print(f'Total time: {end-start:.2f} seconds')

    if plot_distributions:
        visualize_distributions(cells, secretory, ciliary, threshold_secretory)

    if return_distributions:
        return ciliary, secretory

def visualize_distributions(cells, secretory, ciliary, threshold_secretory,
                            save=False, scalebar=True, normals=False,
                            zoom_in=False, zoom=14, x1=940, x2=990, y1=600,
                            y2=640, evs=None, autoscale=False):
    from matplotlib import collections as mc
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm

    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from mpl_toolkits.axes_grid1 import inset_locator

    collection_secretories = np.array([s.coordinates for s in secretory])
    col_secretories = np.zeros((collection_secretories.shape[0], 2, 2))
    col_secretories[:,:,0] = collection_secretories[:,:,0]
    col_secretories[:,:,1] = collection_secretories[:,:,1]
    if len(ciliary) > 0:
        collection_ciliaries = np.array([c.coordinates for c in ciliary])
        col_ciliaries = np.zeros((collection_ciliaries.shape[0], 2, 2))
        col_ciliaries[:,:,0] = collection_ciliaries[:,:,0]
        col_ciliaries[:,:,1] = collection_ciliaries[:,:,1]

    # probability to human readable format
    s = int(threshold_secretory * 100)
    c = str(100 - s)
    s = str(s)

    _, ax1 = plt.subplots(1, 1,figsize=[20,10])

    pct_secretory = 100 / len(cells) * len(secretory)
    pct_ciliary = 100 / len(cells) * len(ciliary)
    ax1.set_title('Hypothetical cell distribution at the ampulla')
    ax1.add_collection(mc.LineCollection(col_secretories, color='red', linewidths=2))
    if len(ciliary) > 0:
        ax1.add_collection(mc.LineCollection(col_ciliaries, color='blue', linewidths=2))
    ax1.legend(ax1.collections[:2], [
        f'''Secretory cells: {len(secretory):,d} ({pct_secretory:.1f}%)
            Ciliary cells: {len(ciliary):,d} ({pct_ciliary:.1f}%)'''])

    if normals:
        collection_secretories_normals = np.array(
            [[s.mp, evmetry.project_point_on_direction(s.mp, s.normal, 1)] for s in secretory])
        col_sec_norm = np.zeros((collection_secretories_normals.shape[0], 2, 2))
        col_sec_norm[:,:,0] = collection_secretories_normals[:,:,0]
        col_sec_norm[:,:,1] = collection_secretories_normals[:,:,1]
        ax1.add_collection(mc.LineCollection(col_sec_norm, color='black', linewidths=1))
        if len(ciliary) > 0:
            collection_ciliaries_normals = np.array(
                [[c.mp, evmetry.project_point_on_direction(c.mp, c.normal, 1)] for c in ciliary])
            col_cil_norm = np.zeros((collection_ciliaries_normals.shape[0], 2, 2))
            col_cil_norm[:,:,0] = collection_ciliaries_normals[:,:,0]
            col_cil_norm[:,:,1] = collection_ciliaries_normals[:,:,1]
            ax1.add_collection(mc.LineCollection(col_cil_norm, color='black', linewidths=1))
    if autoscale:
        ax1.autoscale()
    else:
        if len(ciliary) > 0:
            max_x = np.max([np.max(col_ciliaries[:,:,0]), np.max(col_secretories[:,:,0])])
            max_y = np.max([np.max(col_ciliaries[:,:,1]), np.max(col_secretories[:,:,1])])
        else:
            max_x = np.max(col_secretories[:,:,0])
            max_y = np.max(col_secretories[:,:,1])
        ax1.axis([0, max_x * 1.1, 0, max_y * 1.1])
        ax1.set_aspect('equal')

    if scalebar:
        fontprops = fm.FontProperties(size=12)
        scalebar = AnchoredSizeBar(ax1.transData,
                               48, '100 um', 'lower right', 
                               pad=0.1,
                               frameon=False,
                               size_vertical=1,
                               fontproperties=fontprops)
        ax1.add_artist(scalebar)

    if zoom_in:
        zoom = 14
        axins = inset_locator.zoomed_inset_axes(ax1, zoom, loc=3) # zoom = 6
        axins.add_collection(mc.LineCollection(col_secretories, color='red', linewidths=2))
        if len(ciliary) > 0:
            axins.add_collection(mc.LineCollection(col_ciliaries, color='blue', linewidths=2))
        axins.set_title(f"Zoomed in x{zoom}")

        # sub region of the original image
        #x1, x2, y1, y2 = 940, 990, 600, 640
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)

        plt.xticks(visible=False)
        plt.yticks(visible=False)
        fontprops = fm.FontProperties(size=12)
        scalebar = AnchoredSizeBar(ax1.transData,
                                   22.5* zoom, '50 um', 'lower right', 
                                   pad=0.1,
                                   frameon=False,
                                   size_vertical=1,
                                   fontproperties=fontprops)
        axins.add_artist(scalebar)

        # draw a bbox of the region of the inset axes in the parent axes and
        # connecting lines between the bbox and the inset axes area
        inset_locator.mark_inset(ax1, axins, loc1=2, loc2=4, fc="none", ec="0.5")
        
    if evs is not None:
        # add the evs to the display
        for ev in evs:
            circ = plt.Circle(ev.location, ev.radius_um, color='k')
            ax1.add_artist(circ)
        
    # if save:
    #     fig.savefig('03_ampulla_hypothetical_cell_distributions_'+s+'-'+c+'_s-c_secretory_'+base_filename+'.png', bbox_inches='tight', dpi=300)

def reconstruct_boundaries_from_agents(all_cells, secretory, ciliary, avg_cell_diameter=5):
    """
    Every agent representing the boundaries is independent from each other.
    To compute the concentrations at distance from the boundaries we need to 
    reconstruct those shapes.
    To do so, we need to identify which agents have shared edges and collect them.
    avg_cell_diameter in um
    """
    edges = elements.GridStorage(avg_cell_diameter) 
    for cell in all_cells:
        edges.store(cell)
    