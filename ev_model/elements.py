import numpy as np
import math
import random
import sys
import xml.etree.cElementTree as cET
from collections import UserDict
from ev_model.utilities import geometry as evmetry

this = sys.modules[__name__]

# Dynamic viscosity of water
# source http://www.viscopedia.com/viscosity-tables/substances/water/
water_viscosity = {20: 0.0010016, 21: 0.0009775, 22: 0.0009544, 23: 0.0009321,
                    24: 0.0009107, 25: 0.00089}

# global variables
tempC = 20
this.const_Boltzmann = 1.3806504e-23
this.const_Temperature_K = tempC + 273.15
this.const_water_viscosity = water_viscosity[tempC]

class Cell(evmetry.Segment):
    def __init__(self, _id, p1x, p1y, p2x, p2y, mpx, mpy, diameter, internal):
        super(Cell, self).__init__(p1x, p1y, p2x, p2y, mpx, mpy, internal)
        self._id = _id
        self.radius = diameter / 2.
        self.diameter = diameter
        self.name = None

    @property
    def coordinates(self):
        return [(self.p1[0], self.p1[1]), (self.p2[0], self.p2[1])]

    def get_xml_description(self):
        return """
    <xagent>
        <name>{}</name>
        <id>{:d}</id>
        <p1_x>{:f}</p1_x>
        <p1_y>{:f}</p1_y>
        <p2_x>{:f}</p2_x>
        <p2_y>{:f}</p2_y>
        <x>{:f}</x>
        <y>{:f}</y>
        <direction_x>{:f}</direction_x>
        <direction_y>{:f}</direction_y>
        <direction_x_unit>{:f}</direction_x_unit>
        <direction_y_unit>{:f}</direction_y_unit>
        <direction_length>{:f}</direction_length>
        <normal_x>{:f}</normal_x>
        <normal_y>{:f}</normal_y>
        <unit_normal_x>{:f}</unit_normal_x>
        <unit_normal_y>{:f}</unit_normal_y>
        <normal_length>{:f}</normal_length>
    </xagent>""".format(self.name, self._id, self.p1[0], self.p1[1],
        self.p2[0], self.p2[1], self.mp[0], self.mp[1],
        self.direction[0], self.direction[1], self.direction_unit[0],
        self.direction_unit[1], self.direction_length, self.normal[0],
        self.normal[1], self.unit_normal[0], self.unit_normal[1],
        self.normal_length)

class ExtendedCell():
    def __init__(self, cell):
        self._id = int(cell['id'])
        self.cell = cell
        self.x = cell['x']
        self.y = cell['y']
        self.position = np.array([self.x, self.y])
        self.p1 = np.array([float(cell['p1_x']), float(cell['p1_y'])])
        self.p2 = np.array([float(cell['p2_x']), float(cell['p2_y'])])
        self.p3 = dict()
        self.p4 = dict()
        self.AB = dict()
        self.AD = dict()
        self.BC = dict()
        self.CD = dict()
        self.C1 = dict()
        self.C2 = dict()
        self.C3 = dict()
        self.C4 = dict()
        self.distances_computed = list()
        self.inner = np.vstack([self.p1, self.p2])
        self.outer = dict()
        self.unit_normal = np.array([float(cell['unit_normal_x']),
            float(cell['unit_normal_y'])])
        for tag, value in cell.items():
            if tag not in ['id', 'x', 'y', 'p1_x', 'p1_y', 'p2_x', 'p2_y',
                    'unit_normal_y', 'unit_normal_y']:
                setattr(self, tag, value)

    def compute_points_at_distance(self, distance):
        p4 = self.p1 + self.unit_normal * distance
        p3 = self.p2 + self.unit_normal * distance
        self.add_points_at_distance(distance, p3, p4)

    def add_points_at_distance(self, distance, p3, p4):
        self.p3[distance] = p3
        self.p4[distance] = p4
        self.outer[distance] = np.vstack([p3, p4])
        # works for CW defined shapes, not CCW
        A, B, C, D = self.p1, self.p2, p3, p4
        AB = B - A
        AB[1] *= -1
        AD = D - A
        AD[1] *= -1
        BC = C - B
        BC[1] *= -1
        CD = D - C
        CD[1] *= -1
        self.AB[distance], self.AD[distance] = AB, AD
        self.BC[distance], self.CD[distance] = BC, CD

        self.C1[distance] = -1. * self.yx_plus_xy(AB, A)
        self.C2[distance] = -1. * self.yx_plus_xy(AD, A)
        self.C3[distance] = -1. * self.yx_plus_xy(BC, B)
        self.C4[distance] = -1. * self.yx_plus_xy(CD, C)
        self.distances_computed.append(distance)

    def yx_plus_xy(self, a, b):
        return a[1] * b[0] + a[0] * b[1]

    def within_distance(self, distance, point):
        if not isinstance(point, np.ndarray):
            raise TypeError("m is not a Numpy Array")
        D1 = self.yx_plus_xy(self.AB[distance], point) + self.C1[distance]
        D2 = self.yx_plus_xy(self.AD[distance], point) + self.C2[distance]
        D3 = self.yx_plus_xy(self.BC[distance], point) + self.C3[distance]
        D4 = self.yx_plus_xy(self.CD[distance], point) + self.C4[distance]
        return 0 >= D1 and 0 >= D4 and 0 <= D2 and 0 >= D3

    def __str__(self):
        return f'Extended cell {self._id} at ({self.x:.2f} {self.y:.2f})'

    def __repr__(self):
        return f'Extended cell {self._id} at ({self.x:.2f} {self.y:.2f})'

class SecretoryCell(Cell):
    def get_xml_description(self, target_source_points):
        spts = np.vstack((self.source_points, np.zeros([target_source_points
            - self.source_points.shape[0], 2])))
        return f"""
    <xagent>
        <name>{self.name}</name>
        <id>{self._id:d}</id>
        <p1_x>{self.p1[0]:f}</p1_x>
        <p1_y>{self.p1[1]:f}</p1_y>
        <p2_x>{self.p2[0]:f}</p2_x>
        <p2_y>{self.p2[1]:f}</p2_y>
        <x>{self.mp[0]:f}</x>
        <y>{self.mp[1]:f}</y>
        <direction_x>{self.direction[0]:f}</direction_x>
        <direction_y>{self.direction[1]:f}</direction_y>
        <direction_x_unit>{self.direction_unit[0]:f}</direction_x_unit>
        <direction_y_unit>{self.direction_unit[1]:f}</direction_y_unit>
        <direction_length>{self.direction_length:f}</direction_length>
        <normal_x>{self.normal[0]}</normal_x>
        <normal_y>{self.normal[1]}</normal_y>
        <unit_normal_x>{self.unit_normal[0]}</unit_normal_x>
        <unit_normal_y>{self.unit_normal[1]}</unit_normal_y>
        <normal_length>{self.normal_length}</normal_length>
        <time_to_next_secretion_attempt>{self.time_to_next_secretion_attempt}</time_to_next_secretion_attempt>
        <source_points>{len(self.source_points_locations)}</source_points>
        <source_points_xs>{'f,'.join([str(n) for n in spts[:,0]])+'f'}</source_points_xs>
        <source_points_ys>{'f,'.join([str(n) for n in spts[:,1]])+'f'}</source_points_ys>
        <last_source_point_secreting>{self.last_source_point_secreting}</last_source_point_secreting>
    </xagent>"""

    def __init__(self, _id, p1x, p1y, p2x, p2y, mpx, mpy, diameter, internal):
        super(SecretoryCell, self).__init__(_id, p1x, p1y, p2x, p2y, mpx, mpy, diameter, internal)
        self.name = 'SecretoryCell'
        self.source_points_locations = []
        self.time_to_next_secretion_attempt = 0
        self.last_source_point_secreting = -1

    @property
    def source_points(self):
        return np.array(self.source_points_locations)

class CiliaryCell(Cell):
    def get_xml_description(self):
        return super(CiliaryCell, self).get_xml_description()

    def __init__(self, _id, p1x, p1y, p2x, p2y, mpx, mpy, diameter, internal):
        super(CiliaryCell, self).__init__(_id, p1x, p1y, p2x, p2y, mpx, mpy,
            diameter, internal)
        self.name = 'CiliaryCell'

class EV:
    # to convert:
    # N square metres to square micrometre: multiply N * 1e+12
    # N metres to micrometres: multiply N * 1e+6
    # N micrometres to metres: divide N / 1e+6       
    def mass_from_radius(self):
        """
        Computes the mass as a function of the volume of the EV. Assuming every EV
        is of spherical shape and has the same material density, we can compute the mass (m)
        from the density (p) and volume (v), as follows:  p = m/v -> m = p * v
        To do so, we take the (estimated) mass_of_a_single_100um_EV_in_g as 5.47990700763866E-17
        The volume of a sphere is a function of its diameter v = 4/3 pi * r

        Some values pre-computed for efficiency:
        # the constant part of the sphere volume formula: (4/3)*Pi = 4.1887902047863905
        # the volume_100_nm_ev = 523598.7755982988 nm^3
        # the mass_per_volume_unit = mass_of_a_single_100um_EV_in_g / volume_100_nm_ev
        # the mass_per_volume_unit = 1.0465851455395313e-22

        # (4/3)*PI*mass_per_volume_unit = 4.38392560611093E-22
        mass_p_vol_u_x_4div3_pi = 4.38392560611093E-22
        """
        # Precomputed values
        pi_times_4div3 = 4.1887902047863905
        # mass_per_volume_unit = 1.0465851455395313e-22
        mass_p_vol_u_x_4div3_pi = 4.38392560611093E-22

        # compute the sphere's volume
        self.volume = pi_times_4div3 * self.radius_nm**3
        # compute
        return mass_p_vol_u_x_4div3_pi * self.radius_nm**3

    def __init__(self, **kwargs):
        self._id = kwargs.get('_id', 0)
        self.x_1 = kwargs.get('x_1', 0)  # holds the previous position for this EV
        self.y_1 = kwargs.get('y_1', 0)
        self.vx = kwargs.get('vx', 0)
        self.vy = kwargs.get('vy', 0)
        # placeholders for the values generated for brownian motion
        self.bm_vx = kwargs.get('bm_vx', 0)
        self.bm_vy = kwargs.get('bm_vy', 0)
        # If no location is provided at creation time, the EV is located at the origin and is at rest
        self.x = kwargs.get('x', 0)
        self.y = kwargs.get('y', 0)
        self.base_direction_x = 0
        self.base_direction_y = 0
        self.colour = '()'
        self.radius_um = kwargs.get('radius_um', random.randrange(40, 160, 10) / 1000)
        self.radius_nm = self.radius_um * 1000
        self.radius_m = self.radius_um / 1e+6
        self.mass_g = self.mass_from_radius()
        self.mass_kg = self.mass_g * 1000
        self.mass_ag = self.mass_g / 1e-18

        # Compute the radius-dependent diffusion coefficient defined as
        # D = (boltzmann tempK)/(6 pi fluidViscosity particleRadius)
        self.diffusion_rate_m = ((this.const_Boltzmann * this.const_Temperature_K)
                                 / (6. * math.pi * this.const_water_viscosity * self.radius_m))
        self.diffusion_rate_um = self.diffusion_rate_m * 1e+12
        self.velocity_ums = 0
        self.velocity_ms = 0

    def compute_msd(self, dt, dof):
        # Here, we compute the mean squared displacement which is diffusion-rate dependant.
        # This value is an educated guess of the EV velocity, it appears to be the SD of the MSD
        # and works as the expected displacement of this particle
        self.velocity_ums = math.sqrt(2 * dof * self.diffusion_rate_um * dt)
        self.velocity_ms = math.sqrt(2 * dof * self.diffusion_rate_m * dt)

    @property
    def location(self):
        return np.array([self.x, self.y])

    def get_xml_description(self, degrees_of_freedom, time_in_initial_state):
        return """
    <xagent>
        <name>EV</name>
        <id>{:d}</id>
        <x>{}</x>
        <y>{}</y>
        <x_1>{}</x_1>
        <y_1>{}</y_1>
        <vx>{}</vx>
        <vy>{}</vy>
        <bm_vx>0</bm_vx>
        <bm_vy>0</bm_vy>
        <mass_ag>{:f}f</mass_ag>
        <radius_um>{:f}</radius_um>
        <diffusion_rate_um>{:.25f}</diffusion_rate_um>
        <diff_rate_um_x_twice_dof>{:.25f}</diff_rate_um_x_twice_dof>
        <velocity_ums>{:.10f}</velocity_ums>
        <closest_ev_id>-1</closest_ev_id>
        <closest_ev_distance>100</closest_ev_distance>
        <closest_secretory_cell_id>-1</closest_secretory_cell_id>
        <closest_secretory_cell_distance>100</closest_secretory_cell_distance>
        <closest_ciliary_cell_id>-1</closest_ciliary_cell_id>
        <closest_ciliary_cell_distance>100</closest_ciliary_cell_distance>
        <time_in_initial_state>{:f}</time_in_initial_state>
    </xagent>""".format(self._id,
        f'{self.x:f}' if self.x else '--',
        f'{self.y:f}' if self.y else '--',
        f'{self.x_1:f}' if self.x_1 else '--',
        f'{self.y_1:f}' if self.y_1 else '--',
        f'{self.vx:.10f}' if self.vx else '--',
        f'{self.vy:.10f}' if self.vy else '--',
        self.mass_ag, self.radius_um, self.diffusion_rate_um,
        self.diffusion_rate_um * 2 * degrees_of_freedom,
        self.velocity_ums, time_in_initial_state)

    def print_resume(self, cf=lambda val: f"{val:>9.4f}" if val else '  --.----'):
        print("""       x:{} y:{}
     v x:{} y:{}
    bm x:{} y:{}
   t-1 x:{} y:{}""".format(
        cf(self.x), cf(self.y), cf(self.vx), cf(self.vy),
        cf(self.bm_vx), cf(self.bm_vy), cf(self.x_1),  cf(self.y_1)))
        print('radius_um',self.radius_um)

class GridStorage():
    '''
    Represents a logical grid of n by n cells. Each cell has a size of grid_size
    in each of its dimmensions.
    For optimization purposes, the storage is similar to a sparse matrix.
    Only those cells with content are actually present in the underlying storage
    '''
    def __init__(self, grid_size, minx=None, maxx=None, miny=None, maxy=None):
        self.grid_size = grid_size
        args_list = [minx, maxx, miny, maxy]
        args_none = args_list.count(None)
        self.fixed_size = False
        # index tracks what cell is holding each element in data. We assume all the
        # elements have a unique ID, no clashing should occur.
        # key: element_id
        # value: {'cell_id':0, 'element':object} 
        self.index = {}
        # data is the actual cell-based storage
        self.grid = {}
        self.total_elements = 0
        if 0 < args_none < 4:
            raise ValueError('All range parameters must have values')
        elif args_none == 0:
            self.fixed_size = True
            self.minx, self.maxx, self.miny, self.maxy = args_list
            self.max_col = (self.maxx // self.grid_size) - 1
            self.max_row = (self.maxy // self.grid_size) - 1

    def what_cell(self, x, y):
        """
        Identifies what cell of the grid would contain the pair of co-ordinates
        """
        if self.fixed_size:
            current_cell = [None, None]
            x2 = (x - self.minx) // self.grid_size
            y2 = (y - self.miny) // self.grid_size
            current_cell[0] = 0 if x2 < 0 else x2 if x2 < self.max_col else self.max_col
            current_cell[1] = 0 if y2 < 0 else y2 if y2 < self.max_row else self.max_row

            return tuple(current_cell)
        else:
            current_cell = [x // self.grid_size, y // self.grid_size]
            return tuple(current_cell)

    def nearest_cells_for_coordinates(self, x, y):
        """
        Identifies the eight nearest neighbouring cells to the one that would
        contain the pair of co-ordinates received.
        If cell_index is already present in the grid we return the set
        of non-empy neighbouring cells to cell_index
        """
        cell_index = self.what_cell(x, y)
        return self.nearest_cells_for_cell_index(cell_index)

    def nearest_cells_for_cell_index(self, cell_index):
        min_x = int(cell_index[0] - 1)
        min_y = int(cell_index[1] - 1)
        return cell_index, [c for c in [(x,y)
                for x in range(min_x, min_x + 3)
                for y in range(min_y, min_y + 3)] if c in self.grid]

    def __getitem__(self, cell_idx):
        """
        Returns a dictionary holding the content of a whole cell from the grid.
        Bear in mind the cell stores only a list of element_id. The actual elements
        are stored in the index.
        """
        if cell_idx in self.grid:
            return self.grid[cell_idx]
        else:
            return None

    def store(self, element):
        """
        At the time of insertion, we identify the target cell for storing the element_id based
        on the coordinates ot the element being inserted.
        We also keep an index of the elements stored per cell on the grid. This is useful
        for fetching the elements directly instead of doing a lookup in the grid.
        """
        # cell level index, (x,y)
        cell_idx = (int(element.mp[0]), int(element.mp[1]))
        # grid level index, (x//grid_size, y//grid_size)
        target_cell = self.what_cell(cell_idx[0], cell_idx[1])

        # record what cell would be storing this element
        self.index[element._id] = {'cell_id': target_cell, 'element': element}

        if target_cell in self.grid:
            if cell_idx in self.grid[target_cell]:
                self.grid[target_cell][cell_idx].append(element._id)
            else:
                self.grid[target_cell][cell_idx] = [element._id]
        else:
            self.grid[target_cell] = {cell_idx: [element._id]}
        self.total_elements += 1

    def fetch_element_by_id(self, element_id):
        """
        get the element matching the given ID if such element exists in data
        """
        if element_id in self.index:
            return self.index[element_id]['element']
        else:
            print(f"The provided element_id: {element_id} does not exist in the data set.")
            return None

    def remove_element(self, element):
        """
        Checks if the element exists in the index and fetch the cell_id. Then, use the cell_id
        to find and remove the entry for this element in the grid. Finally, remove the entry
        from the same index.
        No further checking is needed.
        """
        # cell level index
        cell_idx = (int(element.x), int(element.y))
        # grid level index
        target_cell = self.what_cell(cell_idx[0], cell_idx[1])

        if target_cell in self.grid:
            if cell_idx in self.grid[target_cell]:
                for e in range(len(self.grid[target_cell][cell_idx])):
                    if self.grid[target_cell][cell_idx][e] == element:
                        del(self.grid[target_cell][cell_idx][e])
                self.total_elements -= 1
                return True
        return False

