import numpy as np
from ev_model import elements
import unittest

class GridStoreTestMethods(unittest.TestCase):

    def test_what_cell_dynamic(self):
        store = elements.GridStorage(5)
        all_coords = [(2,3), (10,5), (181,200), (10,0)]
        expectations = [(0,0), (2,1), (36,40), (2,0)]
        for (coords, expected) in zip(all_coords, expectations):
            with self.subTest(coords=coords, expected=expected):
                self.assertEqual(store.what_cell(coords[0], coords[1]), expected)
    
    def test_what_cell_fixed(self):
        fixed_store = elements.GridStorage(5, 10, 1000, 10, 500)
        all_coords = [(-100,-100), (-100,200), (19,600), (1010, 15)]
        expectations = [(0,0), (0,38), (1, 99), (199,1)]
        for (coords, expected) in zip(all_coords, expectations):
            with self.subTest(coords=coords, expected=expected):
                self.assertEqual(fixed_store.what_cell(coords[0], coords[1]), expected)
    
    def test_nearest_cells_empty_grid(self):
        empty_store = elements.GridStorage(5)
        all_coords = [(2,3), (10,5), (181,200)]
        expectations = [((0,0),[]), ((2,1),[]), ((36,40),[])]
        for (coords, expected) in zip(all_coords, expectations):
            with self.subTest(coords=coords, expected=expected):
                self.assertEqual(empty_store.nearest_cells(coords[0], coords[1]), expected)

    def test_nearest_cells_cur_cell_not_in_grid_no_neighbours(self):
        store = elements.GridStorage(5)
        store.data = {(0,0):{(1,1):[]}, (2,2):{(10,10):[]}, (2,3):{(14,17):[]}, 
            (3,2):{(15,13):[]}, (4,4):{(20,20):[]}}
        all_coords = [(50,50), (220,18), (10,0)]
        expectations = [((10,10),[]), ((44,3),[]), ((2,0),[])]
        for (coords, expected) in zip(all_coords, expectations):
            with self.subTest(coords=coords, expected=expected):
                self.assertEqual(store.nearest_cells(coords[0], coords[1]), expected)
    
    def test_nearest_cells_cur_cell_not_in_grid_some_neighbours(self):
        store = elements.GridStorage(5)
        store.data = {(0,0):{(1,1):[]}, (1,1):{(9,9):[]}, (2,3):{(14,17):[]}, 
            (3,2):{(15,13):[]}, (4,4):{(20,20):[]}}
        self.assertEqual(store.nearest_cells(12,12), ((2,2),[(1,1),(2,3),(3,2)]))

    def test_nearest_cells_cur_cell_in_grid_no_neighbours(self):
        store = elements.GridStorage(5)
        store.data = {(0,0):{(1,1):[]}, (1,1):{(9,9):[]}, (2,3):{(14,17):[]}, 
            (3,2):{(15,13):[]}, (4,4):{(20,20):[]}}
        self.assertEqual(store.nearest_cells(21,21), ((4,4),[(4,4)]))

    def test_nearest_cells_cur_cell_in_grid_some_neighbours(self):
        store = elements.GridStorage(5)
        store.data = {(0,0):{(1,1):[]}, (1,1):{(9,9):[]}, (2,3):{(14,17):[]}, 
            (3,2):{(15,13):[]}, (4,4):{(20,20):[]}}
        self.assertEqual(store.nearest_cells(21,25), ((4,5),[(4,4)]))

    def test_store_one_element(self):
        store = elements.GridStorage(5)
        ev1 = elements.EV(id=0, x=3, y=3)
        ev2 = elements.EV(id=1, x=1, y=1)
        ev3 = elements.EV(id=2, x=20,y=20)
        self.assertEqual(store.total_elements, 0)
        store.store(ev1)
        store.store(ev2)
        self.assertEqual(store.total_elements, 2)
        store.store(ev3)
        self.assertEqual(store.total_elements, 3)

    def test_get_item(self):
        store = elements.GridStorage(5)
        ev1 = elements.EV(id=0, x=3, y=3)
        ev2 = elements.EV(id=1, x=1, y=1)
        store.store(ev1)
        store.store(ev2)
        self.assertEqual(store[(0,0)], {(3,3):[ev1], (1,1):[ev2]})

if __name__ == '__main__':
    unittest.main()
        