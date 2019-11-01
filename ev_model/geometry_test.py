import numpy as np
import unittest
from ev_model.utilities import geometry

#  w0                   w1                  w2
square_walls = np.array([[[2., .7], [2., 2.]], [[2., 2.], [1., 2.]], [[1., 2.], [1., .7]]])
angled_walls = np.array([[[2.5, .7], [2., 2.]], [[2., 2.], [1., 2.]], [[1., 2.], [0.5, .7]]])

def test_reflect_vector():
    # x0 = np.array([[2., 1.316667]])
    # x1 = np.array([[2., 3.555574]])
    # wall_normal = np.array([[0, 1.]])
    return False

def test_project_point_on_direction():
    p0_1d = np.array([2.0])
    p0_2d = np.array([2.0, 8.0])
    p0_3d = np.array([2.0, 8.0, 1.0])
    p0_4d = np.array([2.0, 8.0, 5.0, 1.0])
    direction_1d = np.array([4.666667])
    direction_2d = np.array([0, 4.666667])
    direction_2d_unit = np.array([0., 1.])
    direction_4d = np.array([0, 4.666667, 5., 7.])
    dist = 0.040000
    new_point_2d = np.array([2., 8.04])
    new_point_3d = np.array([2., 8.04, 1.])

    try:
        _ = project_point_on_direction(p0_1d, direction_2d_unit, dist)
    except ValueError as err:
        if str(err) == 'Source has non-supported shape':
            print('T0 OK, got', err)
        else:
            print('T0 failed, there was another error:', err)
            return False

    try:
        _ = project_point_on_direction(p0_4d, direction_2d_unit, dist)
    except ValueError as err:
        if str(err) == 'Source has non-supported shape':
            print('T1 OK, got', err)
        else:
            print('T1 failed, there is another error:', err)
            return False

    try:
        _ = project_point_on_direction(p0_2d, direction_1d, dist)
    except ValueError as err:
        if str(err) == 'Direction has non-supported shape':
            print('T2 OK, got', err)
        else:
            print('T2 failed, there is another error:', err)
            return False

    try:
        _ = project_point_on_direction(p0_2d, direction_4d, dist)
    except ValueError as err:
        if str(err) == 'Direction has non-supported shape':
            print('T3 OK, got', err)
        else:
            print('T3 failed, there is another error:', err)
            return False

    try:
        _ = project_point_on_direction(p0_1d, direction_4d, dist)
    except ValueError as err:
        if str(err) == 'Direction has non-supported shape' or str(err) == 'Source has non-supported shape':
            print('T4 OK, got', err)
        else:
            print('T4 failed, there is another error:', err)
            return False

    result = project_point_on_direction(p0_2d, direction_2d, dist) == new_point_2d
    if result.all():
        print('T5 passed')
    else:
        print('T5 failed', result)
        return False

    result = project_point_on_direction(p0_3d, direction_2d_unit, dist) == new_point_3d
    if result.all():
        print('T6 passed')
    else:
        print('T6 failed', result)
        return False
    return True

def test_middle_point():
    # test 1
    p1 = np.array([-1, 2])
    p2 = np.array([3, -6])
    mp = np.array([1, -2])
    if not (middle_point(p1, p2) == mp).all():
        return False

    # test 2
    p1 = np.array([2.0, 5.666667])
    p2 = np.array([2.0, 10.333333])
    mp = np.array([2.0, 8.0])
    n_mp = middle_point(p1, p2)
    if not (n_mp == mp).all():
        return False

    return True

def test_two_points_normals_2d():
    # t1
    p1 = np.array([-1, 2])
    p2 = np.array([3, -6])
    two_pn = np.array([0, 0])
    if not (two_points_normals_2d(p1, p2) == two_pn).all():
        return False
    # t2
    p1 = np.array([2., 5.666667])
    p2 = np.array([2., 10.333333])
    two_pn = np.array([0., 1.])
    if not (two_points_normals_2d(p1, p2) == two_pn).all():
        return False
    return True

def all_tests():
    if not test_middle_point():
        raise ValueError('middle_point() failed')
    if not test_two_points_normals_2d:
        raise ValueError('two_points_normals_2d() failed')
    if not test_project_point_on_direction():
        raise ValueError('project_point_on_direction() failed')
    print("All tests passed!")
