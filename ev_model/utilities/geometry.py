import numpy as np


def two_points_distance(p1, p2):
    """
    Computer the Euclidean distance between two points in 2D or 3D
    Geometric interpretation: its the magnitude of the vector p1p2
    :param p1:
    :param p2:
    :return: distance between the points
    """
    p1p2 = p2 - p1
    return np.sqrt(np.dot(p1p2, p1p2))

def vector_orthogonals_2d(vector):
    dx = vector[0]
    dy = vector[1]
    return np.array([[dy, -dx], [-dy, dx]])

def two_points_normals_2d(p1, p2):
    """
    Takes the direction vector from p1 to p2 and applies a 90 degree rotation.

    From https://stackoverflow.com/a/1252738/3830240:
    The matrix representation of the general 2D transformation looks like this:

    x' = x cos(t) - y sin(t)
    y' = x sin(t) + y cos(t)
    where (x,y) are the components of the original vector and (x', y') are the transformed components.

    If t = 90 degrees, then cos(90) = 0 and sin(90) = 1

    :param p1:
    :param p2:
    :return: The normal vectors to p1p2 [[left],[right]]
    """
    p1p2 = p2 - p1

    magnitude = two_points_distance(p1, p2)
    dx = p1p2[0] / magnitude
    dy = p1p2[1] / magnitude
    return np.array([[dy, -dx], [-dy, dx]])

def vector_magnitude(vector):
    return np.sqrt(np.dot(vector, vector))

def unit_vector_2points(vector):
    """
    Normalises the given vector
    :param vector: The n-dimensional vector which needs to be normalised
    :return: the unit vector corresponding to the given vector
    """
    magnitude = vector_magnitude(vector)
    if magnitude == 0:
        print('zero vector')
        return 0, 0
    return vector/magnitude

def middle_point(p1, p2):
    """
    Computes the midpoint laying between two points. Both points must have the same dimension
    :param p1: the first set of coordinates
    :param p2: the second set of coordinates
    :return: the computed coordinates of the point laying at half the distance between p1 and p2
    """
    return (p1 + p2)/2

def reflect_vector(incidence_vector, wall_normal):
    """
    Computes the reflected vector from a wall
    Modified from https://gamedev.stackexchange.com/a/23676
    public static Vector2 Reflect(Vector2 vector, Vector2 normal)
    {
        return vector - 2 * Vector2.Dot(vector, normal) * normal;
    }
    :param incidence_vector:
    :param wall_normal:
    :return:
    """
    return incidence_vector - 2 * np.dot(incidence_vector, wall_normal) * wall_normal

def reflect_vector_2pts(p0, p1, wall_normal):
    return reflect_vector(p1 - p0, wall_normal)

def make_3d(vec):
    tmp = np.zeros(3)
    tmp[:2] = vec
    return tmp

def project_points_on_direction_from_matrix(source, direction, dist):
    return source + direction * dist

def project_point_on_direction(source, direction, dist=1.0, debug_log=False):
    """
    Projects a point from the given source in the given direction and distance.
    :param source: coordinates to use as the starting point for the displacement
    :param direction: direction vector to use for the projection
    :param dist: the distance to displace from the source point
    :return: the computed coordinates of the destination
    """
    if (source.shape == (1, 2) or source.shape == (1, 3)) and (direction.shape == (1, 2) or direction.shape == (1, 3)):
        print('projection delegated to matrix supporting fnction')
        return project_points_on_direction_from_matrix(source, direction, dist)

    # print('s.sh:', source.shape, 'd.sh:', direction.shape)
    if source.shape != (2,) and source.shape != (3,):
        raise ValueError('Source has non-supported shape')
    if direction.shape != (2,) and direction.shape != (3,):
        raise ValueError('Direction has non-supported shape')
    if source.ndim > 2 or direction.ndim > 2:
        raise ValueError('Non supported dimensions')

    # different dimensions but supported shapes, match the dimensions
    if(debug_log):
        print(f'source.ndim: {source.ndim}, direction.ndim: {direction.ndim}')
    if source.ndim != direction.ndim:
        if source.ndim > direction.ndim:
            direction = np.array([direction])
        else:
            source = np.array([source])
        if(debug_log):
            print(f'new dims >> source.ndim: {source.ndim}, direction.ndim: {direction.ndim}')

    if(debug_log):
        print(f'source.shape: {source.shape}, direction.shape: {direction.shape}')
    # same dimensions, check shape
    if source.shape != direction.shape:
        # non matching shapes need to be matched
        # we handle vectors and matrices accordingly
        if source.ndim == 1:
            # 2d vectors must become 3d with 0 valued z-axis
            if len(source) == 2:
                source = make_3d(source)
            else:
                # resize direction
                direction = make_3d(direction)
        else:
            # we are dealing with (1 by n) matrices
            if source.shape[1] == 2:
                # 2d vectors must become 3d with 0 valued z-axis
                source[0] = make_3d(source[0])
            else:
                direction[0] = make_3d(direction[0])
        if(debug_log):
            print(f'NEW SHAPES source.shape: {source.shape}, direction.shape: {direction.shape}')

    # check the direction vector is a unit norm vector (l2-norm)
    magnitude = vector_magnitude(direction)
    if magnitude != 1.0:
        direction = unit_vector_2points(direction)
        magnitude = vector_magnitude(direction)
    if(debug_log):
        print('source:',source, 'direction:',direction, 'magnitude:',magnitude,'dist',dist)
        print('this values produce:',source + (direction * dist))

    return source + (direction * dist)

class Plane:
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        # self.p = d / [a, b, c]

    @classmethod
    def from_normal_and_point(cls, normal, point):
        print(f"plane constructed from p={point} and n={normal} ", end="")
        if normal.shape[0] == 2:
            tmp = np.zeros(3)
            tmp[:2] = normal
            normal = tmp
        norm = np.sqrt(np.dot(normal, normal))
        if norm >= 1.:
            unit_norm = normal / norm
            print(" normalised")
        else:
            unit_norm = normal
            print(" non-normalised")
        a, b, c = unit_norm
        d = -np.dot(unit_norm, point)
        return cls(a, b, c, d)

    @classmethod
    def from_three_points(cls, p1, p2, p3):
        print("Plane constructed from 3 points ", end="")
        v0 = p1 - p2
        v1 = p3 - p2
        norm = np.cross(v1, v0)
        magnitude = np.sum(norm)
        if magnitude > 1.0:
            unit_norm = norm / magnitude
            print("normalised")
        else:
            unit_norm = norm
            print("non-normalised")
        a, b, c = unit_norm
        d = -np.dot(unit_norm, p1)
        return cls(a, b, c, d)

    def distance(self, point):
        return np.dot([self.a, self.b, self.c], point) + self.d

    def __str__(self):
        return f"Ax={self.a:1.3f} By={self.b:1.3f} Cz={self.c:1.3f} D={self.d:1.3f}"

    def get_norm(self):
        return np.array([self.a, self.b, self.c])

class Segment(object):
    def __init__(self, x1, y1, x2, y2, mx, my, internal=False):
        self.p1 = np.array([x1, y1])
        self.p2 = np.array([x2, y2])
        self.mp = np.array([mx, my])
        self.distance_p1p2 = two_points_distance(self.p1, self.p2)
        self.internal = internal
        # The direction we compute here is the p1p2 vector
        self.direction = self.p2 - self.p1
        self.direction_length = vector_magnitude(self.direction)
        self.direction_unit = self.direction / self.direction_length
        self.orthogonals = vector_orthogonals_2d(self.direction)
        self.orthonormals = self.orthogonals / self.direction_length
 
        self.side = 1
        self._normal = self.orthogonals[0 if self.internal else 1]
        self._unormal = self.orthonormals[0 if self.internal else 1]
        self._binormal = self.orthogonals[1 if self.internal else 0]
        self._ubinormal = self.orthonormals[1 if self.internal else 0]
    
    @property
    def location(self):
        return self.mp
    
    @property
    def normal(self):
        """
        The normals of a CCW defined shape 'point outside' the shape.
        :return:
        the normal if the shape is within another
        the binormal otherwise
        """
        return self._normal

    @property
    def binormal(self):
        return self._binormal

    @property
    def unit_normal(self):
        return self._unormal

    @property
    def unit_binormal(self):
        return self._ubinormal

    @property
    def normal_length(self):
        return vector_magnitude(self._normal)
