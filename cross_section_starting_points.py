from scipy import misc
from skimage import measure
from ev_model.utilities import geometry as geo
import matplotlib.pyplot as plt
import cmath


def process_image(img, min_section_vertices=210):
    foreground = (img[:, :, 3] != 0) & (img[:, :, 1] > 240)
    plt.imsave('01_foreground.png', img)

    # identify the image contours
    contours = measure.find_contours(foreground, 0.8)

    # identify those regions with more than n vertices
    greater_than_n = []

    # draw the identified contours with a distinctive colour
    plt.figure()
    fig, (ax) = plt.subplots(1, 1, sharex='col', figsize=(20, 10))
    ax.set_title('Identified contours w/more than ' + str(min_section_vertices) + ' edges highlighted in colour')
    ax.imshow(img, interpolation='nearest', cmap=plt.cm.get_cmap('grey'))
    for n, contour in enumerate(contours):
        if len(contour) > min_section_vertices:
            greater_than_n.append(n)
            ax.plot(contour[:, 1], contour[:, 0], linewidth=2)
    plt.savefig('02_contours.png')
    return [contours, greater_than_n]


def special_contours(contours, greater_than_n, external, internal):
    plt.figure()
    plt.plot(contours[greater_than_n[external]][:, 1], contours[greater_than_n[external]][:, 0])
    for i in internal:
        plt.plot(contours[greater_than_n[i]][:, 1], contours[greater_than_n[i]][:, 0])
    plt.axis([0, 2588, 1858, 0])
    plt.savefig('03_special_contours.png')


def find_vertices(contours, greater_than_n, simplify=True, tolerance=0.25):
    """

    :param contours:
    :param greater_than_n:
    :param simplify:
    :param tolerance:
    :return: a list containing a list of vertices for each sub section / region in the cross section
    """
    vertices = []

    # manually store the first polygons identified
    for i in range(len(greater_than_n)):
        if simplify:
            #  ()
            section_vertices = measure.approximate_polygon(contours[greater_than_n[i]], tolerance)
        else:
            section_vertices = contours[greater_than_n[i]]

        # add this to the global vertices storage
        vertices.append(section_vertices)
    return vertices


def identify_source_points(vertices, use_internal_normal=True):
    source_points = []
    # traverse the vertices on the contour
    n = len(vertices)
    for i in range(n-1):
        p1 = (vertices[i][0], vertices[i][1])
        p2 = (vertices[i+1][0], vertices[i+1][1])

        source_points.append({'p': geo.middle_point(p1, p2), 'normal': geo.unit_vector_2p(
            geo.two_points_normal(p1, p2, use_internal_normal))})
    return source_points


def cells_in_section(section_vertices, cell_width_mean, cell_width_variance=0.0001,
                     oriented_towards_centre=True, secretory_cell_probability=0.5):
    """
    Traverse each pair of vertices forming the section.
    Get the section orientation using its normal
    Compute the number of cells fitting between both vertices.

    :param section_vertices:
    :param cell_width_mean: in millimetres
    :param cell_width_variance: in millimetres
    :param oriented_towards_centre:
    :param secretory_cell_probability: Probability for any cells of being a secretory cell. Usually 0.5
    :return: list of cells in the current section
    """
    cells = list()
    limit = len(section_vertices)
    for i in range(1, limit):
        distance_between_vertices = geo.distance(section_vertices[i-1], section_vertices[i])
        # compute the number of cells to generate
        n_cells = distance_between_vertices / cell_width_mean
    # then check the segment between the last and the first vertices

    return cells


def generate_cells(vertices, invert_orientation, secretory_cell_probability=0.5, cell_diameter_mm2=0.005):
    """
    A typical pig oviduct measures 150mm in length with a diameter of 1 mm.
    Using the formula for the lateral surface area of a cylinder, we get a total area of 706.96 mm^2.
    LSA = 2 * pi * r * h
    (We are not considering the internal folds here)
    Considering a cell size of 5um^2, there would be a total of 141,372,000 cells.

    :param vertices: list with lists of vertices forming each shape of the cross section. One list per shape.
    :param invert_orientation: list of section ids which cells must be oriented facing away from the centre of the shape
    :param secretory_cell_probability: Probability for any cells of being a secretory cell. Usually 0.5
    :param cell_diameter_mm2: target cell diameter in millimetres squared
    :return:
    """
    cell_width = cmath.sqrt(cell_diameter_mm2/cmath.pi) * 2
    coordinates = []
    # As previously stated, contours 12 & 14 lay within contour 11, their normals must point in the other direction
    for i in range(len(vertices)):
        # previously, we were identifying the source points directly
        # now, we need to compute how many
        if i in invert_orientation:
            source_points = identify_source_points(vertices[i], False)
            cells = cells_in_section(vertices[i], cell_width, oriented_towards_centre=False,
                                     secretory_cell_probability=secretory_cell_probability)
        else:
            source_points = identify_source_points(vertices[i])
            cells = cells_in_section(vertices[i], cell_width, secretory_cell_probability=secretory_cell_probability)
        coordinates.append(source_points)

    for i in range(len(coordinates)):
        print('figure:', i, 'n.s.coordinates:', len(coordinates[i]))
    return coordinates


if __name__ == '__main__':
    image = misc.imread('./ampulla_(3-1to3-5_L2-3_tr-amp)_2.png')
    [cs, gtn] = process_image(image)
    must_flip_orientation = [4, 5]
    special_contours(cs, gtn, 3, internal=must_flip_orientation)

    all_vertices = find_vertices(cs, gtn, simplify=True)
    starting_coordinates = generate_cells(all_vertices, invert_orientation=must_flip_orientation)

    # export the coordinates
