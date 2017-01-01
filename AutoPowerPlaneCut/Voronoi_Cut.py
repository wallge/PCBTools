
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
##import numpy as np
##import matplotlib.pyplot as plt
##from scipy.spatial import Voronoi

class Via:
    def __init__(self, net_name, hole_size, pad_size, coord):
        self.net_name = net_name
        self.hole_size = hole_size
        self.pad_size = pad_size
        self.coord = coord


class Coord:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Track:
    def __init__(self, coord0, coord1, width):
        self.coord0 = coord0
        self.coord1 = coord1
        self.width = width

class Arc:
    def __init__(self, center, radius, angle1, angle2, width):
        self.center = center
        self.radius = radius
        self.angle1 = angle1
        self.width = width


def get_vias(f):
    via_list = []
    line_number = 0

    for line in f:
        ##print line_number
        if line == "***BOARD OUTLINE***\n":
            break
        else:

            # split the current line by semicolons and remove the newline character
            fields = line.replace('\n', '').split(';')
            # print fields

            if len(fields) >= 5:
                net_name = fields[0]
                hole_size = float(fields[1])
                pad_size = float(fields[2])
                coord = (float(fields[3]), float(fields[4]))
                new_via = Via(net_name, hole_size, pad_size, coord)
                # print vars(new_via)
                via_list.append(new_via)

        line_number += 1

    return via_list


def get_outline(f):
    board_outline = []
    found_outline = False

    for line in f:
        if line == "***BOARD OUTLINE***\n":
            found_outline = True

        if found_outline:
            # split the current line by semicolons and remove the newline character
            fields = line.replace('\n', '').split(';')
            # print fields
            if fields[0] == "TRACK":
                new_track = Track(Coord(float(fields[1]), float(fields[2])), Coord(float(fields[3]), float(fields[4])), 0)
                board_outline.append(new_track)

            ##elif fields[0] == "ARC":


    return board_outline


def parse_file(file_name):
    f = open(file_name, 'r')
    via_list = get_vias(f)
    f.seek(0)
    pcb_outline = get_outline(f)

    f.close()
    return via_list, pcb_outline





via_list, pcb_outline = parse_file('ItemsList.txt')

points = []
X = []
Y = []
for i in via_list:
    #print type(i.coord), type(i.coord[0])
    points.append([i.coord[0], i.coord[1]])
    X.append(i.coord[0])
    Y.append(i.coord[1])




##if False:
##    # compute Voronoi tesselation
##    vor = Voronoi(points)
##
##    # plot
##    voronoi_plot_2d(vor)
##
##    # colorize
##    for region in vor.regions:
##        if not -1 in region:
##            polygon = [vor.vertices[i] for i in region]
##            plt.fill(*zip(*polygon))
##
##    plt.show()


def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)

# make up data points
#np.random.seed(1234)
#points = np.random.rand(15, 2)

# compute Voronoi tesselation
vor = Voronoi(points)
#print vor.points
#print vor.vertices

#output = open('Vertices.txt', 'w')
#for region in vor.regions:
#    for i in range(0, len(region)-1):
#        point0 = vor.vertices[i]
#        point1 = vor.vertices[i+1]
#        #print i, point0, point1
#        output.write(str(point0[0]) + ' ' + str(point0[1]) + ' ' + str(point1[0]) + ' ' + str(point1[1]) + '\n')
#
#output.close()





# plot
regions, vertices = voronoi_finite_polygons_2d(vor)
#print "--"
#print regions
#print "--"
#print len(vertices)
#print "--"
#print vertices



#output = open('Vertices.txt', 'w')
#for region in vor.regions:
#    for i in range(0, len(region)-1):
#        point0 = vor.vertices[i]
#        point1 = vor.vertices[i+1]
#        #print i, point0, point1
#        output.write(str(point0[0]) + ' ' + str(point0[1]) + ' ' + str(point1[0]) + ' ' + str(point1[1]) + '\n')
#
#output.close()

#"""
output = open('Vertices.txt', 'w')
# colorize
for region in regions:
    polygon = vertices[region]
    print 'polygon:', polygon
    for vertex in polygon:
        output.write(str(round(vertex[0], 2)) + ' ' + str(round(vertex[1], 2)) + ' ')
    output.write('\n')

    plt.fill(*zip(*polygon), alpha=0.4)

output.close()

#"""
#print type(points), type(points[0]), type(points[0][0])
#print points[:][0]
#print points
#plt.plot(points, 'ko')
#X = points[:]
#print X
#print points
#Y = points[1][:]
#print len(X), len(Y)
#print points[:][1]


min_coord = Coord(100000, 100000)
max_coord = Coord(-100000, -100000)

def get_min_coord(min, c):
    if min.x > c.x:
        min.x = c.x

    if min.y > c.y:
        min.y = c.y

def get_max_coord(max, c):
    if max.x < c.x:
        max.x = c.x

    if max.y < c.y:
        max.y = c.y


for item in pcb_outline:
    if item.__class__.__name__ == "Track":
        plt.plot((item.coord0.x, item.coord1.x), (item.coord0.y, item.coord1.y), 'r-')
        get_min_coord(min_coord, item.coord0)
        get_min_coord(min_coord, item.coord1)
        get_max_coord(max_coord, item.coord0)
        get_max_coord(max_coord, item.coord1)

plt.plot(X, Y, 'ko')
##plt.xlim(vor.min_bound[0] - 0.1, vor.max_bound[0] + 0.1)
##plt.ylim(vor.min_bound[1] - 0.1, vor.max_bound[1] + 0.1)


width = max_coord.x - min_coord.x
height = max_coord.y - min_coord.y

plt.xlim(min_coord.x - 0.1*width, max_coord.x + 0.1*width)
plt.ylim(min_coord.y - 0.1*height, max_coord.y + 0.1*height)

#print min_coord.x, min_coord.y
#print max_coord.x, max_coord.y

plt.show()



