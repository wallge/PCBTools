import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon

class LoadPoint:
    def __init__(self, net_name, coord, color):
        self.net_name = net_name
        self.coord = coord
        self.color = color

class PowerPlaneCut:
    def __init__(self, net_name, boundary, color, layer_num=None):
        self.net_name = net_name
        self.boundary = boundary
        self.color = color
        self.layer_num = layer_num


class Via:
    def __init__(self, net_name, hole_size, pad_size, coord):
        self.net_name = net_name
        self.hole_size = hole_size
        self.pad_size = pad_size
        self.coord = coord

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


def get_load_points(f):
    loads = []
    points = []
    color_lut = {}

    for line in f:
        if line == "***BOARD OUTLINE***\n":
            break
        else:

            # split the current line by semicolons and remove the newline character
            fields = line.replace('\n', '').split(';')
            if len(fields) >= 5:
                name = fields[0]
                #hole_size = float(fields[1])
                #pad_size = float(fields[2])
                coord = Point(float(fields[3]), float(fields[4]))
                color = color_lut.get(name)
                if color is None:
                    color_lut[name] = np.random.rand(3)

                load = LoadPoint(name, coord, color)
                loads.append(load)
                #X.append(coord.x)
                #Y.append(coord.y)
                points += coord.coords[:]

    return loads, points, color_lut


def get_outline(f):
    found_outline = False

    board_points = []

    for line in f:
        if line == "***BOARD OUTLINE***\n":
            found_outline = True

        if found_outline:
            # split the current line by semicolons and remove the newline character
            fields = line.replace('\n', '').split(';')
            # print( fields
            if fields[0] == "TRACK":
                board_points += [(float(fields[1]), float(fields[2])), (float(fields[3]), float(fields[4]))]

            ##elif fields[0] == "ARC":

    board_outline = Polygon(board_points).buffer(0)

    return board_outline


def parse_file(file_name):
    f = open(file_name, 'r')
    loads = get_load_points(f)
    f.seek(0)
    pcb_outline = get_outline(f)
    f.close()
    return loads, pcb_outline


def make_cells_finite(vor):
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
    ##if radius is None:
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
        #print( "p1 = ", p1, "all_ridges = ", all_ridges
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

            #far_point = find_pcb_outline_intersection(vor.vertices[v2], direction, pcb_outline)
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



def get_min_coord(_min, c):
    if _min[0] > c.x:
        _min[0] = c.x

    if _min[1] > c.y:
        _min[1] = c.y


def get_max_coord(_max, c):
    if _max[0] < c.x:
        _max[0] = c.x

    if _max[1] < c.y:
        _max[1] = c.y

def save_poly(out_file, poly_coords):
    for vertex in poly_coords:
        out_file.write(str(round(vertex[0], 2)) + ' ' + str(round(vertex[1], 2)) + ' ')
    out_file.write('\n')


def merge_voronoi_cells(loads, voronoi):
    regions, vertices = make_cells_finite(voronoi)
    merge_polys = {}

    for i, region in enumerate(regions):
        poly_vertices = vertices[region]
        polygon = Polygon(poly_vertices)
        load = loads[i]
        net_name = load.net_name

        merge_poly = merge_polys.get(net_name)
        if merge_poly is None:
            merge_polys[net_name] = polygon
        else:
            merge_polys[net_name] = merge_poly.union(polygon)

        #color = color_lut[net_name]
    return merge_polys


def create_cuts(polys_dict, pcb_outline, color_lut):
    plane_cuts = []

    for net_name, poly in polys_dict.items():
        color = color_lut[net_name]
        poly = pcb_outline.intersection(poly)

        if isinstance(poly, Polygon):
            plane_cut = PowerPlaneCut(net_name, poly, color)
            plane_cuts.append(plane_cut)
        else:
            for p in poly:
                plane_cut = PowerPlaneCut(net_name, p, color)
                plane_cuts.append(plane_cut)

    return plane_cuts


class PowerPlaneAssignment:
    def __init__(self, plane_cuts, adjacency, assignments, num_planes):
        self.plane_cuts = plane_cuts
        self.adjacency = adjacency
        self.assignments = assignments
        self.num_planes = num_planes
        #create empty dict for each of the power plane layers
        self.layer_dicts = num_planes * [{}]
        #create a  list of assigned polygons for each layer
        self.assigned = num_planes * [[]]
        #create a list of unassigned polygons for each layer
        self.unassigned = num_planes * [[]]

        for idx, plane_cut in enumerate(plane_cuts):
            #assign the layer number
            net_name = plane_cut.net_name
            layer_num = assignments[idx]
            layer_dict = self.layer_dicts[layer_num]
            cut_indices = layer_dict.get(net_name)
            if cut_indices is None:
                cut_indices = [idx]
                layer_dict[net_name] = cut_indices
            else:
                layer_dict[net_name].append(idx)

            #add this polygon to the list of vacancies on the other pcb layers
            for l in self.num_planes:
                if l == layer_num:
                    self.assigned[l].append(idx)
                else:
                    self.unassigned[l].append(idx)

    #trying to figure out how to do the assignments of the polygons
    def evaluate(self):
        #for each power layer in the board
        for plane_idx in range(self.num_planes):
            #get the indices of assigned cuts
            assigned_idxs = self.assigned[plane_idx]
            # and the indices of unassigned cuts
            unassigned_idxs = self.unassigned[plane_idx]
            while unassigned_idxs:
                for assigned_idx in assigned_idxs:
                    assigned_cut = self.plane_cuts[assigned_idx]
                    adjacent_idxs = self.adjacency[assigned_idx]

                    for adjacent_idx in adjacent_idxs:
                        adjacent_cut = self.plane_cuts[adjacent_idx]




