import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon
import itertools
import voronoi_cut

results, pcb_outline = voronoi_cut.parse_file('ItemsList.txt')
loads, points, color_lut = results

points_dict = {}
loads_dict = {}
for load in loads:
    net_name = load.net_name
    items = points_dict.get(net_name)
    if items is None:
        points_dict[net_name] = load.coord.coords[:]
        loads_dict[net_name] = [load]
    else:
        points_dict[net_name] += load.coord.coords[:]
        loads_dict[net_name].append(load)


net_names = list(points_dict.keys())
num_nets = len(net_names)

configs = [list(i) for i in itertools.product([0, 1], repeat=num_nets)]
config_poly_counts = []

min_poly_count = 1e9
min_config = None

for ci, config in enumerate(configs):
    points_by_layer = [[], []]
    loads_by_layer = [[], []]

    for idx, c in enumerate(config):
        net_name = net_names[idx]
        layer_points = points_dict[net_name]
        points_by_layer[c] += layer_points
        loads_by_layer[c] += loads_dict[net_name]

    if 0 in config:
        # compute Voronoi tesselation
        vor = Voronoi(points_by_layer[0])
        merged_polys = voronoi_cut.merge_voronoi_cells(loads_by_layer[0], vor)
        plane_cuts0 = voronoi_cut.create_cuts(merged_polys, pcb_outline, color_lut)
    else:
        plane_cuts0 = [None]

    if 1 in config:
        # compute Voronoi tesselation
        vor = Voronoi(points_by_layer[1])
        merged_polys = voronoi_cut.merge_voronoi_cells(loads_by_layer[1], vor)
        plane_cuts1 = voronoi_cut.create_cuts(merged_polys, pcb_outline, color_lut)
    else:
        plane_cuts1 = [None]

    poly_count = len(plane_cuts0) + len(plane_cuts1)

    if poly_count <= min_poly_count:
        min_poly_count = poly_count
        min_config = config
        #print(poly_count, config)
        net_count_dict = {}
        if 0 in config:
            for cut in plane_cuts0:
                net_count = net_count_dict.get(cut.net_name)
                if net_count is None:
                    net_count_dict[cut.net_name] = 1
                else:
                    net_count_dict[cut.net_name] += 1

        if 1 in config:
            for cut in plane_cuts1:
                net_count = net_count_dict.get(cut.net_name)
                if net_count is None:
                    net_count_dict[cut.net_name] = 1
                else:
                    net_count_dict[cut.net_name] += 1

        #for key, val in net_count_dict.items():
        #    if val > 1:
        #        print(key, val)


points_by_layer = [[], []]
loads_by_layer = [[], []]
for idx, c in enumerate(min_config):
    net_name = net_names[idx]
    layer_points = points_dict[net_name]
    points_by_layer[c] += layer_points
    loads_by_layer[c] += loads_dict[net_name]

vor0 = Voronoi(points_by_layer[0])
merged_polys0 = voronoi_cut.merge_voronoi_cells(loads_by_layer[0], vor0)
plane_cuts0 = voronoi_cut.create_cuts(merged_polys0, pcb_outline, color_lut)

vor1 = Voronoi(points_by_layer[1])
merged_polys1 = voronoi_cut.merge_voronoi_cells(loads_by_layer[1], vor1)
plane_cuts1 = voronoi_cut.create_cuts(merged_polys1, pcb_outline, color_lut)
plt.subplot(1, 2, 1)

for cut in plane_cuts0:
    poly_vertices = list(cut.boundary.exterior.coords)
    #voronoi_cut.save_poly(output, poly_vertices)
    plt.fill(*zip(*poly_vertices), c=cut.color, alpha=0.4, )
    X = []
    Y = []
    for load in loads_by_layer[0]:
        X.append(load.coord.x)
        Y.append(load.coord.y)

    plt.plot(X, Y, 'ko')

plt.subplot(1, 2, 2)
for cut in plane_cuts1:
    poly_vertices = list(cut.boundary.exterior.coords)
    #voronoi_cut.save_poly(output, poly_vertices)
    plt.fill(*zip(*poly_vertices), c=cut.color, alpha=0.4, )
    X = []
    Y = []
    for load in loads_by_layer[1]:
        X.append(load.coord.x)
        Y.append(load.coord.y)
    plt.plot(X, Y, 'ko')

    '''
    if isinstance(pcb_outline, Polygon):
        X, Y = pcb_outline.exterior.coords.xy
    else:
        X = []
        Y = []
        for p in pcb_outline:
            x, y = p.exterior.coords.xy
            X += x
            Y += y

    plt.plot(X, Y, 'ko')
    '''

plt.show()
#9 [0, 0, 0, 1, 1, 0, 1, 1]

'''
output = open('Vertices.txt', 'w')
for cut in plane_cuts:
    poly_vertices = list(cut.boundary.exterior.coords)
    voronoi_cut.save_poly(output, poly_vertices)
    plt.fill(*zip(*poly_vertices), c=cut.color, alpha=0.4, )

output.close()

X = []
Y = []
for load in loads:
    X.append(load.coord.x)
    Y.append(load.coord.y)

plt.plot(X, Y, 'ko')

plt.show()
'''


