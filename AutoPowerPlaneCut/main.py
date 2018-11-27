import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon

import voronoi_cut


loads, pcb_outline = voronoi_cut.parse_file('ItemsList0.txt')
loads, points, color_lut = loads
# compute Voronoi tesselation
vor = Voronoi(points)

regions, vertices = voronoi_cut.make_cells_finite(vor)
merge_polys = {}

for i, region in enumerate(regions):
    poly_vertices = vertices[region]
    polygon = Polygon(poly_vertices)
    load = loads[i]
    net_name = load.name

    merge_poly = merge_polys.get(net_name)
    if merge_poly is None:
        merge_polys[net_name] = polygon
    else:
        merge_polys[net_name] = merge_poly.union(polygon)

    color = color_lut[net_name]

output = open('Vertices.txt', 'w')

for key, val in merge_polys.items():
    color = color_lut[key]
    val = pcb_outline.intersection(val)
    if isinstance(val, Polygon):
        poly_vertices = list(val.exterior.coords)
        voronoi_cut.save_poly(output, poly_vertices)
        plt.fill(*zip(*poly_vertices), c=color, alpha=0.4, )
    else:
        for poly in val:
            poly_vertices = list(poly.exterior.coords)
            voronoi_cut.save_poly(output, poly_vertices)
            plt.fill(*zip(*poly_vertices), c=color, alpha=0.4, )

output.close()

X = []
Y = []
for load in loads:
    X.append(load.coord.x)
    Y.append(load.coord.y)

plt.plot(X, Y, 'ko')

plt.show()



