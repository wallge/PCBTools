import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon

import voronoi_cut


loads, pcb_outline = voronoi_cut.parse_file('ItemsList0.txt')
loads, points, color_lut = loads
# compute Voronoi tesselation
vor = Voronoi(points)
merged_polys = voronoi_cut.merge_voronoi_cells(loads, vor)
plane_cuts = voronoi_cut.create_cuts(merged_polys, pcb_outline, color_lut)

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



