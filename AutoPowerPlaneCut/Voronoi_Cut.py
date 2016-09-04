
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d


class Via:
    def __init__(self, net_name, hole_size, pad_size, coord):
        self.net_name = net_name
        self.hole_size = hole_size
        self.pad_size = pad_size
        self.coord = coord



def parse_file(file_name):
    f = open(file_name, 'r')
    via_list = []

    for line in f:
        #split the current line by semicolons and remove the newline character
        fields = line.replace('\n', '').split(';')
        #print fields
        net_name = fields[0]
        hole_size = float(fields[1])
        pad_size = float(fields[2])
        coord = (float(fields[3]), float(fields[4]))
        new_via = Via(net_name, hole_size, pad_size, coord)
        #print vars(new_via)
        via_list.append(new_via)

    f.close()
    return via_list





via_list = parse_file('ItemsList.txt')

points = []
for i in via_list:
    points.append(i.coord)



# compute Voronoi tesselation
vor = Voronoi(points)

# plot
voronoi_plot_2d(vor)

# colorize
for region in vor.regions:
    if not -1 in region:
        polygon = [vor.vertices[i] for i in region]
        plt.fill(*zip(*polygon))

plt.show()






