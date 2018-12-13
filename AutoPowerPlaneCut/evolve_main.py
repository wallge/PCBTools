import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon
import itertools
import voronoi_cut
import random

from deap import base
from deap import creator
from deap import tools

results, pcb_outline = voronoi_cut.parse_file('ItemsList1.txt')
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


creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()
# Attribute generator
toolbox.register("attr_bool", random.randint, 0, 1)
# Structure initializers
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, num_nets)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
# CXPB  is the probability with which two individuals
#       are crossed
#
# MUTPB is the probability for mutating an individual
CXPB, MUTPB = 0.5, 0.2

def eval(individual):
    points_by_layer = [[], []]
    loads_by_layer = [[], []]

    for idx, c in enumerate(individual):
        net_name = net_names[idx]
        layer_points = points_dict[net_name]
        points_by_layer[c] += layer_points
        loads_by_layer[c] += loads_dict[net_name]

    if 0 in individual:
        # compute Voronoi tesselation
        vor = Voronoi(points_by_layer[0])
        merged_polys = voronoi_cut.merge_voronoi_cells(loads_by_layer[0], vor)
        plane_cuts0 = voronoi_cut.create_cuts(merged_polys, pcb_outline, color_lut)
    else:
        plane_cuts0 = [None]

    if 1 in individual:
        # compute Voronoi tesselation
        vor = Voronoi(points_by_layer[1])
        merged_polys = voronoi_cut.merge_voronoi_cells(loads_by_layer[1], vor)
        plane_cuts1 = voronoi_cut.create_cuts(merged_polys, pcb_outline, color_lut)
    else:
        plane_cuts1 = [None]

    poly_count = len(plane_cuts0) + len(plane_cuts1)

    return poly_count, ##1.0 / poly_count,


toolbox.register("evaluate", eval)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)

random.seed(64)
pop = toolbox.population(n=200)
print("Start of evolution")
fitnesses = list(map(toolbox.evaluate, pop))
for ind, fit in zip(pop, fitnesses):
    ind.fitness.values = fit

print("  Evaluated %i individuals" % len(pop))

# Extracting all the fitnesses of
fits = [ind.fitness.values[0] for ind in pop]

# Variable keeping track of the number of generations
g = 0

# Begin the evolution
while g < 20:
    # A new generation
    g = g + 1
    print("-- Generation %i --" % g)

    # Select the next generation individuals
    offspring = toolbox.select(pop, len(pop))
    # Clone the selected individuals
    offspring = list(map(toolbox.clone, offspring))

    # Apply crossover and mutation on the offspring
    for child1, child2 in zip(offspring[::2], offspring[1::2]):

        # cross two individuals with probability CXPB
        if random.random() < CXPB:
            toolbox.mate(child1, child2)

            # fitness values of the children
            # must be recalculated later
            del child1.fitness.values
            del child2.fitness.values

    for mutant in offspring:

        # mutate an individual with probability MUTPB
        if random.random() < MUTPB:
            toolbox.mutate(mutant)
            del mutant.fitness.values

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
    fitnesses = map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    print("  Evaluated %i individuals" % len(invalid_ind))

    # The population is entirely replaced by the offspring
    pop[:] = offspring

    # Gather all the fitnesses in one list and print the stats
    fits = [ind.fitness.values[0] for ind in pop]

    length = len(pop)
    mean = sum(fits) / length
    sum2 = sum(x * x for x in fits)
    std = abs(sum2 / length - mean ** 2) ** 0.5

    print("  Min %s" % min(fits))
    print("  Max %s" % max(fits))
    print("  Avg %s" % mean)
    print("  Std %s" % std)

print("-- End of (successful) evolution --")

best_ind = tools.selBest(pop, 1)[0]
print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))


points_by_layer = [[], []]
loads_by_layer = [[], []]
for idx, c in enumerate(best_ind):
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
    center = cut.boundary.representative_point()
    plt.annotate(
        cut.net_name,
        xy=(center.x, center.y), xytext=(-20, 20),
        textcoords='offset points', ha='right', va='bottom',
        #bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    X = []
    Y = []
    last_item = None
    for load in loads_by_layer[0]:

        X.append(load.coord.x)
        Y.append(load.coord.y)

    plt.plot(X, Y, 'ko')

plt.subplot(1, 2, 2)
for cut in plane_cuts1:
    poly_vertices = list(cut.boundary.exterior.coords)
    #voronoi_cut.save_poly(output, poly_vertices)
    plt.fill(*zip(*poly_vertices), c=cut.color, alpha=0.4, )
    center = cut.boundary.representative_point()
    plt.annotate(
        cut.net_name,
        xy=(center.x, center.y), xytext=(-20, 20),
        textcoords='offset points', ha='right', va='bottom',
        # bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    X = []
    Y = []
    for load in loads_by_layer[1]:
        X.append(load.coord.x)
        Y.append(load.coord.y)
    plt.plot(X, Y, 'ko')

plt.show()
