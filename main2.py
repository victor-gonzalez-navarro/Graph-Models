import networkx as nx
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import collections
from random import randint
from scipy.special import comb


# ------------------------------------------------------------------------------------------------------------ FUNCTIONS
def plot_degree_distr(d,title):
    degreeCount = collections.Counter(d)
    deg, cnt = zip(*degreeCount.items())
    cnt = cnt/np.sum(cnt)
    plt.bar(deg, cnt, width=0.60, color='b')
    plt.title(title)
    extra = max(cnt)+max(cnt)*0.05
    plt.ylim((0, extra))
    plt.xlabel('Degree')
    plt.ylabel('Propbability')
    return extra
    #plt.show()

type = 1 # Type of graph I want to create

# --------------------------------------------------------------------------------------- Erdös-Rényi random graph model
if type == 1:
    N = 1000
    p = 0.5  # high value indicates more edges
    G = nx.Graph()
    scala = np.linspace(0, N - 1, N)
    G.add_nodes_from(scala)
    for i in range(0, N):
        for j in range(i + 1, N):
            if random.uniform(0, 1) > (1 - p):
                G.add_edge(i, j)

    # Plot degree distribution
    di = list(G.degree)
    degree = [dix[1] for dix in di]
    plt.subplot(121)
    extra = plot_degree_distr(degree,'Experimental Degree distribution')

    # Ground truth
    int_min = min(degree)
    int_max = max(degree)
    axix = np.linspace(int_min, int_max, int(int_max - int_min + 1))
    axiy = [comb(N - 1, item) * (p ** item) * ((1 - p) ** (N - 1 - item)) for item in axix]
    plt.subplot(122)
    plt.bar(axix, axiy, width=0.80, color='g')
    plt.title('Ground Truth Degree distribution')
    plt.xlabel('Degree')
    plt.ylabel('Probability')
    plt.ylim((0,extra))
    plt.suptitle('Erdös-Rényi random graph model (N='+str(N)+',p='+str(p)+')',fontsize=16)
    plt.show()

    # # Plot the graph
    # nx.draw(G, with_labels=False)
    # plt.draw()
    # plt.show()

    # Plot degree distribution networkX
    # G2 = nx.erdos_renyi_graph(N, p, seed=None, directed=False)
    # di = list(G2.degree)
    # degree = [dix[1] for dix in di]
    # plt.subplot(132)
    # plot_degree_distr(degree, 'Degree distribution NetworkX (N=' + str(N) + ',p=' + str(p) + ')')

    # # Experimental results
    # num_edges = G.number_of_edges()
    # print('The number of edges is: ' + str(num_edges))
    # # Theoretical results
    # print('The expected number of edges is: ' + str(p * N * (N - 1) / 2))
    # error = ((abs(num_edges - (p * N * (N - 1) / 2))) * 100) / (p * N * (N - 1) / 2)
    # print('The relative error over 100 is: ' + str(error))


# ----------------------------------------------------------------------------------- Watts-Strogatz “small-world” model
# There is an error on the algorithm explained in wikipedia
elif type == 2:
    N = 50
    K = 4
    beta = 0.5
    G = nx.Graph()
    scala = np.linspace(0, N - 1, N)
    G.add_nodes_from(scala)
    for i in range(0, N):
        for j in range(i + 1, N):
            value = (abs(i - j)) % (N - K / 2)
            if (value >= 0) and (value <= K / 2):
                G.add_edge(i, j)

    for ni in range(0, N):
        for nj in range(ni + 1, ni + 1 + int(K / 2)):
            if random.uniform(0, 1) < beta:
                connected = tuple(G._adj[ni]) + (ni,)
                right_node = nj % N
                G.remove_edge(ni, right_node)
                possible_connection = [e for e in list(scala) if e not in connected]
                nk = randint(0, len(possible_connection))
                G.add_edge(ni, nk)

    pos = dict()
    degree = 0
    for i in range(0, N):
        pos[i] = np.array([math.cos(degree * (math.pi / 180)), math.sin(degree * (math.pi / 180))])
        degree = degree + 360 / N

    nx.draw(G, with_labels=True, pos=pos, node_size=100)
    plt.draw()
    plt.show()


# --------------------------------------------------------------------------------------- Barabási-Albert algorithm (BA)
elif type == 3:
    N = 50
    Nini = 10
    m = 4
    G = nx.Graph()
    scala = np.linspace(0, Nini - 1, Nini)
    G.add_nodes_from(scala)
    for i in range(0, Nini):
        for j in range(i + 1, Nini):
            G.add_edge(i, j)

    # Add nodes to the network sequentially
    for i in range(Nini,N):
        G.add_node(i)
        degree = [0]
        # Recalculation of the probability distribution
        for itx in range(len(G.degree)-1):
            degree = degree + [degree[itx]+G.degree[itx],]
        sumatotal_degree = degree[-1]

        connected_nodes = []
        for midx in range(m):
            value = random.uniform(0, 1)*sumatotal_degree
            nodev = (np.digitize([value], degree)[0])-1
            if nodev not in connected_nodes:
                connected_nodes = connected_nodes + [nodev,]
                G.add_edge(nodev, i)


# --------------------------------------------------------------------------------------------- Configuration Model (CM)
elif type == 4:
    pois = int(input('Do you want a Poisson or a Power Law degree distribution? Insert a number (0-POISS, '
                    '1-POWER LAW): '))

    N = 50
    G = nx.Graph()
    scala = np.linspace(0, N - 1, N)
    G.add_nodes_from(scala)

    # Poisson Distribution
    if pois == 0:
        lambd = 5
        degree = np.random.poisson(lambd, N)
    # Powe Law Distribution
    elif pois == 1:
        mini = 1
        maxi = N-1
        gamma = 2.2
        g = 1-gamma # gamma = 1-g  -->g-1 = -gamma
        r = np.random.random(N)
        ag, bg = mini ** g, maxi ** g
        degree = np.round((ag + (bg - ag) * r) ** (1. / g))

    # Plot degree distribution
    plot_degree_distr(degree,'Degree distribution')

    # Compute acumulative degree
    degreeacum = [0]
    for itx in range(len(degree)):
        degreeacum = degreeacum + [degreeacum[itx] + degree[itx],]
    sumatotal_degree = degreeacum[-1]

    for nodeval in range(N):
        tries = 0
        while (degree[nodeval] != 0) and (tries <100):
            value = random.uniform(0, 1)*sumatotal_degree
            rv = (np.digitize([value], degreeacum)[0])-1
            if (rv != nodeval) and (not(G.has_edge(rv, nodeval))):
                G.add_edge(nodeval, rv)
                degree[nodeval] = degree[nodeval] -1
                degree[rv] = degree[rv] -1

                # Compute acumulative degree
                degreeacum = [0]
                for itx in range(len(degree)):
                    degreeacum = degreeacum + [degreeacum[itx] + degree[itx], ]
                sumatotal_degree = degreeacum[-1]
            tries = tries + 1
