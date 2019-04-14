import networkx as nx
import random
import numpy as np

from random import randint

# ------------------------------------------------------------------------------------------------------ HYPERPARAMETERS
type = 1  # Type of graph I want to create (Ptions: 1, 2, 3, 4)

N = 100

# ******************************************* Erdös-Rényi random graph model
p = 0.2  # high value indicates more edges
# ******************************************* Watts-Strogatz “small-world” model
K = 8
beta = 0.5  # Probability of rewiring
# ******************************************* Barabási-Albert algorithm (BA)
Nini = 10  # number of initial nodes
m = 3  # m defines the average degree
# ******************************************* Configuration Model (CM)
pois = 0  # 0=Poisson, 1=Power Law
lambd = 2  # Parameter of Possion
gamma = 2.7  # Parameter of Power Law


# --------------------------------------------------------------------------------------- Erdös-Rényi random graph model
if type == 1:
    G = nx.Graph()
    scala = np.linspace(0, N - 1, N)
    G.add_nodes_from(scala)
    for i in range(0, N):
        for j in range(i + 1, N):
            if random.uniform(0, 1) > (1 - p):
                G.add_edge(i, j)


# ----------------------------------------------------------------------------------- Watts-Strogatz “small-world” model
if type == 2:
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

# --------------------------------------------------------------------------------------- Barabási-Albert algorithm (BA)
if type == 3:
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
if type == 4:
    G = nx.Graph()
    scala = np.linspace(0, N - 1, N)
    G.add_nodes_from(scala)

    # Generate samles from Poisson Distribution
    if pois == 0:
        degree = []
        for ite in range(N):
            newval = np.random.poisson(lambd, 1)[0]
            if newval >= 1:
                degree = degree + [newval,]
        degree = np.array(degree)
        nbins = 6
    # Generate samples from Power Law Distribution
    elif pois == 1:
        x = np.linspace(start=1, stop=5000, num=5000)
        C = 100
        y = C*(x**(-gamma))
        degreeacum = [0]
        for itx in range(len(y)):
            degreeacum = degreeacum + [degreeacum[itx] + y[itx], ]
        sumatotal_degree = degreeacum[-1]
        degree = []
        for nodeval in range(N):
            value = random.uniform(0, 1) * sumatotal_degree
            degree = degree + [np.digitize([value], degreeacum)[0], ]

    # Compute acumulative degree
    degreeacum = [0]
    for itx in range(len(degree)):
        degreeacum = degreeacum + [degreeacum[itx] + degree[itx],]
    sumatotal_degree = degreeacum[-1]
    for nodeval in range(N):
        tries = 0
        while (degree[nodeval] != 0) and (tries <100000):
            value = random.uniform(0, 1)*sumatotal_degree
            rv = (np.digitize([value], degreeacum)[0])-1
            if (rv != nodeval) and (not(G.has_edge(rv, nodeval))):
                G.add_edge(nodeval, rv)
                degree[nodeval] = degree[nodeval] -1
                degree[rv] = degree[rv]-1
                # Compute acumulative degree
                degreeacum = [0]
                for itx in range(len(degree)):
                    degreeacum = degreeacum + [degreeacum[itx] + degree[itx], ]
                sumatotal_degree = degreeacum[-1]
            tries = tries + 1
    # Remove nodes with degree 0
    remove = [key for key,value in dict(G.degree).items() if value == 0]
    G.remove_nodes_from(remove)
