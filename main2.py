import networkx as nx
import random

from random import randint
from scipy.special import comb
from auxiliary_functions import *


# ------------------------------------------------------------------------------------------------------------ MAIN CODE
type = 4  # Type of graph I want to create

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


# ----------------------------------------------------------------------------------- Watts-Strogatz “small-world” model
# There is an error on the algorithm explained in wikipedia
elif type == 2:
    N = 1000
    K = 8
    beta = 0.2
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

    # Plot degree distribution
    di = list(G.degree)
    degree = [dix[1] for dix in di]
    plt.subplot(121)
    extra = plot_degree_distr(degree,'Experimental Degree distribution')

    # Ground truth
    int_min = min(degree)
    int_max = max(degree)
    axix = np.linspace(int_min, int_max, int(int_max - int_min + 1))
    axiy = []
    for item in axix:
        suma = 0
        for it in range(0,(int(min((item-K/2),(K/2)))+1)):
            suma = suma + ((comb(K/2, it)*((1-beta)**(it))*(beta**((K/2)-it))*((beta*(K/2))**(item-(K/2)-it))*(np.exp(
                -beta*K/2)))/(math.factorial(item-(K/2)-it)))
        axiy = axiy + [suma,]
    plt.subplot(122)
    plt.bar(axix, axiy, width=0.60, color='g')
    plt.title('Ground Truth Degree distribution')
    plt.xlabel('Degree k')
    plt.ylabel('Fraction of vertices Pk having degree k')
    plt.ylim((0,extra))
    plt.suptitle('Watts-Strogatz graph model (N='+str(N)+',k='+str(K)+',p='+str(beta)+')',fontsize=16)
    plt.show()

    # # Plot the graph
    # pos = dict()
    # degree = 0
    # for i in range(0, N):
    #     pos[i] = np.array([math.cos(degree * (math.pi / 180)), math.sin(degree * (math.pi / 180))])
    #     degree = degree + 360 / N
    # nx.draw(G, with_labels=True, pos=pos, node_size=100)
    # plt.draw()
    # plt.show()

# --------------------------------------------------------------------------------------- Barabási-Albert algorithm (BA)
elif type == 3:
    N = 1000
    Nini = 10
    m = 4  # m defines the average degree
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

    # Plot degree distribution Linear Scale
    di = list(G.degree)
    degree = [dix[1] for dix in di]
    plt.subplot(131)
    extra = plot_degree_distr(degree, 'Experimental Degree distribution')
    # Ground truth Linear Scale
    int_min = min(degree)
    int_max = max(degree)
    axix = np.linspace(int_min, int_max, int(int_max - int_min + 1))
    axiy = [2*m*(m+1)/(item*(item+1)*(item+2)) for item in axix]
    plt.subplot(133)
    plt.bar(axix, axiy, width=0.80, color='g')
    plt.title('Ground Truth Degree distribution')
    plt.xlabel('Degree k')
    plt.ylabel('Fraction of vertices Pk having degree k')
    plt.ylim((0, extra))
    plt.suptitle('Barabási-Albert graph model (N=' + str(N) + ',m=' + str(m) + ')', fontsize=16)
    # Right Logaritmic degree distribution
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
    plt.subplot(132)
    a,b,anch = ploting_logscale(degree_sequence)
    plt.bar(a, b, width=anch, align='edge', edgecolor="black", color='b', lw=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Degree k')
    plt.title('Experimental Degree distribution (log-scale)')
    plt.show()


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
    # Power Law Distribution
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

    # Plot degree distribution Linear Scale
    di = list(G.degree)
    degree = [dix[1] for dix in di]
    plt.subplot(121)
    extra = plot_degree_distr(degree, 'Experimental Degree distribution')
    if pois == 0:
        distri = 'Poisson'
    else:
        distri = 'Power-Law'
    plt.suptitle('Configuration Model (CM) (N=' + str(N) + ',Distribution=' + distri + ')', fontsize=16)
    # Right Logaritmic degree distribution
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
    plt.subplot(122)
    a,b,anch = ploting_logscale(degree_sequence)
    plt.bar(a, b, width=anch, align='edge', edgecolor="black", color='b', lw=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Degree k')
    plt.title('Experimental Degree distribution (log-scale)')
    plt.show()
