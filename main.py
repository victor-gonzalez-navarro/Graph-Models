import networkx as nx
import matplotlib.pyplot as plt
import collections
import numpy as np
import pandas as pd
import math


def ploting_linearscale(deg, pdf, fig):
    pdfnew = ();
    degnew = ();
    counter = 0
    for i in range(max(deg), min(deg) - 1, -1):
        if i == deg[counter]:
            pdfnew = pdfnew + (pdf[counter],)
            counter = counter + 1
        else:
            pdfnew = pdfnew + (0,)
        degnew = degnew + (i,)
    ccdfnew = tuple([np.sum(pdfnew[0:i + 1]) for i in range(len(degnew))])
    ploting2_linearscale(deg, pdf, 1, 'Degree PDF (Linear-scale)', 'Fraction of vertices Pk having degree k', 1, degnew)
    ploting2_linearscale(degnew, ccdfnew, 2, 'Degree CCDF (Linear-scale)',
                         'Fraction of vertices Pk having degree k or greater', 1, degnew)


def ploting2_linearscale(x, y, numsubplot, title, ylab, anch, xnew):
    ax = fig.add_subplot(1, 2, numsubplot)
    ax.bar(x, y, width=anch, align='edge', color='b', edgecolor="black", lw=2)
    # ax.yaxis.grid()
    plt.title(title)
    # plt.xlim([3,100])
    plt.ylabel(ylab)
    plt.xticks(xnew)
    plt.xlabel('Degree k')
    # plt.grid()


def ploting_logscale(x, y, fig):
    nbins = 10;
    pdflog = np.zeros(nbins)
    lin = np.linspace(math.log10(min(x)), math.log10(max(x) + 1), nbins + 1)
    logx = [math.log10(vall) for vall in x]
    num = 0
    for vall in logx:
        interval = 0
        stop = False
        for vall2 in lin:
            if stop == False:
                if vall - vall2 < 0:
                    stop = True
                interval = interval + 1
            if stop == False and vall2 == lin[-1]:
                interval = nbins
        pdflog[interval - 2] = pdflog[interval - 2] + y[num]
        num = num + 1
    pdflog = pdflog / (np.sum(y))
    x = tuple(lin);
    y = tuple(pdflog)
    x10 = [10 ** xv for xv in x]
    anch = ()
    for i in range(0, len(x) - 1):
        anch = anch + ((x10[i + 1] - x10[i]),)
    ccdf2 = tuple([np.sum(y[i:len(y)]) for i in range(len(y))])
    ploting2_logscale(x10[:-1], tuple(y), 1, 'Degree PDF (Log-scale)', 'Fraction of vertices Pk having degree k', anch)
    ploting2_logscale(x10[:-1], tuple(ccdf2), 2, 'Degree CCDF (Log-scale)',
                      'Fraction of vertices Pk having degree k or greater', anch)


def ploting2_logscale(x, y, numsubplot, title, ylab, anch):
    ax = fig.add_subplot(1, 2, numsubplot)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.bar(x, y, width=anch, align='edge', color='b', edgecolor="black", lw=2)
    plt.title(title)
    plt.ylabel(ylab)
    plt.xlabel('Degree k')
    plt.grid()


# ------------------------------------------------------------------------------------------------------------ MAIN CODE
# graph_names =  ['real/airports_UW']
graph_names = ['real/airports_UW']
#graph_names = ['model/ER1000k8']
# graph_names = ['toy/20x2+5x2','toy/circle9','toy/graph3+1+3','toy/graph3+2+3','toy/grid-p-6x6','toy/rb25','toy/star','toy/wheel']
# graph_names = graph_names + ['model/256_4_4_2_15_18_p','model/256_4_4_4_13_18_p','model/BA1000','model/ER1000k8','model/ER5000-kmed8','model/homorand_N1000_K4_0','model/homorand_N1000_K6_0','model/rb125','model/SF_1000_g2.5','model/SF_1000_g2.7','model/SF_1000_g3.0','model/SF_500_g2.7','model/ws1000','model/ws2000']
# graph_names = graph_names + ['real/airports_UW','real/dolphins','real/PGP','real/zachary_unwh']
information = np.zeros((len(graph_names), 9))
dict_ploting = dict()

for i in range(len(graph_names)):
    # for i in range(1):
    G = nx.read_pajek('A1-networks/' + graph_names[i] + '.net')
    # G = nx.random_powerlaw_tree(2000, gamma=3, seed=None, tries=10000)

    # Number of nodes
    Nnodes = len(G.node)

    # Plotting
    if Nnodes < 0:
        nx.draw(G, with_labels=True, font_weight='bold')
        plt.show()
    information[i, 0] = Nnodes

    # Number of edges
    Nedges = len(G.edges)
    information[i, 1] = Nedges

    # Minimum, maximum and average degree
    avg_degree = 0
    nodesdegree = [n[1] for n in G.degree]
    min_degree = min(nodesdegree)
    max_degree = max(nodesdegree)
    avg_degree = np.mean(nodesdegree)
    information[i,2] = min_degree; information[i,3] = max_degree; information[i,4] = avg_degree

    # Clustering Coefficient
    # dict_information = dict()
    # for key1,value1 in G._adj.items():
    #    connections = ()
    #    for key2,value2 in value1.items():
    #        connections = connections + (key2,)
    #    dict_information[key1] = connections

    # CC = 0
    # for vertex,conected_vertices in dict_information.items():
    #    #print('**'+vertex+'**')
    #    dv = len(conected_vertices)
    #    if dv != 1:
    #        nv = 0
    #        for connected_vertex in conected_vertices:
    #            for connected_vertex2 in conected_vertices:
    #                if connected_vertex2 in dict_information[connected_vertex]:
    #                    nv = nv + 1
    #        nv = nv/2
    #        CCv = (2*nv)/(dv*(dv-1))
    #        CC = CC + CCv
    #    else:
    #        CCv = 0
    # CC = CC / Nnodes

    Gundirected = nx.Graph(G)
    CC = nx.average_clustering(Gundirected)
    information[i, 5] = CC

    # Assortativity
    if min_degree == max_degree:
        # print('The Assortativity is: 1')
        assorta = 1
    else:
        assorta = nx.degree_assortativity_coefficient(Gundirected)
    information[i, 6] = assorta

    # Average Shortest Path
    avg_shortestpath = nx.average_shortest_path_length(Gundirected)
    information[i, 7] = avg_shortestpath

    # Diameter
    diameter = nx.diameter(Gundirected)
    information[i, 8] = diameter

    # Degree distribution
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())
    pdf = tuple([val / np.sum(cnt) for val in cnt])

    # Right Linear degree distribution and Complementary Cumulative Distribution Function
    fig = plt.figure(figsize=(20, 4))
    ploting_linearscale(deg, pdf, fig)

    # Right Logaritmic degree distribution and Complementary Cumulative Distribution Function
    fig = plt.figure(figsize=(20, 4))
    ploting_logscale(deg, cnt, fig)
    plt.show()

    dict_ploting[i] = graph_names[i] + '.net'

# Print table
df = pd.DataFrame(information)
df.columns = ['Num. nodes', 'Num. edges', 'Min. degree', 'Max. degree', 'Avg. degree', 'Avg. clust. coef.',
              'Assortativity', 'Avg. path length', 'Diameter']
cols = ['Num. nodes', 'Num. edges', 'Min. degree', 'Max. degree', 'Diameter']
df[cols] = df[cols].applymap(np.int64)
pd.options.display.float_format = '{:,.4f}'.format
df.rename(index=dict_ploting, inplace=True)