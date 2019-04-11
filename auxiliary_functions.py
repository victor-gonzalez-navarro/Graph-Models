import numpy as np
import math
import matplotlib.pyplot as plt
import collections


# ------------------------------------------------------------------------------------------------------------ FUNCTIONS
def plot_degree_distr(d,title):
    degreeCount = collections.Counter(d)
    deg, cnt = zip(*degreeCount.items())
    cnt = cnt/np.sum(cnt)
    plt.bar(deg, cnt, width=0.60, color='b')
    plt.title(title)
    extra = max(cnt)+max(cnt)*0.1
    plt.ylim((0, extra))
    plt.xlabel('Degree k')
    plt.ylabel('Fraction of vertices Pk having degree k')
    return extra


def ploting_logscale(degree_sequence):
    degreeCount = collections.Counter(degree_sequence)
    x, y = zip(*degreeCount.items())

    nbins = 10
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
    x = tuple(lin)
    y = tuple(pdflog)
    x10 = [10 ** xv for xv in x]
    anch = ()
    for i in range(0, len(x) - 1):
        anch = anch + ((x10[i + 1] - x10[i]),)
    return x10[:-1],tuple(y),anch
