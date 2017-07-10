import matplotlib.pyplot as plt
import math


# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure(1)
# ax = Axes3D(fig)

def Fitness(quantity, dict):
    """

    :param quantity: 
    :param dict: 
    :return: 
    """

    mu = dict['mu']
    sigma = dict['sigma']
    comp=0
    for quant in range(len(quantity)):
        comp += abs(quantity[quant] - mu[quant]) / sigma[quant]
    if comp < 1*len(quantity):
        return True
    else:
        return False


def Enqueue(Q, elem):
    Q.append(elem)


def Dequeue(Q):
    elem = Q[0]
    Q.remove(Q[0])
    return elem


def bfs(graph, graph_dat, zone, tol, header, criteria):
    """
    Breadth first search of graph to cluster cells/reactors. 
    :param graph: graph of cells/reactors
    :param graph_dat: dictionary of data of cells/reactors
    :param zone: 
    :param tol: list of tolerances for splitting criteria
    :param header: list of header strings with index in graph = list.index+1
    :param criteria: List of splitting criteria
    :return: graph of clustered reactors
            cluster={'graph':{}, 'info':{}}
    """
    print header
    count = 1
    startpt = 1  # int(len(graph)/2)
    cluster = {}
    graph[startpt]['color'] = 'gray'
    graph[startpt]['distance'] = 0
    graph[startpt]['predecessor'] = 1
    Q = []
    avg=[]
    for cri in criteria:
        q_max = []
        q_min = []
        for z in graph_dat:
            if header.index(cri) + 1 in graph_dat[z]:
                quantities = graph_dat[z][header.index(cri) + 1]
                q_max.append(max(quantities))
                q_min.append(min(quantities))
        avg.append(tol[criteria.index(cri)] * (max(q_max) - min(q_min)))
    Enqueue(Q, startpt)
    # searching for cell zone in which cell with index 1 is present
    for i in zone:
        if zone[i][0] == 12:
            if startpt >= zone[i][1] and startpt <= zone[i][2]:
                break
    # initialising cluster data dictionary: same format as graph_dat
    zid = i
    cluster['info'] = {zid: {}}
    for j in graph_dat[i]:
        cluster['info'][zid][j] = [graph_dat[zid][j][0]]
    # print data_in
    mu=[]
    for cri in criteria:
        mu.append(cluster['info'][zid][header.index(cri) + 1][0])
    cluster["graph"] = {1: {"adjacent": [], 'color': 'white', 'distance': 'inf', \
                            'predecessor': '', 'faces': [], 'cells': [startpt], \
                            'sigma': avg, 'mu': mu, 'n': 1}}
    neighbours = {1: []}
    within = {}
    num_cell=0
    # BFS starts
    while Q != []:

        u = Dequeue(Q)

        for element in graph[u]['adjacent']:
            if element != 0:
                if graph[element]['color'] == 'white':
                    graph[element]['color'] = 'gray'
                    graph[element]['distance'] = graph[u]['distance'] + 1
                    graph[element]['predecessor'] = u
                    Enqueue(Q, element)
        # check for fitness
        # search for zone in which cell exists
        # print u
        for i in zone:
            if zone[i][0] == 12:
                if u >= zone[i][1] and u <= zone[i][2]:
                    break
        q1=[]
        for cri in criteria:
            spot = header.index(cri)
            q1.append(graph_dat[i][spot + 1][u - zone[i][1]])
        if (Fitness(q1, cluster['graph'][count]) and u != startpt):
            n = cluster['graph'][count]['n']
            mu = cluster['graph'][count]['mu']
            cluster['graph'][count]['n'] = n + 1
            cluster['graph'][count]['mu'] = [(q1[ind] + n * mu[ind]) / (n + 1) for ind in range(len(q1))]
            cluster['graph'][count]['sigma'] = avg
            #searchin for zone
            for i in zone:
                if zone[i][0] == 12:
                    if u >= zone[i][1] and u <= zone[i][2]:
                        break
            ##updating cluster quantities
            # Cell Volume
            vol = header.index('Cell Volume') + 1
            v1 = cluster['info'][zid][vol][count - 1]
            v2 = graph_dat[i][vol][u - zone[i][1]]
            v = v1 + v2
            cluster['info'][zid][vol][count - 1] = v
            # Mass
            density = header.index('Density') + 1
            d1 = cluster['info'][zid][density][count - 1]
            d2 = graph_dat[i][density][u - zone[i][1]]
            m1 = d1 * v1
            m2 = d2 * v2
            m = m1 + m2
            # Density
            index = header.index('Density') + 1
            cluster['info'][zid][index][count - 1] = m / (v)
            # Pressure
            index = header.index('Static Pressure') + 1
            p1 = cluster['info'][zid][index][count - 1]
            p2 = graph_dat[i][index][u - zone[i][1]]
            P = (m1 * p1 + m2 * p2) / m
            cluster['info'][zid][index][count - 1] = P
            # Enthalpy
            index = header.index('Specific Heat (Cp)') + 1
            cp1 = cluster['info'][zid][index][count - 1]
            cp2 = graph_dat[i][index][u - zone[i][1]]
            index = header.index('Static Temperature') + 1
            T1 = cluster['info'][zid][index][count - 1]
            T2 = graph_dat[i][index][u - zone[i][1]]
            H1 = m1 * cp1 * T1
            H2 = m2 * cp2 * T2
            H = H1 + H2
            index = header.index('Enthalpy') + 1
            cluster['info'][zid][index][count - 1] = H
            if H < 0:
                print H
                print count
                print H1
                print H2
            # Moles
            index = header.index('Static Temperature') + 1
            T1 = cluster['info'][zid][index][count - 1]
            T2 = graph_dat[i][index][u - zone[i][1]]
            R = 8314
            n1 = (p1 + 1e5) * v1 / (R * T1)
            n2 = (p2 + 1e5) * v2 / (R * T2)
            n = n1 + n2
            # cluster['info'][zid][index][count - 1] = n1 +n2
            # Cp
            index = header.index('Specific Heat (Cp)') + 1
            Cp = (H / (m * (P + 1e5))) * ((n * R) / v)
            cluster['info'][zid][index][count - 1] = Cp
            # Temperature
            index = header.index('Static Temperature') + 1
            T = H / (m * Cp)
            cluster['info'][zid][index][count - 1] = T
            # T1 = cluster['info'][zid][index][count - 1]
            # T2 = graph_dat[i][index][u - zone[i][1]]
            # T = (m1 * T1 + m2 * T2) / m
            # cluster['info'][zid][index][count - 1] = T
            # Species
            for h in header:
                if h.find('Mole Fraction') != -1:
                    index = header.index(h) + 1
                    X1 = cluster['info'][zid][index][count - 1]
                    X2 = graph_dat[i][index][u - zone[i][1]]
                    cluster['info'][zid][index][count - 1] = (n1 * X1 + n2 * X2) / n
            # X-Coordinate
            index = header.index('X-Coordinate') + 1
            x1 = cluster['info'][zid][index][count - 1]
            x2 = graph_dat[i][index][u - zone[i][1]]
            x = (m1 * x1 + m2 * x2) / m
            cluster['info'][zid][index][count - 1] = x
            # Y-Coordinate
            index = header.index('Y-Coordinate') + 1
            y1 = cluster['info'][zid][index][count - 1]
            y2 = graph_dat[i][index][u - zone[i][1]]
            y = (m1 * y1 + m2 * y2) / m
            cluster['info'][zid][index][count - 1] = y
            # Z-Coordinate
            index = header.index('Z-Coordinate') + 1
            z1 = cluster['info'][zid][index][count - 1]
            z2 = graph_dat[i][index][u - zone[i][1]]
            z = (m1 * z1 + m2 * z2) / m
            cluster['info'][zid][index][count - 1] = z
            ##Adding faces: set() allows only one occurence of common faces
            s1 = set(cluster['graph'][count]['faces'])
            s2 = set(graph[u]['faces'])
            # union of sets
            cluster['graph'][count]['faces'] = list(s1 | s2)
            # Adding cells
            c1 = cluster['graph'][count]['cells']
            c2 = graph[u]['cells']
            cluster['graph'][count]['cells'] = c1 + c2
        elif u != startpt:
            # creating new cluster
            count += 1
            # initialising cluster data dictionary: same format as graph_dat
            zid = i
            for j in graph_dat[i]:
                try:
                    cluster['info'][zid][j].append(graph_dat[zid][j][u - zone[zid][1]])
                except KeyError:
                    cluster['info'][zid] = {j: [graph_dat[zid][j][u - zone[zid][1]]]}
            mu=[]
            for cri in criteria:
                data_list = cluster['info'][zid][header.index(cri) + 1]
                mu.append(data_list[len(data_list) - 1])
            cluster['graph'][count] = {"adjacent": [], 'color': 'white', \
                                       'distance': 'inf', 'predecessor': '', 'faces': graph[u]['faces'],
                                       'cells': graph[u]['cells'], \
                                       'sigma': avg, 'mu': mu, 'n': 1}

        try:
            neighbours[count] = list(set(graph[u]['adjacent'] + neighbours[count]))
        except KeyError:
            neighbours[count] = graph[u]['adjacent']
        if u in neighbours[count]:
            neighbours[count].remove(u)
        if 0 in neighbours[count]:
            neighbours[count].remove(0)

        within[u] = count
        # check for adjacent cells in neighbouring clusters
        """for g in neighbours.keys():
            if u in neighbours[g]:
                cluster["graph"][g]["adjacent"] = list(set(cluster["graph"][g]["adjacent"] + [count]))
                cluster["graph"][count]["adjacent"] = list(set(cluster["graph"][count]["adjacent"] + [g]))"""
        graph[u]['color'] = 'black'
        #num_cell+=1
        #print "num_cell=", num_cell
    for g in neighbours:
        for c in neighbours[g]:
            if g!=within[c]:
                cluster["graph"][g]["adjacent"].append(within[c])

    return cluster
