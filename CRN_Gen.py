from BFS import bfs
import pickle
import time
#from ReadMesh import readmesh
#from datfilepy import readdata
from Para_plot import*
def Generate():
    nCRN=2500
    t0=time.time()
    t1=time.clock()
    #readmesh()
    #readdata("CFDsol.dat")
    with open("mesh.pkl", 'rb') as f:
        graph = pickle.load(f)
    with open("zone.pkl", 'rb') as f:
        zone = pickle.load(f)
    with open("bc2f.pkl", 'rb') as f:
        bc = pickle.load(f)
    with open("pfaces.pkl", 'rb') as f:
        pfaces = pickle.load(f)
    print "Completed reading mesh"
    with open("data.pkl", 'rb') as f:
        graph_dat = pickle.load(f)
    with open("header.pkl", 'rb') as f:
        header = pickle.load(f)
    print "Completed reading data"

    red1=1
    red2=1
    criteria=['Static Temperature', 'Mole fraction of oh']
    tol = [0.01,0.01]
    increment=[0.05,0.05]
    print 'criteria:',criteria
    #while len(graph)>nCRN:
    for loop in range(300):
        length=len(graph)
        red1=red2
        print "Number of Clusters=",length
        if length<nCRN:
            break
        cluster=bfs(graph,graph_dat,zone,tol,header,criteria)

        del graph
        del graph_dat
        graph=cluster["graph"]
        graph_dat=cluster["info"]
        if loop>0:
            red2= length-len(graph)
            if red2==0:
                tol = [t + increment[tol.index(t)] for t in tol]
            elif red2<red1:
                tol = [t + increment[tol.index(t)] for t in tol]
                print "tol=",tol


    t01=time.time()-t0
    t11=time.clock() - t1
    print "time to cluster0=",t01
    print "time to cluster1=",t11
    with open('graph2plot.pkl','wb') as file:
        pickle.dump(graph,file,pickle.HIGHEST_PROTOCOL)
    with open('graphdata.pkl','wb') as file:
        pickle.dump(graph_dat,file,pickle.HIGHEST_PROTOCOL)
    #dict=graph2para(graph_dat, header)
    para3D(graph_dat, header)
    t02 = time.time() - t01
    t12 = time.clock() - t11
    print "time to para0=", t02
    print "time to para1=", t12
    #with open('cluster.pkl','wb') as file:
    #    pickle.dump(cluster,file,pickle.HIGHEST_PROTOCOL)
    return 0

if __name__=="__main__":
    Generate()