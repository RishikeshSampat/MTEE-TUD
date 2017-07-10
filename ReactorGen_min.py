import math
import pickle
import time
import cantera as ct
from numpy import dot, array, add, matrix, linalg, subtract, ndarray, matmul, float64
import SaveReactors as sr

class Transport(object):
    def __init__(self,Reactors):
        self.Reactors=Reactors

def advance_solve(X, PSR, dt, gas_obj):
    num_spec=len(X)
    error=0
    gas_obj.TPX=PSR.thermo.T, PSR.thermo.P, X
    net = ct.ReactorNet([PSR])
    xin = PSR.thermo.X
    net.advance(dt)
    xfin = PSR.thermo.X
    func = [abs(xfin[n] - xin[n]) for n in range(num_spec)]
    for num in func:
        error+=num
    #print "func val=",error
    #returns summation of errors of all species over reactor
    #print "Working"
    return error
def constraint_func(X, PSR, dt, gas_obj):
    return (sum(X)-1)
def Gen(energy):
    print "start"
    with open("graph2plot.pkl", 'rb') as f:
        graph = pickle.load(f)
    with open("header.pkl", 'rb') as f:
        header = pickle.load(f)
    with open("zone.pkl", 'rb') as f:
        zone = pickle.load(f)

    print "start program"
    PSR=[]
    MFC=[]
    Res=[]
    key=[]

    #species taken from gri-3.0
    species=['H2','H','O','O2','OH','H2O','HO2','H2O2','C','CH',\
             'CH2','CH2(S)','CH3','CH4','CO','CO2','HCO','CH2O',\
             'CH2OH','CH3O','CH3OH','C2H','C2H2','C2H3','C2H4','C2H5',\
             'C2H6','HCCO','CH2CO','HCCOH','N','NH','NH2','NH3','NNH',\
             'NO','NO2','N2O','HNO','CN','HCN','H2CN','HCNN','HCNO','HOCN',\
             'HNCO','NCO','N2','AR','C3H7','C3H8','CH2CHO','CH3CHO']
    with open("facearea.pkl", 'rb') as f:
        facearea = pickle.load(f)
    with open("graphdata.pkl", 'rb') as f:
        data = pickle.load(f)
    with open("data_std.pkl", 'rb') as f:
        data_std = pickle.load(f)
    with open("pfaces.pkl", 'rb') as f:
        periodic = pickle.load(f)
    print header
    for i in data[19].keys():
        s1 = header[i-1].find('Mole fraction of')
        s2 = header[i - 1].find('<')
        s3 = header[i - 1].find('>')
        name=header[i-1][len('Mole fraction of '):]
        if s1!=-1 and s2==-1 and s3==-1 and (name.upper() in species):
            key.append(name)
    #creating ideal reactors
    print "Creating Reactors"
    g = ct.Solution('gri30.cti')
    volume_total = 0
    for cell in graph.keys():
        T=data[19][header.index('Static Temperature')+1][cell-1]#data_std[19][3][cell-1]#3 is xfile.h index for temperature
        P=data[19][header.index('Static Pressure')+1][cell-1]+ct.one_atm#data_std[19][1][cell-1]+ct.one_atm#1 is xfile.h index for pressure
        X={j.upper():data[19][header.index('Mole fraction of '+j)+1][cell-1] for j in key}
        if T<0:
            print "Negative temp=",T
            print "In cell=",cell
            print "Out of=",len(graph.keys())
            print "Adjacent=",graph[cell]['adjacent']
        g.TPX=T,P,X
        psr=ct.IdealGasReactor(g, energy=energy)
        vol=data[19][header.index('Cell Volume')+1][cell-1]
        #print vol
        psr.volume=vol
        PSR.append(psr)
        volume_total += vol
        del psr
    #generating mass flow controllers
    #also generating walls for heat transfer
    print "Generating MFCs"
    Reserv={}
    res_num=0
    mfc_flow=0
    interior_faces=0
    cells_react=[]
    valves=[]
    MassImb = [0]*len(PSR)
    for cell1 in graph:
        cells_react.append(len(graph[cell1]['cells']))
        graph[cell1]['mf']=[]
        #neighbouring clusters
        for link in graph[cell1]['adjacent']:
            if link==0:
                print graph[cell1]['adjacent']
            #list of  common faces between reactors
            s1 = set(graph[cell1]['faces'])
            s2 = set(graph[link]['faces'])
            face=list(s1 & s2)
            interior_faces+=len(face)

            area=0
            mflow=0
            v=array([0,0,0])
            for f in face:
                #c0 is LHS cell of face
                c0 = facearea[f][2][0]
                mflux = data_std[20][18][f - 1]
                #if cell on LHS of face in reactor then proceed
                if c0 in graph[cell1]['cells']:
                    #print mflux
                    #print facearea[f]
                    mflow += mflux#*facearea[f][0]
                    #area+=facearea[f][0]
                #else:
                 #   mflow += -mflux
            #print mflow
            graph[cell1]['mf'].append(mflow)
            #creating reservoir with adjacent reactor conditions
            g.TPX = PSR[link - 1].thermo.TPX
            Res_temp = ct.Reservoir(g)
            res_num += 1
            Reserv[res_num]=[Res_temp,PSR[link - 1]]

            MassImb[cell1-1]=MassImb[cell1-1]+mflow

            if mflow>0:
                #condition when mass flow going out of current reactor
                power = math.log10(PSR[cell1 - 1].volume)
                #valve created at exit of current reactor
                val = ct.Valve(PSR[cell1 - 1], Res_temp, K=10 ** power)
                valves.append(val)
                if energy=='on':
                    wall=ct.Wall(PSR[cell1 - 1], Res_temp,U=0, Q=0)
                #creating reservoir with current reactor conditions
                g.TPX = PSR[cell1 - 1].thermo.TPX
                Res_temp = ct.Reservoir(g)
                res_num += 1
                Reserv[res_num] = [Res_temp, PSR[cell1 - 1]]
                #MFC created at inlet of adjacent reactor
                mfc = ct.MassFlowController(Res_temp, PSR[link - 1], mdot = mflow)
                MFC.append(mfc)
                mfc_flow += mflow
                if energy=='on':
                    wall=ct.Wall(Res_temp, PSR[link - 1],U=0, Q=0)


            else:
                #condition when mass flow enterin reactor
                #MFC created at inlet of current reactor
                mfc = ct.MassFlowController(Res_temp, PSR[cell1 - 1], mdot = -mflow)
                MFC.append(mfc)
                mfc_flow += mflow
                if energy=='on':
                    wall=ct.Wall(Res_temp, PSR[cell1 - 1],U=0, Q=0)
                #creating reservoir with current reactor conditions
                g.TPX = PSR[cell1 - 1].thermo.TPX
                Res_temp = ct.Reservoir(g)
                res_num += 1
                Reserv[res_num] = [Res_temp, PSR[cell1 - 1]]
                power = math.log10(PSR[link - 1].volume)
                #valve created at exit of adjacent reactor
                val = ct.Valve(PSR[link - 1], Res_temp, K=10 ** power)
                # MFC.append(mfc)
                # mfc_flow+=mflow
                valves.append(val)
                if energy == 'on':
                    wall = ct.Wall(PSR[link - 1], Res_temp, U=0, Q=0)



    #Boundary conditions
    MFC_inlet=[]
    MFC_periodic=[]
    Exhaust=[]
    PC=[]
    mflow_tot=0
    inlet_faces=0
    f_area_tot=0
    OutletReact=[]
    InletReact=[]
    print "Solving BCs"
    for cell in graph:
        for face in graph[cell]['faces']:
            for z in zone:
                if zone[z][0]==13 and face>=zone[z][1] and face<=zone[z][2]:
                    #mass flow inlet
                    # requires mass flow controller and reservoir
                    #print zone[z]
                    if zone[z][3]==20:
                        #print data[z].keys()
                        inlet_faces+=1
                        f_area=facearea[face][0]
                        #f_area=data[z][header.index('Face Area Magnitude')+1][face-zone[z][1]]
                        f_area_tot+=abs(f_area)
                        mflux = data_std[z][18][face - zone[z][1]]
                        mflow = abs(mflux)
                        if mflow==0:
                            print face," ",z
                        mflow_tot+=mflux
                        T = data_std[z][3][face - zone[z][1]]
                        P = data_std[z][1][face - zone[z][1]]+ct.one_atm
                        X = {j.upper(): data[19][header.index('Mole fraction of ' + j) + 1]\
                            [cell - zone[19][1]] for j in key}
                        g.TPX = T, P, X
                        res=ct.Reservoir(g)
                        Res.append(res)
                        mfc=ct.MassFlowController(res,PSR[cell - 1],mdot = mflow)
                        MFC_inlet.append(mfc)
                        MassImb[cell - 1] = MassImb[cell - 1] + mflow
                        InletReact.append(cell-1)
                    #periodic boundary
                    if zone[z][3]==12:
                        face_area=facearea[face][0]
                        mflux=data_std[z][18][face-zone[z][1]]#18:mflux id
                        #mflow=mflux*face_area
                        mflow=mflux
                        if face in periodic:
                            pair=periodic[face]#second face in pair, i.e flow towards it
                            for cell_check in graph:
                                if pair in graph[cell_check][face]:
                                    g.TPX = PSR[cell_check - 1].thermo.TPX
                                    Res_temp = ct.Reservoir(g)
                                    res_num += 1
                                    Reserv[res_num] = [Res_temp, PSR[cell_check - 1]]
                                    MassImb[cell - 1] = MassImb[cell - 1] + mflow
                                    if mflow > 0:
                                        # mfc=ct.MassFlowController(PSR[cell1-1][0],Res_temp,mdot=mflow)
                                        power = math.log10(PSR[cell - 1].volume)
                                        val = ct.Valve(PSR[cell - 1], Res_temp, K=10 ** power)
                                        # MFC.append(mfc)
                                        # mfc_flow+=mflow
                                        valves.append(val)
                                        g.TPX = PSR[cell - 1].thermo.TPX
                                        Res_temp = ct.Reservoir(g)
                                        res_num += 1
                                        Reserv[res_num] = [Res_temp, PSR[cell - 1]]
                                        mfc = ct.MassFlowController(Res_temp, PSR[cell_check - 1], mdot= mflow)
                                        MFC.append(mfc)
                                        mfc_flow += mflow

                                    else:
                                        mfc = ct.MassFlowController(Res_temp, PSR[cell - 1], mdot=-mflow)
                                        MFC.append(mfc)
                                        mfc_flow += mflow
                                        g.TPX = PSR[cell - 1].thermo.TPX
                                        Res_temp = ct.Reservoir(g)
                                        res_num += 1
                                        Reserv[res_num] = [Res_temp, PSR[cell - 1]]
                                        power = math.log10(PSR[cell_check - 1].volume)
                                        val = ct.Valve(PSR[cell_check - 1], Res_temp, K=10 ** power)
                                        # MFC.append(mfc)
                                        # mfc_flow+=mflow
                                        valves.append(val)
                                    break

                    #wall
                    if zone[z][3]==3:
                        #Res_env
                        T = 300
                        P = ct.one_atm
                        g.TP = T, P
                        Res_env = ct.Reservoir(g)
                        w=ct.Wall(PSR[cell - 1], Res_env, U=0, Q=0)

                    #outlet
                    if zone[z][3]==36:
                        mflux = data_std[z][18][face - zone[z][1]]  # 18:mflux id
                        mflow = mflux
                        MassImb[cell - 1] = MassImb[cell - 1] + mflow
                        T = data[19][header.index('Static Temperature') + 1][cell - zone[19][1]]
                        P = data[19][header.index('Static Pressure') + 1][cell - zone[19][1]]+ct.one_atm
                        X = {j.upper(): data[19][header.index('Mole fraction of ' + j) + 1] \
                            [cell - zone[19][1]] for j in key}
                        #T=300
                        #P=ct.one_atm
                        #X={'O2':1,'N2':3.76}
                        g.TPX = T, P, X
                        exhaust = ct.Reservoir(g)
                        power = math.log10(PSR[cell-1].volume)
                        pc=ct.Valve(PSR[cell-1],exhaust,K=10**power)
                        Exhaust.append(exhaust)
                        PC.append(pc)
                        OutletReact.append(cell-1)


    print "Saving initial state"

    sr.SaveData(cells_react,'ReactorCells')
    print "Total inlet flow=",mflow_tot
    print "Inlet faces=",inlet_faces
    print "Inlet area total=",f_area_tot
    print "Total Interior faces=", interior_faces
    print "MFC total flow=",mfc_flow
    print "Total Volume=", volume_total
    print "Solving network"
    #net=ct.ReactorNet(PSR)

    #print "Abs Tol",net.atol
    #net.verbose=True
    Reserv[1][0].insert(Reserv[1][1].thermo)
    print "Advancing"
    dt=0.001
    t0 = time.time()
    t1 = time.clock()
    error_log=[]
    error_prev=1
    bound = []

    for s in PSR[0].thermo.X:

        bound.append((0, 1))

    for iter in range(10000):
        print "Iteration ", iter

        # Solving individual reactors
        for iter in range(1000):
            error = 0
            denom = 0
            X_initial = []
            func_initial = []
            for r in PSR:
                net = ct.ReactorNet([r])


                tolx=0.00001
                try:
                    e_internal = 0
                    xin = r.thermo.X
                    net.advance(dt)
                    xfin = r.thermo.X
                    func = [(xfin[n] - xin[n]) for n in range(len(xin))]
                    e = [abs(i) for i in func]
                    for num in range(len(e)):
                        e_internal += e[num]
                        denom += xfin[num]


                    """if (e_internal)<tolx:
                        #print "Internal iter=",iter
                        break"""
                except(RuntimeError):
                    print "Reactor num=", PSR.index(r)
                    print r.thermo.X
                    print r.thermo.P
                    print r.thermo.T
                X_initial.append(xin)
                func_initial.append(func)
                e = [abs(i) for i in func]
                for num in range(len(e)):
                    error += e[num]
                    denom += xfin[num]
            # updating corresponding reservoirs
            for num in Reserv:
                Reserv[num][0].insert(Reserv[num][1].thermo)
            # Residual analysis
            error = error / denom
            error_reduc = error / error_prev
            error_prev = error
            #print error
            #print "int_reduc=", error_reduc
            if error_reduc > 0.95:
                #print "error=", error
                break
            error_log.append(error)

        # Residual analysis
        print error
        print "reduc=", error_reduc
        if error < 1:
            print "error=", error

            break
        error_log.append(error)

    print "Solution completed"
    t01 = time.time() - t0
    t11 = time.clock() - t1
    print "time to cluster0=", t01
    print "time to cluster1=", t11
    sr.Save(PSR, header, data, MassImb, OutletReact, InletReact)
    import matplotlib.pyplot as plt
    plt.figure(2)
    plt.plot(error_log)
    plt.title('Convergence error d2f')
    plt.xlabel('Iteration')
    plt.ylabel('Error')
    plt.show()

    with open('errorlog.pkl','wb') as file:
        pickle.dump(error_log,file,pickle.HIGHEST_PROTOCOL)


if __name__=="__main__":


    Gen(energy='off')