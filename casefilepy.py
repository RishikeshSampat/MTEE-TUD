#!/usr/bin python

import sys
import re
from struct import*
import binascii
from copy import copy
from numpy import zeros, array, sign, cross, dot, ones, arctan, sin, cos, pi, mod, sqrt,linalg,subtract
import pickle
import math

print "Converting from ANSYS Fluent format graph"

# Use regular expressions to identify sections and tokens found in a fluent file
re_dimline  = re.compile(r"\(2\s(\d)\)")
re_comment  = re.compile(r"\(0\s.*")
re_zone0    = re.compile(r"\(10\s\(0\s(\w+)\s(\w+)\s(\d+)\s(\d+)\)\)")
re_zone     = re.compile(r"\(10\s\((\w+)\s(\w+)\s(\w+)\s(\d+)\s(\d)\)(\(|)")
re_face0    = re.compile(r"\(13(\s*)\(0\s+(\w+)\s+(\w+)\s+(0|0 0)\)\)")
re_face     = re.compile(r"\(13(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\)(\s*)(\(|)")
re_periodic = re.compile(r"\(18.*\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\).*\(")
re_pfaces   = re.compile(r"((^\s)|)(\w+)(\s*)(\w+)")
re_cells0   = re.compile(r"\(12(\s*)\(0(\s+)(\w+)(\s+)(\w+)(\s+)(0|0 0)\)\)")
re_cells    = re.compile(r"\(12.*\((\w+)\s+(\w+)\s+(\w+)\s+(\d+)\s+(\d+)\)\)")
re_cells2   = re.compile(r"\(12(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\)(\s*)(\(|)")
re_zones    = re.compile(r"\((45|39)\s+\((\d+)\s+(\S+)\s+(\S+).*\)\((.*|[0-9]+[\.]*[0-9]*)\)\)")
re_parant   = re.compile(r"(\s*\)(\s*)|\s*\)\)(\s*)|\s*\(\s*)")
re_cell_part_id=re.compile(r"\(40(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s*\)")
# The fluent mesh (the .msh file) is basically stored as a list of nodes, and then a
# list of faces for each zone of the mesh, the interior and the boundaries.

# Declare som maps that will be built when reading in the lists of nodes and faces:
cell_map = {}               # Maps cell id with nodes
cell_face_map = {}          # Maps cell id with faces
face_cell_map = {}          # Maps face id with two cells
face_list = []              # List of faces [[id, 2-4 nodes, 2 connecting cells and type]]
face_map = {}               # For each cell a dictionary with key=face and val=local face number
neighbour={}                # neighbouring cells
nodes = {}                  # Maps node number to coordinates
pfaces={}
face_area_list={}
face_node={}                #maps face id to nodes

# Information about connectivity and boundaries
boundary_cells2face = {}         # {cellnum:{zone_id:[face1,face2..],..},..}

# Some global values
num_cells = {}              # Total number of cells in different zones
zone={}                     # {id:[cell(12),first index, last index], id:[face(13),fid,lid,type]..}
zones = {}                  # zone information
zone_number_of_faces = {}   # number of faces for each zone

def read_periodic(ifile):
    """Scan past periodic section"""
    while 1:
        line = ifile.readline()
        a = re.search(re_pfaces, line)
        if a:
            pfaces[a.group(1)]=a.group(2)
            continue
        break

def read_zone_nodes(dim, Nmin, Nmax, ifile):
    """Scan lines for nodes and return in an array."""
    line = ifile.readline()
    readline = False

    if re.search(re_parant, line): # check for initial paranthesis
        readline = True
        #dummy = lines.pop(0)
    for i in range(Nmin, Nmax + 1):
        if readline:
            line = ifile.readline()
        string=''
        if re.search(re_parant,line):
            for x in line:
                if x!=')':
                    string=string+x
            line=string
        readline = True
        nodes[i] = [float(x) for x in line.split()]

def face_area():
    """
    calculate face area by splitting polygon into triangles. Assumption that
    nodes in cyclic order in list. Starting node kept as fixed vertex. Other two vertices 
    selected by traversing thorugh list.
    :return: 
    """
    xcoord = []
    ycoord = []
    zcoord = []
    for i in face_node:
        nd=face_node[i]
        if nd==[]:
            continue
        pt1=array(nodes[nd[0]])
        area=0
        for j in range(1,len(nd)-3+2):
            pt2=array(nodes[nd[j]])
            pt3=array(nodes[nd[j+1]])
            a=linalg.norm(pt1-pt2)
            b=linalg.norm(pt2-pt3)
            c=linalg.norm(pt1-pt3)
            s=(a+b+c)/2
            area+=math.sqrt(s*(s-a)*(s-b)*(s-c))
        v1 = subtract(pt1, pt3)
        v2 = subtract(pt3, pt2)
        v = cross(v1, v2)
        face_area_list[i] = [area, v, face_cell_map[i]]

def read_faces(zone_id, Nmin, Nmax, bc_type, face, ifile):
    """Read all faces and create cell_face_map + some boundary maps."""
    line = ifile.readline()
    readline = False
    if re.search(re_parant, line): # check for initial paranthesis
        readline = True

    ls = []
    face_index=Nmin
    for i in range(Nmin, Nmax + 1):
        if readline:
            line = ifile.readline()
        readline = True
        string=''
        lineorig=line
        if re.search(re_parant,line):
            for x in line:
                if x!=')':
                    string=string+x
                #if x==')':
                   # print "last line",i
            line=string
        ln = line.split()
        if face == 0:
            nd = int(ln[0]) # Number of nodes
            nds = [int(x, 16) for x in ln[1:(nd + 1)]]
            cells = [int(x, 16) for x in ln[(nd + 1):]]
            face_index=i
        else:
            nd = face
            nds = [int(x, 16) for x in ln[:nd]]
            cells = [int(x, 16) for x in ln[nd:]]

        face_node[i]=nds
        face_cell_map[i]=cells

        if min(cells) == 0: # A boundary zone
            num=max(cells)
            if num in boundary_cells2face:
                if zone_id in boundary_cells2face[num]:
                    boundary_cells2face[num][zone_id].append(i)
                else:
                    boundary_cells2face[num][zone_id]=[i]
            else:
                boundary_cells2face[num]={zone_id:[i]}

        for c in cells:
            if c > 0 and nds!=[]:
                if not c in cell_map:
                    cell_map[c] = copy(nds)
                    neighbour[c] = {}
                    neighbour[c]['adjacent']= copy(cells)
                    neighbour[c]['faces'] = [i]
                else:
                    cell_map[c] = list(set(cell_map[c] + nds))
                    neighbour[c]['adjacent']= list(set(neighbour[c]['adjacent'] + cells))
                    neighbour[c]['faces'].append(i)
                neighbour[c]['adjacent'].remove(c)
                neighbour[c]['color']= 'white'
                neighbour[c]['distance'] = 'inf'
                neighbour[c]['predecessor'] = ''
                neighbour[c]['cells'] = [c]



def scan_fluent_mesh(ifile):
    """Scan fluent mesh and generate numerous maps."""
    # Warning! Not yet tested for multiple interior zones
    dim = 0
    one = 0
    num_faces = 0
    while 1:
        line = ifile.readline()
        if len(line) == 0:
            print 'Finished reading file\n'
            break

        if dim == 0: # Dimension usually comes first
            a = re.search(re_dimline, line)
            if a:
                print 'Reading dimensions\n'
                dim = int(a.group(1))
                print 'Mesh is ' + str(dim) + 'D\n'
                continue

        if one == 0: # The total number of nodes
            a = re.search(re_zone0, line)
            if a:
                print 'Reading zone info\n'
                one, num_vertices, dummy1, dummy2 = int(a.group(1)), \
                     int(a.group(2), 16), int(a.group(3), 16), int(a.group(4))
                continue

        a = re.search(re_zone, line) # Nodes
        if a:
            zone_id, first_id, last_id, dummy1, dummy2 = int(a.group(1), 16), \
                int(a.group(2), 16), int(a.group(3), 16), int(a.group(4)), \
                int(a.group(5))
            print 'Reading ', last_id - first_id + 1,' nodes in zone ', zone_id + 1, '\n'
            read_zone_nodes(dim, first_id, last_id, ifile)
            continue

        a = re.search(re_zones,line) # Zone info
        if a:
            print 'Reading zone ', line
            dummy, zone_id, zone_type, zone_name, radius =  \
                       int(a.group(1)), int(a.group(2)),  a.group(3), \
                       a.group(4), a.group(5)
            zones[zone_id] = [zone_type, zone_name, radius]
            continue

        a = re.search(re_cells0, line) # Get total number of cells/elements
        if a:
            print 'Reading cell info ', line
            first_id, tot_num_cells = int(a.group(3),16), int(a.group(5), 16)
            continue

        a = re.search(re_cells,line) # Get the cell info.
        if a:
            zone_id, first_id, last_id, bc_type, element_type = \
                int(a.group(1), 16), int(a.group(2), 16), int(a.group(3), 16), \
                int(a.group(4), 16), int(a.group(5), 16)
            print 'Reading ', last_id - first_id + 1,' cells in zone ', zone_id, '\n'
            if last_id == 0:
                raise TypeError("Zero elements!")
            num_cells[zone_id] = [first_id, last_id, bc_type, element_type]
            zone[zone_id]=[12,first_id,last_id]
            continue

        a = re.search(re_cells2,line) # Get the cell info.
        if a:
            print line
            raise TypeError("Wrong cell type. Can only handle one single cell type")

        a = re.search(re_face0, line)
        if a:
            print 'Reading total number of faces\n', line
            num_faces = int(a.group(3),16)
            continue

        a = re.search(re_face, line)
        if a:
            print 'Reading faces ', line
            print line
            zone_id, first_id, last_id, bc_type, face_type = \
                 int(a.group(2), 16), int(a.group(3), 16), int(a.group(4),16), \
                 int(a.group(5), 16), int(a.group(6), 16)
            #print zone_id
            #print first_id
            #print last_id
            read_faces(zone_id, first_id, last_id, bc_type, face_type, ifile)
            zone_number_of_faces[zone_id] = last_id - first_id + 1
            zone[zone_id]=[13,first_id,last_id,bc_type]
            continue

        a = re.search(re_cell_part_id, line)
        if a:
            print 'Reading cell partition id ', line
            line = ifile.readline()
            if re.search(re_parant, line): # check for initial paranthesis
                while(1):
                    #read through lines
                  line = ifile.readline()
                  if re.search(re_parant, line): # check for closing paranthesis
                      break
            continue

        a = re.search(re_periodic, line)
        if a:
            print 'Reading periodic connectivity\n', line
            read_periodic(ifile)
            continue

        #print 'Line = ',line
        if any([re.search(st, line) for st in (re_parant, re_comment)]) or \
                                                             not line.strip():
            #print [re.search(st, line) for st in (re_parant, re_comment)]
            continue

        # Should not make it here
        #print 'Line = ',line
        #raise IOError('Something went wrong reading fluent mesh.')
def convert(fluentmesh):
    """Converts a fluent mesh to a mesh format that can be used by FEniCS.

         fluentmesh = fluent mesh (*.msh file)

    """
    ofilename = fluentmesh[:-4]
    ifile  = open(fluentmesh, "r")
    scan_fluent_mesh(ifile)
    ifile.close()
    print "Writing Files"
    with open('mesh.pkl','wb') as file:
        pickle.dump(neighbour,file,pickle.HIGHEST_PROTOCOL)
    with open('zone.pkl','wb') as file:
        pickle.dump(zone,file,pickle.HIGHEST_PROTOCOL)
    with open('bc2f.pkl','wb') as file:
        pickle.dump(boundary_cells2face,file,pickle.HIGHEST_PROTOCOL)
    with open('pfaces.pkl','wb') as file:
        pickle.dump(pfaces,file,pickle.HIGHEST_PROTOCOL)

    face_area()
    with open('facearea.pkl', 'wb') as file:
        pickle.dump(face_area_list, file, pickle.HIGHEST_PROTOCOL)
    print "Finished Writing Files"
    return 0


if __name__=='__main__':
    neighb=convert("CFD_CRN.cas")
