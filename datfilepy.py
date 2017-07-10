#!/usr/bin python
import sys
import re
from struct import*
import binascii
from copy import copy
from numpy import zeros, array, sign, cross, dot, ones, arctan, sin, cos, pi, mod, sqrt
import pickle


print "Reading ANSYS Fluent data file"

# Use regular expressions to identify sections and tokens found in a fluent file

re_data     = re.compile(r"\(314(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\)(\s*)(\(|)")
re_data_std = re.compile(r"\(300(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\)(\s*)(\(|)")
re_parant   = re.compile(r"(\s*\)(\s*)|\s*\)\)(\s*)|\s*\(\s*)")
re_comment  = re.compile(r"\(0\s.*")
re_cfn      = re.compile(r"\(cell-function-list")
# The fluent mesh (the .msh file) is basically stored as a list of nodes, and then a
# list of faces for each zone of the mesh, the interior and the boundaries.

# Declare som maps that will be built when reading in the lists of nodes and faces:
Data={}                     #{zone id:{var type:[...],var2:[..]},z id 2:{..}}
Data_std={}                 #{zone id:{var type:[...],var2:[..]},z id 2:{..}}
# Information about connectivity and boundaries
boundary_cells = {}         # List of cells attached to a boundary facet. Key is zone id

# Some global values
num_cells = {}              # Total number of cells in different zones
zones = {}                  # zone information
zone_number_of_faces = {}   # number of faces for each zone
header=[]
def read_header(line):
    ln=line.split('(1')
    count=1
    for i in ln:
        str = i.split('"')
        if len(str)>2:
            header.append(str[3])
            #print str[3]

def read_data(type_id, zone_id, size, n_time_levels, n_phases, Nmin, Nmax, ifile):
    """Read data"""
    line = ifile.readline()
    readline = False
    #if re.search(re_parant, line): # check for initial paranthesis
        #readline = True

    ls = []
    ls_sub=[]
    #print Nmax

    for i in range(Nmin, Nmax + 1):
        if i==Nmax:
            print "last line"
        if readline:
            line = ifile.readline()
        readline = True
        string=''
        if re.search(re_parant,line):
            for x in line:
                if x!=')' and x!='(':
                    string=string+x
                if x==')' or x=='(':
                    print line
                    print i
            line=string
        if len(line)>0:
            #print "line=",line
            ln = line.split()
            ls_sub=[float(x) for x in ln]
            ls.append(ls_sub[0])


    if zone_id in Data:
        Data[zone_id][type_id]=ls
    else:
        Data[zone_id]={type_id:ls}

    print "read data complete"

def read_data_std(type_id, zone_id, size, n_time_levels, n_phases, Nmin, Nmax, ifile):
    """Read data"""
    line = ifile.readline()
    readline = False
    #if re.search(re_parant, line): # check for initial paranthesis
        #readline = True

    ls = []
    ls_sub=[]
    #print Nmax

    for i in range(Nmin, Nmax + 1):
        if i==Nmax:
            print "last line"
        if readline:
            line = ifile.readline()
        readline = True
        string=''
        if re.search(re_parant,line):
            for x in line:
                if x!=')' and x!='(':
                    string=string+x
                if x==')' or x=='(':
                    print line
                    print i
            line=string
        if len(line)>0:
            #print "line=",line
            ln = line.split()
            ls_sub=[float(x) for x in ln]
            ls.append(ls_sub[0])


    if zone_id in Data_std:
        Data_std[zone_id][type_id]=ls
    else:
        Data_std[zone_id]={type_id:ls}

    print "read standard data complete"



def scan_dat(ifile):
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
        a = re.search(re_cfn, line)
        if a:
            print line
            read_header(line)
            sample=line
        a = re.search(re_data, line)
        if a:
            print 'Reading faces ', line
            #print line
            type_id, zone_id, size, n_time_levels, n_phases, first_id, last_id = \
                 int(a.group(2)), int(a.group(3)), int(a.group(4)), \
                 int(a.group(5)), int(a.group(6)), int(a.group(7)),\
                 int(a.group(8),)
            #print "type=",a.group(2)
            #print "zone id=",a.group(3)
            #print "Nmax=",int(a.group(8))
            read_data(type_id, zone_id, size, n_time_levels, n_phases, first_id, last_id, ifile)
            zone_number_of_faces[zone_id] = last_id - first_id + 1
            continue

        a = re.search(re_data_std, line)
        if a:
            print 'Reading faces ', line
            # print line
            type_id, zone_id, size, n_time_levels, n_phases, first_id, last_id = \
                int(a.group(2)), int(a.group(3)), int(a.group(4)), \
                int(a.group(5)), int(a.group(6)), int(a.group(7)), \
                int(a.group(8), )
            # print "type=",a.group(2)
            # print "zone id=",a.group(3)
            # print "Nmax=",int(a.group(8))
            read_data_std(type_id, zone_id, size, n_time_levels, n_phases, first_id, last_id, ifile)
            zone_number_of_faces[zone_id] = last_id - first_id + 1
            continue


        #print 'Line = ',line
        if any([re.search(st, line) for st in (re_parant, re_comment)]) or \
                                                             not line.strip():
            #print [re.search(st, line) for st in (re_parant, re_comment)]
            continue

        # Should not make it here
        #print 'Line = ',line
        #raise IOError('Something went wrong reading fluent mesh.')
def readdata(datafile):
    """Converts a fluent mesh to a mesh format that can be used by FEniCS.

         fluentmesh = fluent mesh (*.msh file)

    """
    ofilename = datafile[:-4]
    ifile  = open(datafile, "r")
    scan_dat(ifile)
    ifile.close()
    header.append('Static Temperature')
    #xfile Static temperature id=3
    for zone in Data:
        # xfile Static temperature id=3
        if 3 in Data_std[zone]:
            Data[zone][header.index('Static Temperature')+1]=Data_std[zone][3]
        #xfile Enthalpy id=4
        if 4 in Data_std[zone]:
            Data[zone][header.index('Enthalpy')+1]=Data_std[zone][4]
    with open('data.pkl', 'wb') as file:
        pickle.dump(Data, file, pickle.HIGHEST_PROTOCOL)
    with open('data_std.pkl', 'wb') as file:
        pickle.dump(Data_std, file, pickle.HIGHEST_PROTOCOL)
    with open('header.pkl', 'wb') as file:
        pickle.dump(header, file, pickle.HIGHEST_PROTOCOL)
    return 0


if __name__=='__main__':
    readdata("CFD_CRN.dat")
    for i in Data:
        print "zone ",i
        for j in Data[i]:
            print j,
        print "\n"

