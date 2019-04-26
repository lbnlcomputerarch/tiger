#***************************************************************************************************
#                        
# File Name:             qubo_gen.py 
# Application Version:   v0.1 
# Application Developer: Anastasiia Butko (LBNL)
#                        
# Software:              Task assIGnment mappER (TIGER)
# Author:                Anastasiia Butko (LBNL)
# Website:               https://github.com/lbnlcomputerarch/tiger
#                        
# The copyright information of this software can be found in the file COPYRIGHT. 
#
#*************************************************************************************************
from qubo_gen import *
from qmasm_gen import *
from sol_gen import *
from gen_support import *

import random
import math
import os
import sys
import ast
import re

import networkx as nx
import numpy as np
#import matplotlib.pyplot as plt
import collections
import fileinput

from random import randint
#from matplotlib import pylab as pl
from itertools import izip_longest, ifilter

def generate_qbsolv(TCG, SG, ARC, AppName, ProcNum, RowNum, ColNum, fileID, adjustedLevels, prevSol):
	# Write the qmasm file
	file = open(AppName + fileID + ".qubo", 'w')
	print("Generating single .qubo file...")
	file.write("p qubo 0 0 0 0\n")

	QubitsNum = len(SG)*ProcNum
	CouplersNum = 0
	# weight = [-2, 2]
	# strength = [-1, 1]

	penalty = 500
	diff = 100
	div = 10

	CostRange = []
	compCostRange = []
	commCostRange = []

	for i, task in enumerate(SG.nodes()):
		taskCost = SG.node[task]['weight']
		for j, core in enumerate(sorted(ARC.nodes())): #range(ProcNum):
			# Write costs per each node
			compCost = float(ARC.node[core]['weight'])
			compCostRange.append(taskCost * compCost)

	print("Computation cost range: [" + str(min(compCostRange)) + ", " + str(max(compCostRange)) + "]")
	if (max(compCostRange)-min(compCostRange) > 0):
		koef1 = abs(diff/(max(compCostRange)-min(compCostRange)))
		koef2 = abs((-diff*max(compCostRange) - 0.0*min(compCostRange))/(max(compCostRange)-min(compCostRange)))
	else:
		koef1 = 0 #abs(2/max(compCostRange))
		koef2 = diff #0

	TaskQubitDict = collections.OrderedDict()
	for i, task in enumerate(sorted(SG.nodes())):
		taskCost = SG.node[task]['weight']
		for j, core in enumerate(sorted(ARC.nodes())): #range(ProcNum):
			# Write costs per each node
			compCost = float(ARC.node[core]['weight'])
			TaskQubitDict.update({task: (i*ProcNum)})
	 		file.write(str(j + i*ProcNum) + " " + str(j + i*ProcNum) + " " + str( ((taskCost * compCost * koef1) - koef2)/div ) + "\n")

	# Analyze input edges from the outside
	input_edges = TCG.in_edges(sorted(SG.nodes()))
	missed_edges = 0
	extra_qubits = 0
	src_list = {}
	
	for link in input_edges:
		if link not in SG.edges() and prevSol:
			src = link[0]
			if src not in src_list:
				src_list.update({src: QubitsNum})
				file.write(str(QubitsNum) + " " + str(QubitsNum) + " " + str( -penalty ) + "\n")
				QubitsNum += 1
			missed_edges += 1
		elif link not in SG.edges():
			src = link[0]
			if src not in src_list:
				src_list.update({src: QubitsNum})
	
	extra_qubits = len(src_list)
	print("Number of missed edges: " + str(missed_edges) + " Number of extra qubits: " + str(extra_qubits))

	# Allocate couplers per each edge
	for idx, link in enumerate(SG.edges()):
		src = link[0]
		dst = link[1]
		edgeCost = SG.adj[src][dst]['weight']
		for i in sorted(ARC.nodes()):
			for j in sorted(ARC.nodes()):
				if i!=j:
					commCost = nx.dijkstra_path_length(ARC, source=i, target=j, weight='weight')
				else:
					commCost = 0.0
				commCostRange.append(edgeCost * commCost)

	for link in input_edges:
		if link not in SG.edges() and prevSol:
			src = link[0]
			dst = link[1]
			edgeCost = TCG.edge[src][dst]['weight']
			for i in sorted(ARC.nodes()):
				if prevSol[src] != i:
					commCost = nx.dijkstra_path_length(ARC, source=prevSol[src], target=i, weight='weight')
				else:
	 				commCost = 0.0
	 			commCostRange.append(edgeCost * commCost)	

	if commCostRange:
		print("Communication cost range: [" + str(min(commCostRange)) + ", " + str(max(commCostRange)) + "]")

		koef3 = abs(diff/(max(commCostRange)-min(commCostRange)))
		koef4 = abs((-diff*max(commCostRange) - 0.0*min(commCostRange))/(max(commCostRange)-min(commCostRange)))

		for idx, link in enumerate(SG.edges()):
			src = link[0]
			dst = link[1]
			edgeCost = SG.adj[src][dst]['weight']
			for i in sorted(ARC.nodes()):
				for j in sorted(ARC.nodes()):
					A = TaskQubitDict[src] + (i[0]*ColNum + i[1])
					B = TaskQubitDict[dst] + (j[0]*ColNum + j[1])
					if i!=j:
						commCost = nx.dijkstra_path_length(ARC, source=i, target=j, weight='weight')
					else:
						commCost = 0.0
					file.write(str(A) + " " + str(B) + " " + str( ((edgeCost * commCost * koef3) - koef4)/div) + "\n") #NOTE: koef1 & koef2 are used!
					CouplersNum += 1

		for link in input_edges:
			if link not in SG.edges() and prevSol:
				src = link[0]
				dst = link[1]
				edgeCost = TCG.edge[src][dst]['weight']
				for i in sorted(ARC.nodes()):
					A = src_list[src]
					B = TaskQubitDict[dst] + (i[0]*ColNum + i[1])

					if prevSol[src] == i:
						commCost = 0.0
					else:
						commCost = nx.dijkstra_path_length(ARC, source=prevSol[src], target=i, weight='weight')
					A, B = B, A
					file.write(str(A) + " " + str(B) + " " + str( ((edgeCost * commCost * koef3) - koef4)/div) + "\n") #NOTE: koef1 & koef2 are used!
					CouplersNum += 1

	# Allocate couplers per each level
	for lvl in adjustedLevels:
	 	lvlTasksNum = len(adjustedLevels[lvl])
		for i in adjustedLevels[lvl]:
	 		VerCopNum = vertical_couplings(TaskQubitDict.keys().index(i), penalty, file, ProcNum)
	 		CouplersNum += VerCopNum
	 	if lvlTasksNum > 1:
	 		for i, itemA in enumerate(adjustedLevels[lvl]):
	 			for j in range(i, len(adjustedLevels[lvl])):
	 				if i != j:
	 					HorCopNum = horizontal_couplings(TaskQubitDict[adjustedLevels[lvl][i]], TaskQubitDict[adjustedLevels[lvl][j]], penalty, file, ProcNum)
	 					CouplersNum += HorCopNum

	print("Number of qubits: "+ str(QubitsNum) + " Number of couplers: " + str(CouplersNum) + "\n")
	file.close()

	for line in fileinput.FileInput(AppName + fileID + ".qubo", inplace=True):
		if line == "p qubo 0 0 0 0\n":
			sys.stdout.write("p qubo 0 " + str(QubitsNum) + " " + str(QubitsNum) + " " + str(CouplersNum) + "\n")
		else:
			sys.stdout.write(line)

	fileinput.close()

	# if not missed_edges:
	# 	extra_qubits = 0

	return TaskQubitDict, extra_qubits