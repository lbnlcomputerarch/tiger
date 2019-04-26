#***************************************************************************************************
#                        
# File Name:             sol_gen.py 
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
import collections

from random import randint
from itertools import izip_longest, ifilter

def task_per_node_check(solution, ProcNum):
	valid = 1
	valid_answer = 1
	for i in range(0, len(solution)/ProcNum):
		valid = 1 if sum(solution[(i*ProcNum):(i*ProcNum + ProcNum)]) == 1 else 0
		#print(solution[(i*ProcNum):(i*ProcNum + ProcNum)], valid)
		if valid == 0:
			valid_answer = 0
			solution[(i*ProcNum):(i*ProcNum + ProcNum)] = one_hot(randint(0, ProcNum-1), ProcNum)
	return valid_answer, solution

def one_hot(i, ProcNum):
    a = np.zeros(ProcNum, 'uint8')
    a[i] = 1
    return a


def build_solution(TCG, SG, ARC, AppName, ProcNum, RowNum, ColNum, QubitsNum, fileID, TaskQubitDict, prevSol):
	# Read the ARC file into ARC graph
	lines = open(AppName + fileID + ".sol", 'r').read().split("\n")
	headerEnd = 4

	sol_lines = open(AppName + fileID + ".sn", 'r').read().split("\n")
	numberOfSolutions = int(sol_lines[1])
	dwAnswerList = []
	
	for i in range(0, numberOfSolutions):
		start = i * (headerEnd + QubitsNum + 1) + headerEnd
		solutionLines = lines[start:start + QubitsNum]
		dwAnswer = [None]*len(SG.nodes())*ProcNum
		for i, line in enumerate(solutionLines):
			idx = int((' '.join(line.split())).split(" ")[0])
			value = int(ast.literal_eval((' '.join(line.split())).split(" ")[2]))
			if idx < len(dwAnswer):
				dwAnswer[idx] = value

		if task_per_node_check(dwAnswer, ProcNum)[0]:
			dwAnswerList.append(dwAnswer)

	optimalSolution = []
	for dwAnswer in dwAnswerList:
		compSolution = 0
		for i, task in enumerate(SG.nodes()):
			taskCost = SG.node[task]['weight']
			for j, core in enumerate(sorted(ARC.nodes())):
				# Write costs per each node
				compCost = float(ARC.node[core]['weight'])
				compSolution += dwAnswer[j+i*ProcNum]*taskCost*compCost

		commSolution = 0
		for idx, link in enumerate(SG.edges()):
			src = link[0]
			dst = link[1]
			edgeCost = SG.edge[src][dst]['weight']
			for i in sorted(ARC.nodes()):
				for j in sorted(ARC.nodes()):
					if i!=j:
						A = TaskQubitDict[src] + (i[0]*ColNum + i[1])
						B = TaskQubitDict[dst] + (j[0]*ColNum + j[1])
						commCost = nx.dijkstra_path_length(ARC, source=i, target=j, weight='weight')
						commSolution += dwAnswer[A]*dwAnswer[B]*edgeCost*commCost

		input_edges = TCG.in_edges(sorted(SG.nodes()))

		extraCommCost = 0
		for link in input_edges:
			if link not in SG.edges() and prevSol:
				src = link[0]
				dst = link[1]
				edgeCost = TCG.adj[src][dst]['weight']
				for i in sorted(ARC.nodes()):
					B = TaskQubitDict[dst] + (i[0]*ColNum + i[1])
					if prevSol[src] != i:
						commCost = nx.dijkstra_path_length(ARC, source=prevSol[src], target=i, weight='weight')
					else:
						commCost = 0.0
					extraCommCost += dwAnswer[B]*edgeCost*commCost

		commSolution += extraCommCost

		solution = compSolution + commSolution
		if not optimalSolution:
			optimalSolution.append(solution)
			optimalSolution.append(compSolution)
			optimalSolution.append(commSolution)
			optimalSolution.append(dwAnswer)
		else:
			if optimalSolution[0] > solution:
				optimalSolution[0] = solution
				optimalSolution[1] = compSolution
				optimalSolution[2] = commSolution
				optimalSolution[3] = dwAnswer

	# File for valid solution number tracking
	file = open(AppName + ".vs", 'ab')	
	if not dwAnswerList:
		print("No valid solution found \n")
		file.write("0 ")
		file.close()
		# Generate random solution
		dwAnswer = [] #[None]*len(TCG.nodes())*ProcNum
		for i in range(0, len(SG.nodes())):
			dwAnswer.extend(one_hot(randint(0, ProcNum-1), ProcNum))
			optimalSolution.append(0)
			optimalSolution.append(0)
			optimalSolution.append(0)
			optimalSolution.append(dwAnswer)
	else:
		print("Total number of valid solutions: " + str(len(dwAnswerList)))
		print(optimalSolution[3])
		print("Computation cost: " + str(optimalSolution[1]) + "\n" + "Communication cost: " + str(optimalSolution[2]) + "\n" + "Total cost: " + str(solution) + "\n")
		file.write("1 ")
		file.close()
	return optimalSolution

def get_dwAnswer(solutions, SG, TaskQubitDict, QubitsNum, ProcNum):
	dwAnswer = [None]*QubitsNum
	NumSG = len(SG)
	for i in range(0, NumSG):
		for task in SG[i].nodes():
			subAnswer = solutions[i][3][TaskQubitDict[i][task]:TaskQubitDict[i][task] + ProcNum]
			dwAnswer[task*ProcNum:task*ProcNum + ProcNum] = subAnswer
	
	return dwAnswer

def qbsolv_sol(AppName, fileID, extraQubits):
	lines = open(AppName + fileID + ".sol", 'r').read().split("\n")
	dwAnswer = [int(n) for n in list(lines[1])]
	if extraQubits:
		end = len(dwAnswer)-extraQubits
		del dwAnswer[end:]
	return dwAnswer

def combine_solutions(TCG, ARC, dwAnswer, ProcNum, RowNum, ColNum):
	CompCostTCG = 0
	for i, task in enumerate(TCG.nodes()):
		taskCost = TCG.node[task]['weight']
		for j, core in enumerate(sorted(ARC.nodes())):
			# Write costs per each node
			compCost = float(ARC.node[core]['weight'])
			CompCostTCG += dwAnswer[j+i*ProcNum]*taskCost*compCost

	CommConstTCG = 0
	for idx, link in enumerate(TCG.edges()):
		src = link[0]
		dst = link[1]
		edgeCost = TCG.adj[src][dst]['weight']
		for i in sorted(ARC.nodes()):
			for j in sorted(ARC.nodes()):
				if i!=j:
					A = src*ProcNum + (i[0]*ColNum + i[1])
					B = dst*ProcNum + (j[0]*ColNum + j[1])
					commCost = nx.dijkstra_path_length(ARC, source=i, target=j, weight='weight')
					CommConstTCG += dwAnswer[A]*dwAnswer[B]*edgeCost*commCost

	finalSolutionTCG = CompCostTCG + CommConstTCG

	print(dwAnswer)
	print("\n")
	print("TCG combined solution: " + str(finalSolutionTCG) + " Comp/Comm: " + str(CompCostTCG) + "/" + str(CommConstTCG))

def record(solution, AppName):
	file = open(AppName + ".os", 'ab')
	file.writelines( "%s " % item for item in solution)
	file.write("\n")
	file.close()	

def restore_previous(SG, TaskQubitDict, AppName, ProcNum, ARC):
	file = open(AppName + ".os", 'r')
	lines = file.read().split("\n")
	prevSol = {}
	for idx, line in enumerate(lines):
		dwAnswer = [int(n) for n in line.split()]
		if dwAnswer:
			prevSol.update(previous_solution(dwAnswer, SG[idx], TaskQubitDict[idx], ProcNum, ARC))
	file.close()
	return prevSol

def previous_solution(solution, SG, TaskQubitDict, ProcNum, ARC):
	prevSol = {}
	for task in SG.nodes():
		subAnswer = solution[TaskQubitDict[task]: TaskQubitDict[task] + ProcNum]
		pos = ARC.nodes()[randint(0, ProcNum-1)] #workaround, has to be fixed
		for i, core in enumerate(sorted(ARC.nodes())):
			if subAnswer[i]:
				pos = core		
		prevSol.update({task: pos})
	return prevSol

def restore(AppName):
	file = open(AppName + ".os", 'r')
	lines = file.read().split("\n")
	finalSolution = []
 	for idx, line in enumerate(lines):
 		dwAnswer = [int(n) for n in line.split()]
 		if dwAnswer:
 			finalSolution.extend(dwAnswer)
 	file.close()
 	return finalSolution

def valid_percentage(AppName):
	file = open(AppName + ".vs", 'r')
	line = file.read()
	vsList = [int(n) for n in line.split()]
	vs = sum(vsList)
	print("Percentage of valid solutions found: " + str(vs*100.0/len(vsList)) + "%")
	file.close()
