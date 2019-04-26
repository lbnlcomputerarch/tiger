#***************************************************************************************************
#                        
# File Name:             tiger.py 
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

if __name__ == "__main__":

	TCG_file = sys.argv[1] 	#"example1.ram"
	ARC_file = sys.argv[2] 	#"mesh_2x2_homo.arc"

	print("Creating Architecture Graph ...")
	# Read the ARC file into ARC graph
	lines = open(ARC_file, 'r').read().split("\n")
	headerStart = 14
	headerEnd = 17

	# Read TCG file info
	# global Topology, ProcNum, RowNum, ColNum, AppName
	Topology = (' '.join(lines[headerStart].split())).split(" ")[1]
	ProcNum = int((' '.join(lines[headerStart+1].split())).split(" ")[1])
	RowNum = int((' '.join(lines[headerStart+2].split())).split(" ")[1])
	ColNum = int((' '.join(lines[headerStart+3].split())).split(" ")[1])
	

	ProcCost = (' '.join(lines[headerStart+5].split())).split(" ")[2].split(',')
	LinkCost = float((' '.join(lines[headerStart+6].split())).split(" ")[2])

	# Create the architecture graph
	ARC = nx.grid_2d_graph(RowNum, ColNum)
	for idx, nodeID in enumerate(sorted(ARC.nodes())):
		ARC.add_node(nodeID, weight=ProcCost[idx])
	for edgeID in ARC.edges():
		ARC.add_edge(edgeID[0], edgeID[1], weight=LinkCost)

	print("Creating Task Communication Graph ...")
	# Read the TCG file into TCG graph
	lines = open(TCG_file, 'r').read().split("\n")
	headerStart = 19
	headerEnd = headerStart+3

	# Read TCG file info
	AppName = lines[headerStart - 1 ].split(" ")[1]
	TaskNum = int(lines[headerStart].split(" ")[0])
	EdgeNum = int(lines[headerStart].split(" ")[1])

	TaskLines = lines[headerEnd:(headerEnd + TaskNum)]
	EdgeLines = lines[(headerEnd + TaskNum):]

	# Create the task graph
	TCG = nx.DiGraph()
	for i, line in enumerate(TaskLines):
		taskCost = int(line.split(" ")[3])
		TCG.add_node(i, weight=taskCost)

	for line in EdgeLines:
		src = int(line.split(" ")[1])
		dst = int(line.split(" ")[2])
		edgeCost = int(line.split(" ")[3])
		TCG.add_edge(src, dst, weight=edgeCost)

	root = []
	for taskID, degree in TCG.in_degree:
		if degree == 0:
			root.append(taskID)
	levels = level_groups(TCG, root, TaskNum)
	adjustedLevels = adjust_levels(TCG, levels, ProcNum)

	QubitsNum = TaskNum*ProcNum
	QubitsLim = 30

	NumSG = 0
	LpSG = 1

	opCode = sys.argv[3]

	if opCode == "i":
		SG, SL = devide_graph(TCG, LpSG, adjustedLevels)
		NumSG = len(SG)

		print("Number of qubits required: " + str(QubitsNum) + " Qubit limit: " + str(QubitsLim) + "\n")
		print("Number of sub-graphs required: " + str(NumSG))
		print("Number of adjusted levels: " + str(len(adjustedLevels)))
		print("Numer of levels per sub-graph: " + str(LpSG) + "\n")

	if opCode == "qb":
		prevSol = {}
		TaskQubitDict = []
		TaskQubitDict = generate_qbsolv(TCG, TCG, ARC, AppName, ProcNum, RowNum, ColNum, '', adjustedLevels, prevSol)
		print("Done.")

	if opCode == "qbs":
		dwAnswer = qbsolv_sol(AppName, '', 0)
		val = task_per_node_check(dwAnswer, ProcNum)
		if val[0]:
			combine_solutions(TCG, ARC, dwAnswer, ProcNum, RowNum, ColNum)
		else:
			combine_solutions(TCG, ARC, val[1], ProcNum, RowNum, ColNum)
			print("Solution is not valid.")

	if opCode == "dqb":
		prevSol = {}
		TaskQubitDict = []
		SG, SL = devide_graph(TCG, LpSG, adjustedLevels)
		NumSG = len(SG)
		for i in range(0, NumSG):
			TaskQubitDict.append(generate_qbsolv(TCG, SG[i], ARC, AppName, ProcNum, RowNum, ColNum, str(i), SL[i], prevSol)[0])

	if opCode == "dqbs":
		prevSol = {}
		TaskQubitDict = []
		SG, SL = devide_graph(TCG, LpSG, adjustedLevels)
		NumSG = len(SG)
		for i in range(0, NumSG):
			TaskQubitDict.append(generate_qbsolv(TCG, SG[i], ARC, AppName, ProcNum, RowNum, ColNum, str(i), SL[i], prevSol)[0])
		finalSolution = []
		for i in range(0, NumSG):
			optimalSolution = []
			optimalSolution.append(0)
			optimalSolution.append(0)
			optimalSolution.append(0)
			optimalSolution.append(qbsolv_sol(AppName, str(i), 0))
			finalSolution.append(optimalSolution)
		dwAnswer = get_dwAnswer(finalSolution, SG, TaskQubitDict, QubitsNum, ProcNum)
		combine_solutions(TCG, ARC, dwAnswer, ProcNum, RowNum, ColNum)
		if not task_per_node_check(dwAnswer, ProcNum)[0]:
			print("Solution is not valid")

	if opCode == "ddqb":
		it = int(sys.argv[4])
		SG, SL = devide_graph(TCG, LpSG, adjustedLevels)
		NumSG = len(SG)
		if it == 0:
			# Generate current qmasm 
			prevSol = {}
			TaskQubitDict = generate_qbsolv(TCG, SG[it], ARC, AppName, ProcNum, RowNum, ColNum, str(it), SL[it], prevSol)
		elif it == NumSG:
			# Restore previous solutions to get their placements
			prevSol = {}
			TaskQubitDictList = []
			for i in range(0, it):
				TaskQubitDictList.append(task_to_qubit(SG[i], ProcNum))
			prevSolList = restore_previous(SG, TaskQubitDictList, AppName, ProcNum, ARC)
			# Build previous solution sub-graph and solution dictionary
			SubQubitsNum = len(SG[it-1].nodes())*ProcNum
			TaskQubitDict, extra_qubits = generate_qbsolv(TCG, SG[it-1], ARC, AppName, ProcNum, RowNum, ColNum, str(it-1), SL[it-1], prevSol)
			solution = qbsolv_sol(AppName, str(it-1), extra_qubits)
			val = task_per_node_check(solution, ProcNum)
			record(solution, AppName)	
		else:
			# Restore previous solutionst o get their placements
			if it > 1:
				TaskQubitDictList = []
				for i in range(0, it):
					TaskQubitDictList.append(task_to_qubit(SG[i], ProcNum))
				prevSolList = restore_previous(SG, TaskQubitDictList, AppName, ProcNum, ARC)
			else:
				prevSolList = {}

			# Build previous solution sub-graph and solution dictionary
			prevSol = {}
			SubQubitsNum = len(SG[it-1].nodes())*ProcNum
			TaskQubitDict, extra_qubits = generate_qbsolv(TCG, SG[it-1], ARC, AppName, ProcNum, RowNum, ColNum, str(it-1), SL[it-1], prevSol)
			solution = qbsolv_sol(AppName, str(it-1), extra_qubits)
			task_per_node_check(solution, ProcNum)
			record(solution, AppName)

			# Restore previous solutions to get their placements
			TaskQubitDictList = []
			for i in range(0, it):
				TaskQubitDictList.append(task_to_qubit(SG[i], ProcNum))
			prevSolList = restore_previous(SG, TaskQubitDictList, AppName, ProcNum, ARC)

			# Generate current qubo 
			TaskQubitDict, extra_qubits = generate_qbsolv(TCG, SG[it], ARC, AppName, ProcNum, RowNum, ColNum, str(it), SL[it], prevSolList)

	if opCode == "ddqbs":
		dwAnswer = restore(AppName)
		val = task_per_node_check(dwAnswer, ProcNum)
		if not val[0]:
			print("Solution is not valid")
			combine_solutions(TCG, ARC, val[1], ProcNum, RowNum, ColNum)
		else:
			combine_solutions(TCG, ARC, dwAnswer, ProcNum, RowNum, ColNum)

	if opCode == "q":
		prevSol = {}
		TaskQubitDict = []
		TaskQubitDict = generate_qmasm(TCG, ARC, AppName, ProcNum, RowNum, ColNum, '', adjustedLevels, prevSol, TCG)

	if opCode == "qs":
		prevSol = {}
		TaskQubitDict = []
		TaskQubitDict = generate_qmasm(TCG, ARC, AppName, ProcNum, RowNum, ColNum, '', adjustedLevels, prevSol, TCG)[0]
		build_solution(TCG, TCG, ARC, AppName, ProcNum, RowNum, ColNum, QubitsNum, '', TaskQubitDict, prevSol)

	if opCode == "dq":
		prevSol = {}
		TaskQubitDict = []
		SG, SL = devide_graph(TCG, LpSG, adjustedLevels)
		NumSG = len(SG)
		for i in range(0, NumSG):
			TaskQubitDict.append(generate_qmasm(SG[i], ARC, AppName, ProcNum, RowNum, ColNum, str(i), SL[i], prevSol, TCG))

	if opCode == "dqs":
		prevSol = {}
		TaskQubitDict = []
		SG, SL = devide_graph(TCG, LpSG, adjustedLevels)
		NumSG = len(SG)
		for i in range(0, NumSG):
			TaskQubitDict.append(generate_qmasm(SG[i], ARC, AppName, ProcNum, RowNum, ColNum, str(i), SL[i], prevSol)[0], TCG)
		finalSolution = []
		for i in range(0, NumSG):
			SubQubitsNum = len(SG[i].nodes())*ProcNum
			finalSolution.append(build_solution(TCG, SG[i], ARC, AppName, ProcNum, RowNum, ColNum, SubQubitsNum, str(i), TaskQubitDict[i], prevSol))
		dwAnswer = get_dwAnswer(finalSolution, SG, TaskQubitDict, QubitsNum, ProcNum)
		combine_solutions(TCG, ARC, dwAnswer, ProcNum, RowNum, ColNum)
		valid_percentage(AppName)

	if opCode == "ddq":
		it = int(sys.argv[4])
		SG, SL = devide_graph(TCG, LpSG, adjustedLevels)
		NumSG = len(SG)
		if it == 0:
			# Generate current qmasm 
			prevSol = {}
			TaskQubitDict = generate_qmasm(SG[it], ARC, AppName, ProcNum, RowNum, ColNum, str(it), SL[it], prevSol, TCG)
		elif it == NumSG:
			# Restore previous solutions to get their placements
			prevSol = {}
			TaskQubitDictList = []
			for i in range(0, it):
				TaskQubitDictList.append(task_to_qubit(SG[i], ProcNum))
			prevSolList = restore_previous(SG, TaskQubitDictList, AppName, ProcNum, ARC)
			# Build previous solution sub-graph and solution dictionary
			SubQubitsNum = len(SG[it-1].nodes())*ProcNum
			TaskQubitDict, extra_qubits = generate_qmasm(SG[it-1], ARC, AppName, ProcNum, RowNum, ColNum, str(it-1), SL[it-1], prevSol, TCG)
			solution = build_solution(TCG, SG[it-1], ARC, AppName, ProcNum, RowNum, ColNum, SubQubitsNum + extra_qubits, str(it-1), TaskQubitDict, prevSolList)[3]
			record(solution, AppName)	
		else:
			# Restore previous solutionst o get their placements
			if it > 1:
				TaskQubitDictList = []
				for i in range(0, it):
					TaskQubitDictList.append(task_to_qubit(SG[i], ProcNum))
				prevSolList = restore_previous(SG, TaskQubitDictList, AppName, ProcNum, ARC)
			else:
				prevSolList = {}

			# Build previous solution sub-graph and solution dictionary
			prevSol = {}
			SubQubitsNum = len(SG[it-1].nodes())*ProcNum
			TaskQubitDict, extra_qubits = generate_qmasm(SG[it-1], ARC, AppName, ProcNum, RowNum, ColNum, str(it-1), SL[it-1], prevSol, TCG)
			solution = build_solution(TCG, SG[it-1], ARC, AppName, ProcNum, RowNum, ColNum, SubQubitsNum + extra_qubits, str(it-1), TaskQubitDict, prevSolList)[3]
			record(solution, AppName)

			# Restore previous solutions to get their placements
			TaskQubitDictList = []
			for i in range(0, it):
				TaskQubitDictList.append(task_to_qubit(SG[i], ProcNum))
			prevSolList = restore_previous(SG, TaskQubitDictList, AppName, ProcNum, ARC)

			# Generate current qmasm 
			TaskQubitDict = generate_qmasm(SG[it], ARC, AppName, ProcNum, RowNum, ColNum, str(it), SL[it], prevSolList, TCG)

	if opCode == "ddqs":
		dwAnswer = restore(AppName)
		combine_solutions(TCG, ARC, dwAnswer, ProcNum, RowNum, ColNum)
		valid_percentage(AppName)	