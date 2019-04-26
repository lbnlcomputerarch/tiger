#***************************************************************************************************
#                        
# File Name:             gen_support.py 
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

def longest_path_length(G, src, dst):
	paths = nx.all_simple_paths(G, src, dst)
	check = list(nx.all_simple_paths(G, src, dst))
	if check:
		return len(max(paths, key = len))-1
	else:
		return 0

def level_groups(G, root, size):
	levels = collections.OrderedDict()
	levels.update({0: root})
	for taskID in G:
		if taskID not in root:
			levelID = 0
			for rootIdx in root:
				level = longest_path_length(G, rootIdx, taskID)
				levelID = level if level > levelID else levelID
			if (levelID in levels):
				levels[levelID].append(taskID)
			else:
				levels.update({levelID: [taskID]})
	return levels

def adjust_levels(G, levels, procNum):
	levelShift = 0
	adjustedLevels = collections.OrderedDict()
	print("\nAdjusting levels ...")
	for il, lvl in enumerate(levels):
		if len(levels[lvl]) > procNum:
			subLevelsNum = len(levels[lvl])/procNum
			remainderTasks = len(levels[lvl])%procNum
		
			print("Num of sublevels: " + str(subLevelsNum) + " RemainderTasks: " + str(remainderTasks))
			for idx in range(0, subLevelsNum):
				adjustedLevels.update({(lvl+levelShift): levels[lvl][(idx*procNum):(idx*procNum + procNum)]})
				levelShift += 1
			# if remainderTasks > 0:
			# 	adjustedLevels.update({(lvl+levelShift): levels[lvl][(subLevelsNum*procNum):(subLevelsNum*procNum + procNum)]})
			if remainderTasks > 0:
				if (il < len(levels)-1):
					if len(levels[lvl+1]) < procNum:
						bufA = []
						bufB = []
						freeSpace = procNum - len(levels[lvl+1])
						for subTask in range(0, remainderTasks):
							depFlag = 0
							for nextLevelTask in range(0, len(levels[lvl+1])):
								if (G.has_edge(levels[lvl][subLevelsNum*procNum+subTask], levels[lvl+1][nextLevelTask])):
									depFlag = 1
									break	
							if((depFlag != 1) and (len(bufA) < freeSpace)):
								bufA.append(levels[lvl][subLevelsNum*procNum+subTask])
							else:
								bufB.append(levels[lvl][subLevelsNum*procNum+subTask])					
						for taskID in bufA:
							levels[lvl+1].append(taskID)
						if bufB:
							adjustedLevels.update({(lvl+levelShift): bufB})
						else:
							levelShift -= 1
					else:
						adjustedLevels.update({(lvl+levelShift): levels[lvl][(subLevelsNum*procNum):(subLevelsNum*procNum + procNum)]})
				else:
					adjustedLevels.update({(lvl+levelShift): levels[lvl][(subLevelsNum*procNum):(subLevelsNum*procNum + procNum)]})
			else:
				levelShift -= 1
		else:
			adjustedLevels.update({(lvl+levelShift): levels[lvl]})
	print("\n")
	return adjustedLevels

def vertical_couplings(taskID, penalty, file, ProcNum):
	number = 0
	for i in range(0, ProcNum):
		for j in range(i, ProcNum):
			if i!=j:
				A = taskID*ProcNum + i
				B = taskID*ProcNum + j
				file.write(str(A) + " " + str(B) + " " + str(penalty) + "\n")
				number += 1
	return number

def horizontal_couplings(taskID_0, taskID_1, penalty, file, ProcNum):
	number = 0
	for p in range(0, ProcNum):
		A = taskID_0 + p
		B = taskID_1 + p
		file.write(str(A) + " " + str(B) + " " + str(penalty) + "\n")
		number += 1
	return number

def task_to_qubit(SG, ProcNum):
	TaskQubitDict = collections.OrderedDict()
	for i, task in enumerate(sorted(SG.nodes())):
		TaskQubitDict.update({task: (i*ProcNum)})
	return TaskQubitDict

def devide_graph(TCG, LpSG, adjustedLevels):
	SG = []
	SL = []

	chunks = [adjustedLevels.iteritems()]*LpSG
	g = (dict(ifilter(None, v)) for v in izip_longest(*chunks))
	SL = list(g)

	for sub_level in SL:
		SGlevels = []
		for key, value in sub_level.iteritems():
			SGlevels.extend(value)
		SG.append(TCG.subgraph(SGlevels))

	return SG, SL