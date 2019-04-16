TIGER: Topology-aware task assIGnment mappER
========

## Abstract ##
Optimal mapping of a parallel code’s communication graph is increasingly
important as both system size and heterogeneity increase.
However, the topology-aware task assignment problem is an NPcomplete
graph isomorphism problem. Existing task scheduling
approaches are either heuristic or based on physical optimization
algorithms, providing different speed and solution quality tradeoffs.
Ising machines such as quantum and digital annealers have
recently become available offering an alternative hardware solution
to solve certain types of optimization problems. We propose
an algorithm that allows expressing the problem for such machines
and a domain specific partition strategy that enables to solve larger
scale problems. TIGER - topology-aware task assignment mapper
tool - implements the proposed algorithm and automatically integrates
task- communication graph and an architecture graph into
the quantum software environment. We use D-Wave’s quantum
annealer to demonstrate the solving algorithm and evaluate the
proposed tool flow in terms of performance, partition efficiency
and solution quality. Results show significant speed-up of the tool
flow and reliable solution quality while using TIGER together with
the proposed partition.

---
## Copyright ##

*Topology-aware task-assIGnment mappER (TIGER)* Copyright (c) 2019, The
Regents of the University of California, through Lawrence Berkeley National
Laboratory (subject to receipt of any required approvals from the U.S.
Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative
works, and perform publicly and display publicly, and to permit other to do
so.
