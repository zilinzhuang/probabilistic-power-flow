# probabilistic-power-flow
A probabilistic power-flow demo, with back/forward sweep algorithm, Monte Carlo random sampling, Latin hypercube sampling and modified algorithm for looped networks

This program, written in 2021, is very easy to understand and can be used by anyone, whether you are a beginner or a student interested in power flow/distributed power generation/power systems, or a researcher with professional expertise. It has the following functions:

1. Processing for PQ, PV, PI, and PQ(V) bus types.
Considering the different constant (controlled) quantities at different buses, the variables are updated via compensation, or through multiple iterations. For PV buses in particular, this paper provides a detailed derivation of a simplified method for calculating reactive power compensation.

2. Based on research into the operation and control methods of wind turbine generators and photovoltaic power generation systems, this paper establishes power-flow calculation models for various types of distributed generation (DG).

3. A detailed comparison is made between two sampling approaches used in simulation-based probabilistic power-flow analysis under uncertainty: Monte Carlo Simulation (MCS) and the Latin Hypercube Sampling (LHS) method.

4. With the backward/forward sweep (BFS) method as the core, this paper proposes several probabilistic power-flow algorithms for distribution networks with different types of DG integration, and develops the corresponding programs.

5. The classical BFS method cannot handle meshed (looped) networks; therefore, the algorithm is modified to address this limitation.
