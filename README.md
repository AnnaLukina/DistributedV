# DistributedV
Distributed Control of a V-Formation

SAC
RESEARCH-ARTICLE
Distributed adaptive-neighborhood control for stochastic reachability in multi-agent systems

Publication:
SAC '19: Proceedings of the 34th ACM/SIGAPP Symposium on Applied ComputingApril 2019Pages 914–921https://doi.org/10.1145/3297280.3297370
0
 
ABSTRACT
We present DAMPC, a distributed, adaptive-horizon and adaptive-neighborhood algorithm for solving the stochastic reachability problem in multi-agent systems, in particular, flocking modeled as a Markov decision process. At each time step, every agent first calls a centralized, adaptive-horizon model-predictive control (AMPC) algorithm to obtain an optimal solution for its local neighborhood. Second, the agents derive the flock-wide optimal solution through a sequence of consensus rounds. Third, the neighborhood is adaptively resized using a flock-wide cost-based Lyapunov function. This way DAMPC improves efficiency without compromising convergence. We evaluate DAMPC's performance using statistical model checking. Our results demonstrate that, compared to AMPC, DAMPC achieves considerable speed-up (two-fold in some cases) with only a slightly lower rate of convergence. The smaller average neighborhood size and lookahead horizon demonstrate the benefits of the DAMPC approach for stochastic reachability problems involving any controllable multi-agent system that possesses a cost function.
