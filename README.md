# Simulation Code of our paper in IEEE TMC: ''Joint Optimization of Base Station Clustering and Service Caching in User-Centric MEC''
**Abstract**: Edge service caching can effectively reduce the delay or bandwidth overhead for acquiring and initializing applications.
To address single-base station (BS) transmission limitation and serious edge effect in traditional cellular-based edge service caching
networks, in this paper, we proposed a novel user-centric edge service caching framework where each user is jointly provided with edge
caching and wireless transmission services by a specific BS cluster instead of a single BS. To minimize the long-term average delay
under the constraint of the caching cost, a mixed integer non-linear programming (MINLP) problem is formulated by jointly optimizing
the BS clustering and service caching decisions. To tackle the problem, we propose JO-CDSD, an efficiently joint optimization algorithm
based on Lyapunov optimization and generalized benders decomposition (GBD). In particular, the long-term optimization problem can
be transformed into a primal problem and a master problem in each time slot that is much simpler to solve. The near-optimal clustering
and caching strategy can be obtained through solving the primal and master problem alternately. Extensive simulations show that the
proposed joint optimization algorithm outperforms other algorithms and can effectively reduce the long-term delay by at most 93.75%
and caching cost by at most 53.12%.

Running Environment: Matlab R2021b

This paper is available on https://ieeexplore.ieee.org/document/10275092

## Single-User Scenario

**main_time_slot.m** ：Performance comparison between the proposed algorithm and several comparison algorithms at different time slots.

**main_cluster_size.m** : Performance comparison under different base station cluster sizes. 

**main_cost_th.m** : Performance comparison under different caching cost thresholds.

## Multi-User Scenario

**main_time_slot.m** ：Performance comparison between the proposed algorithm and several comparison algorithms at different time slots.

**main_cluster_size.m** : Performance comparison under different base station cluster sizes. 

**main_cost_th.m** : Performance comparison under different caching cost thresholds.


