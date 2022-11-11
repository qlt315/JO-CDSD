function X = uplink_optimal(clus)
global costq_t_opti;
global data_t;
global ch;
global xi;
global V;
global K;
global cost_th;
global N_agents;
global delay_edge;
global delay_bkb;
global clus_all;
global theta;

delay_ucn = data_t / log2(1 + calcu_sinr(clus_all, clus, ch));
func = @(X)costq_t_opti * (clus * (xi' * X)' - cost_th) +...
    V * (delay_ucn + theta * delay_edge + (1-theta)*delay_bkb);
X_init = zeros(K,N_agents);
VLB = zeros(K,N_agents);
VUB = ones(K,N_agents);
options = optimoptions('fmincon','MaxIter',100000,'MaxFunEvals',100000);
[X, ~] = fmincon(func,X_init,[],[],[],[],VLB,VUB,@(X)consPrimal(X, clus), options);
end