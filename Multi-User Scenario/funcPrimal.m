function f = funcPrimal(X, clus_u)
global costq_t;
% global req_type
global data_t;
global ch;
global xi;
global V;
% global K;
global cost_th;
% global N_agents;
global clus_all;
global theta;
global B0_size;

delay_ucn = sum(data_t / log2(1 + calcu_sinr(clus_all, clus_u, ch))) / B0_size;
f = costq_t * (clus_u * (xi' * X)' - cost_th) + V / B0_size * delay_ucn + V * theta; 
end