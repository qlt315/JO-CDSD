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
global delay_edge;
global delay_bkb
global clus_all;
global theta;

delay_ucn = data_t / log2(1 + calcu_sinr(clus_all, clus_u, ch));
f = costq_t * (clus_u * (xi' * X)' - cost_th) + V * (delay_ucn + theta * delay_edge + (1-theta) * delay_bkb); 
end