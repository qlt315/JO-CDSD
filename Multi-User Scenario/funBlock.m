function fval = funBlock(X, X_all, ind_bs, clus, optimal_theta)
global costq_t_block;
global data_t;
global ch;
global xi;
global V;
global cost_th;
global clus_all;
global B0_size;

theta = optimal_theta;
X_all(:,ind_bs) = X;
delay_ucn = sum(data_t / log2(1 + calcu_sinr(clus_all, clus, ch))) / B0_size;
fval = costq_t_block * (clus * (xi' * X_all)' - cost_th) +...
    V / B0_size * delay_ucn + V * theta; 
end