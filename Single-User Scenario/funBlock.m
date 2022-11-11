function fval = funBlock(X, X_all, ind_bs, clus)
global costq_t_block;
global data_t;
global ch;
global xi;
global V;
global cost_th;
global delay_edge;
global delay_bkb;
global clus_all;
global theta;

X_all(:,ind_bs) = X;
delay_ucn = data_t / log2(1 + calcu_sinr(clus_all, clus, ch));
fval = costq_t_block * (clus * (xi' * X_all)' - cost_th) +...
    V * (delay_ucn + theta * delay_edge + (1 - theta) * delay_bkb);
end