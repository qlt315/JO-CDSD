function max_cut = max_bender_cut(fea_cut, infea_cut, clus_u)
global costq_t;
global cost_th;
global data_t;
global ch;
global theta;
global N_agents;
global xi;
global V;
global clus_all;
global B0_size;

delay_ucn = sum(data_t) / (B0_size * log2(1 + calcu_sinr(clus_all, clus_u, ch)));
% feasible cut
n = size(fea_cut,1);
G = zeros(n,1);
for i = 1:n
    cach_policy = fea_cut{i,1};
    lambda = fea_cut{i,2};
    [g, geq] = consPrimal(cach_policy, clus_u);
    G(i) = costq_t * (clus_u * (xi' * cach_policy)' - cost_th) +...
        V / B0_size * delay_ucn + V * theta + [g geq] * lambda;
end
max_cut = max(G);
% infeasible cut
n = size(infea_cut,1);
Gin = zeros(n,1);
for i = 1:n
    cach_policy = infea_cut{i,1};
    lambda = infea_cut{i,2};
    [gin, geqin] = consInfeaPrimal(cach_policy, clus_u);
    Gin(i) = [gin(2:N_agents*2+1) geqin] * lambda;
end
if max(Gin) > 0  % L~ should be less than 0
    max_cut = inf;
end
end