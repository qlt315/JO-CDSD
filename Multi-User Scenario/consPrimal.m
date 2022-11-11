function [G, Geq] = consPrimal(X, clus_u)
global size_se;
global comp_se;
global C_bs;
global S_bs;
global req_type;
global N_agents;
global theta;
global B0_size;
global K;
global delay_edge;
global delay_bkb;

G(1:N_agents) = (size_se' * X - clus_u .* S_bs')';
G(1+N_agents:N_agents*2) = (comp_se' * X - clus_u .* C_bs')';
delay_pro = 0;
for kk=1:K
    delay_pro = delay_pro + max(clus_u .* (req_type(kk) * X(kk,:))) * (delay_edge(kk) - delay_bkb(kk)) + delay_bkb(kk);
end
delay_pro = delay_pro / B0_size;
Geq = theta - delay_pro;
end