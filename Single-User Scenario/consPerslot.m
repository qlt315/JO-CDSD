function [G, Geq] = consPerslot(X, clus_u)
global size_se;
global comp_se;
global C_bs;
global S_bs;
global N_agents;
global cost_th;
global xi;

G(1:N_agents) = (size_se' * X - clus_u .* S_bs')';
G(1+N_agents:N_agents*2) = (comp_se' * X - clus_u .* C_bs')';
G(N_agents*2+1) = clus_u * (xi' * X)' - cost_th;
Geq = [];
end