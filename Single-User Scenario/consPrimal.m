function [G, Geq] = consPrimal(X, clus_u)
global size_se;
global comp_se;
global C_bs;
global S_bs;
global req_type;
global N_agents;
global theta;

G(1:N_agents) = (size_se' * X - clus_u .* S_bs')';
G(1+N_agents:N_agents*2) = (comp_se' * X - clus_u .* C_bs')';
Geq = max(clus_u .* (req_type' * X)) - theta;
end