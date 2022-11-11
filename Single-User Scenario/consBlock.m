function [G, Geq] = consBlock(X, X_all, ind_bs, clus)
global size_se;
global comp_se;
global C_bs;
global S_bs;
global req_type;
global N_agents;
global theta;

X_all(:,ind_bs) = X;
G(1:N_agents) = (size_se' * X_all - clus .* S_bs')';
G(1+N_agents:N_agents*2) = (comp_se' * X_all - clus .* C_bs')';
Geq = max(clus .* (req_type' * X_all)) - theta;
end