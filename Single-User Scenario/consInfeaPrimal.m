function [G, Geq] = consInfeaPrimal(X, clus_u)
global N_agents;
global size_se;
global comp_se;
global C_bs;
global S_bs;
global req_type;
global theta;

inter_f(1:N_agents,1) = size_se' * X - clus_u .* S_bs';
inter_f(N_agents+1:N_agents*2,1) = comp_se' * X - clus_u .* C_bs';
inter_f(N_agents*2+1,1) = theta - max(clus_u .* (req_type' * X));
alpha = max(inter_f);
G(1) = -alpha;
G(2:N_agents+1) = (size_se' * X - clus_u .* S_bs' - ones(1,N_agents) * alpha)';
G(N_agents+2:2*N_agents+1) = (comp_se' * X - clus_u .* C_bs' - ones(1,N_agents) * alpha)';
Geq = max(clus_u .* (req_type' * X)) + alpha - theta;
end