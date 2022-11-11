function f = funcInfeaPrimal(X, clus_u)
global K;
global N_agents;
global size_se;
global comp_se;
global C_bs;
global S_bs;
global req_type;
global theta;
global B0_size;
global delay_edge;
global delay_bkb;

inter_f(1:N_agents,1) = size_se' * X - clus_u .* S_bs';
inter_f(N_agents+1:N_agents*2,1) = comp_se' * X - clus_u .* C_bs';
delay_pro = 0;
for kk=1:K
    delay_pro = delay_pro + max(clus_u .* (req_type(kk) * X(kk,:))) * (delay_edge(kk) - delay_bkb(kk)) + delay_bkb(kk);
end
delay_pro = delay_pro / B0_size;
inter_f(N_agents*2+1,1) = delay_pro - theta;
f = max(inter_f);
end