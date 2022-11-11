function [opti_solu, opti_val] = solveMaster(fea_cut, infea_cut, init_clus_u)
global B;
global N_agents;
global clus_all;
global typical_user;
global A;
%% exhaustive search
min_solu = zeros(1,N_agents);
min_val = inf; % 初始化最优值
init_ind = 1;
% ensure that the intial point is feasible
while init_ind < 1024
    node_solu = dec2bin(init_ind)-'0';
    if size(node_solu,2) < N_agents
        less_len = N_agents - size(node_solu,2);
        prex = zeros(1,less_len);
        node_solu = [prex, node_solu];
    end
    clus_all(typical_user,:) = zeros(1,N_agents);
    if sum(node_solu) > B || max(sum(clus_all,1) + node_solu) > A
        init_ind = init_ind + 1;
        continue;
    end
    cur_val = max_bender_cut(fea_cut, infea_cut, node_solu);
    if min_val > cur_val
        min_val = cur_val;
        min_solu = node_solu;
    end
    init_ind = init_ind + 1;
end
%% gibbs
% [~, F_cur] = fmincon(@(d_0) d_0,0,[],[],[],[],0,[],@(d_0) bender_cut(d_0, fea_cut, infea_cut, feasible, clus_u)); 
cosi_max = 0.8;
cosi_min = 0.01;
cosi = cosi_max;
iter_max = 2000;
cosi_step = 2 * (cosi_max - cosi_min) / iter_max;
dev = 0.5;
tau = 1;
val = zeros(1,N_agents);
min_solu_gibbs = [];
min_val_gibbs = inf;
clus_all(typical_user,:) = zeros(1, N_agents); % clear the state of typical user
clus_u = init_clus_u;
F_cur = max_bender_cut(fea_cut, infea_cut, clus_u);
if F_cur < min_val_gibbs
    min_val_gibbs = F_cur;
    min_solu_gibbs = clus_u;
end
%% start the iteration
while tau <= iter_max
    m = randi(N_agents);
    clus_u_new = clus_u;
    clus_u_new(1,m) = 1 - clus_u(1,m); % get the new state
    if sum(clus_u_new) == 0 || sum(clus_u_new) > B || sum(clus_all(:,m))+clus_u_new(1,m)>A
        val(1,tau) = min_val_gibbs;
        tau = tau + 1;
        continue;
    end
    if sum(clus_u_new) == 0
        F_new = F_cur + dev; % 如果是全0的解，允许以一定概率接受，以便于下一次的状态转移
    else
        F_new = max_bender_cut(fea_cut, infea_cut, clus_u_new);
    end
    if F_new == inf
        val(1,tau) = min_val_gibbs;
        tau = tau + 1;
        continue;
    end
    if F_new < min_val_gibbs
        min_val_gibbs = F_new;
        min_solu_gibbs = clus_u_new;
    end
    % 状态转移
    yita = 1 / (1 + exp(min(dev, (F_new - F_cur)) / cosi_max));
    if rand() <= yita
        clus_u = clus_u_new;
        F_cur = F_new;
    end
%     if F_cur < min_val_gibbs
%         min_val_gibbs = F_cur;
%         min_solu_gibbs = clus_u;
%     end
    val(1,tau) = min_val_gibbs;
    tau = tau + 1;
%     if tau <= iter_max / 2
%         cosi = cosi - cosi_step;
%     end
end
opti_solu = min_solu_gibbs;
opti_val = min_val_gibbs;
% figure(); 
% plot(val);
% hold on;
% plot(min_val*ones(size(val)));
% title('convergence of gibbs');
% legend('gibbs','exhaustive search');
% xlabel('iterations');
% ylabel('target value');
% grid on;
end