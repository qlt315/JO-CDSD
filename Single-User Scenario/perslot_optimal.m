function [clus, cach_policy] = perslot_optimal()
%% 穷举法寻找当前时隙每一个分簇下的最优存储策略，比较得到组合最优解。
global costq_t;
global req_type
global data_t;
global ch;
global xi;
global V;
global K;
global cost_th;
global N_agents;
global delay_edge;
global delay_bkb;
global size_se;
global comp_se;
global S_bs;
global C_bs;
global clus_all;
global typical_user;
global B;
global A;

opti_clus = zeros(1,N_agents);
opti_cach_policy = zeros(K,N_agents);
ind = 1;
min_val = inf;
while ind < 1024
    cur_clus = dec2bin(ind)-'0';
    if size(cur_clus,2) < N_agents
        less_len = N_agents - size(cur_clus,2);
        prex = zeros(1,less_len);
        cur_clus = [prex, cur_clus];
    end
    clus_all(typical_user,:) = zeros(1,N_agents);
    if sum(cur_clus) > B || max(sum(clus_all,1) + cur_clus) > A
        ind = ind + 1;
        continue;
    end
    %% 求解当前分簇下的优化问题
    delay_ucn = data_t / log2(1 + calcu_sinr(clus_all, cur_clus, ch));
    % 和instant的主要区别在于这个优化目标，为在cost约束下的delay
    func = @(X)delay_ucn + max(cur_clus .* (req_type' * X)) * (delay_edge - delay_bkb) + delay_bkb;
    X_init = zeros(K,N_agents);
    VLB = zeros(K,N_agents);
    VUB = ones(K,N_agents);
    options = optimoptions('fmincon','MaxIter',100000,'MaxFunEvals',100000);
    [cur_cach_policy, cur_val] = fmincon(func,...
        X_init,[],[],[],[],VLB,VUB,@(X) consPerslot(X, cur_clus), options);
    if min_val > cur_val
        min_val = cur_val;
        opti_clus = cur_clus;
        opti_cach_policy = cur_cach_policy;
    end
    ind = ind + 1;
end
clus = opti_clus;
cach_policy = opti_cach_policy;
end