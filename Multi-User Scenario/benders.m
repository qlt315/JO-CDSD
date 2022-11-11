function [opti_clus_u, opti_cach_policy, flag] = benders()
%% global parameters
global K;
global N_agents;
global clus_all;
global typical_user;
global B;
global A;

%% parameters of Bender's decomposition
epsilon_gap = 1e-4; % stop criterion
iter_max = 100;
init_value = linspace(1,1023,1023);
init_ind = randi(1023);
feasible = 0;
flag = 1;
% ensure that the intial point is feasible
while feasible == 0
    clus_u = dec2bin(init_value(init_ind))-'0';
    init_value(init_ind) = [];
    if size(init_value,2)==0
        flag = 0;
        opti_clus_u = [];
        opti_cach_policy = [];
        return;
    end
    if size(clus_u,2) < N_agents
        less_len = N_agents - size(clus_u,2);
        prex = zeros(1,less_len);
        clus_u = [prex, clus_u];
    end
    clus_all(typical_user,:) = zeros(1,N_agents);
    if sum(clus_u) > B || max(sum(clus_all,1) + clus_u) > A
        init_ind = randi(size(init_value,2));
        continue;
    end
    [~, ~, feasible, ~] = solvePrimal(clus_u);
    init_ind = randi(size(init_value,2));
end
%% calcu an init clus_u solution
% fea_cut = [];
% infea_cut = [];
% new_cut{1,1} = init_primal;
% new_cut{1,2} = init_lambda;
% fea_cut = [fea_cut; new_cut];
% init_clus_u = clus_u;
% [clus_u, ~] = solveMaster(fea_cut, infea_cut, init_clus_u);
%% start benders with the init clus_u
LBD = -1 * inf;
UBD = inf;
fea_cut = [];
infea_cut = [];
upbound = [];
lowbound = [];
tau = 1;
while tau <= iter_max
    %% slove the primal problem
    [primal_solu, lambdaFea, feasible, opti_val_primal] = solvePrimal(clus_u);
    if feasible == 0  % add an infeasible cut
        [primal_solu, lambdaInfea] = solveInfeaPrimal(clus_u);
        new_cut{1,1} = primal_solu;
        new_cut{1,2} = lambdaInfea;
        infea_cut = [infea_cut; new_cut];
    else
        % stop criteria
        if opti_val_primal < inf && opti_val_primal - LBD < epsilon_gap && LBD < inf
            upbound = [upbound; opti_val_primal];
            lowbound = [lowbound; LBD];
            break;
        else % add a feasible cut
            new_cut{1,1} = primal_solu;
            new_cut{1,2} = lambdaFea;
            fea_cut = [fea_cut; new_cut];
            if UBD > opti_val_primal
                UBD = opti_val_primal;
                opti_clus_u = clus_u;
                opti_cach_policy = primal_solu;
            end
        end
    end        
    %% slove the master problem
    init_clus_u = clus_u;
    [clus_u, opti_val_master] = solveMaster(fea_cut, infea_cut, init_clus_u);
    % update LBD
    if LBD < opti_val_master && opti_val_master < inf
        LBD = opti_val_master;
    end
    upbound = [upbound; UBD];
    lowbound = [lowbound; LBD];
    % stop criteria
    if UBD - LBD < epsilon_gap && LBD < inf
        break;
    end
    tau = tau + 1;
end
% figure();
% plot(upbound);
% hold on;
% plot(lowbound);
% grid on;
% xlabel('iterations');
% ylabel('target value');
% legend('upper bound', 'lower bound');
end