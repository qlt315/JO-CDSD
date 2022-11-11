clc;
clear;
%% global variables
global version;
% global H_t_single;
% global Q_t;
% global cost_th;
% 时隙不变全局变量
global V;
global xi;
global size_se;
global comp_se;
global C_bs;
global S_bs;
global K;
global N_agents;
global B;
global A;
global P;
global U;
global typical_user;
global cost_th;
% 时隙变化全局变量
global ch;
global clus_all;
global costq_t;
global costq_t_instant;
global costq_t_opti;
global costq_t_block;
global req_type;
global req_type_ind;
global delay_edge;
global delay_bkb;
global data_t;
global workld_t;
global theta;
% global R;
% rng(5,'twister');
V = 5;
%% parameters of BSs
% 当基站被分配给某个用户的基站簇,该基站的所有天线都会接收用户信号，不以天线区分
B = 3; % 基站簇的大小约束
B_max = 5;
A = 3; % 每个基站天线个数
P = 1; % user's transmit power
N_agents = 10; %number of BSs, i.e., M
% C_bs = randi([2,5],N_agents,1);
% S_bs = randi([2,5],N_agents,1);%computing and storage of BSs
C_bs = 3*ones(N_agents,1);
S_bs = 3*ones(N_agents,1);%computing and storage of BSs
U = 4; %number of users
typical_user = randi(U); % 指定typical用户
% sinr_th = 1;
%% parameters of task
K = 6; %number of service types
% size_se = randi([2,5],K,1); %size vector of services
% comp_se = randi([2,5],K,1); %computing requirement of services
% cost_coef = rand(K,1); %cost coefficient of caching services
size_se = 3 * ones(K,1);
comp_se = 3 * ones(K,1);
data_se = linspace(1,K,K);
workld_se = linspace(0.1,0.1*K,K);
cost_coef = linspace(0.1,0.1*K,K);
R = 0.05; %data rate of backbone
xi = cost_coef' .* size_se;
cost_th = 2;

T = 30;
t_step = 5;
%% 初始化
time = datetime;
version = [num2str(time.Year) num2str(time.Month) num2str(time.Day) num2str(time.Hour)];
% new_folder = sprintf('%s', ['output\' version 'gibbs_conver']);
% mkdir(new_folder);
% new_folder = sprintf('%s', ['output\' version 'admm_conver']);
% mkdir(new_folder);

% 初始化
% proposed
theta = 1;
C_t = zeros(1,N_agents);
X_t = zeros(K, N_agents);
costq_last = 0.1; 
costq_t = 0.1;
delay_T = zeros(2,T);
delay_uplink_T = zeros(2,T);
cach_cost_T = zeros(2,T);
% comparison:上行链路时延最优基站簇划分
C_opti_t = zeros(1,N_agents);
X_t_opti = zeros(K, N_agents);
costq_last_opti = 0.1; 
costq_t_opti = 0.1;
delay_T_opti = zeros(2,T);
delay_uplink_T_opti = zeros(2,T);
cach_cost_T_opti = zeros(2,T);
% comparison: instant optimal
C_t_instant = zeros(1,N_agents);
X_t_instant = zeros(K, N_agents);
costq_last_instant = 0.1; 
costq_t_instant = 0.1;
delay_T_instant = zeros(2,T);
delay_uplink_T_instant = zeros(2,T);
cach_cost_T_instant = zeros(2,T);
% comparison: block descent
C_t_block = zeros(1,N_agents);
X_t_block = zeros(K, N_agents);
costq_last_block = 0.1; 
costq_t_block = 0.1;
delay_T_block = zeros(2,T);
delay_uplink_T_block = zeros(2,T);
cach_cost_T_block = zeros(2,T);

% 多时隙的状态存储
% data_T = rand(1,T)*10;
% workld_T = rand(1,T);
req_type_T = gener_req(T, K, 0.5);
ch_T = zeros(N_agents*A, U, T);
clus_all_T = zeros(U, N_agents, T);
S_bs_T = zeros(N_agents,T);
C_bs_T = zeros(N_agents,T);
for t = 1:T
    ch = sqrt(1/2) * (randn(N_agents*A, U) + 1i * randn(N_agents*A, U));
    ch_T(:,:,t) = ch;
    S_bs_T(:,t) = round(rand(N_agents,1)) * 3;
    C_bs_T(:,t) = round(rand(N_agents,1)) * 3;
end
for t = 1:T
    [clus_all_T(:,:,t),~] = BScluster(sum(ch_T(:,:,t),3));
end
% marker_list = ['o' '*' 'd' 'v' '+'];
% figure(1); %delay
% figure(2); %uplink
% figure(3); %cost
delay_B = zeros(1,B_max);
delay_uplink_B = zeros(1,B_max);
cach_cost_B = zeros(1,B_max);
delay_B_opti = zeros(1,B_max);
delay_uplink_B_opti = zeros(1,B_max);
cach_cost_B_opti = zeros(1,B_max);
delay_B_instant = zeros(1,B_max);
delay_uplink_B_instant = zeros(1,B_max);
cach_cost_B_instant = zeros(1,B_max);
delay_B_block = zeros(1,B_max);
delay_uplink_B_block = zeros(1,B_max);
cach_cost_B_block = zeros(1,B_max);
for B = 1:B_max% begin long term optimization 
%     marker = marker_list(B);
    % 初始化
    % proposed
    C_t = zeros(1,N_agents);
    X_t = zeros(K, N_agents);
    costq_last = 0.1; 
    costq_t = 0.1;
    delay_T = zeros(2,T);
    delay_uplink_T = zeros(2,T);
    cach_cost_T = zeros(2,T);
    % comparison
    C_opti_t = zeros(1,N_agents);
    X_t_opti = zeros(K, N_agents);
    costq_last_opti = 0.1; 
    costq_t_opti = 0.1;
    delay_T_opti = zeros(2,T);
    delay_uplink_T_opti = zeros(2,T);
    cach_cost_T_opti = zeros(2,T);
    
    C_t_instant = zeros(1,N_agents);
    X_t_instant = zeros(K, N_agents);
    costq_last_instant = 0.1; 
    costq_t_instant = 0.1;
    delay_T_instant = zeros(2,T);
    delay_uplink_T_instant = zeros(2,T);
    cach_cost_T_instant = zeros(2,T);
    
    C_t_block = zeros(1,N_agents);
    X_t_block = zeros(K, N_agents);
    costq_last_block = 0.1; 
    costq_t_block = 0.1;
    delay_T_block = zeros(2,T);
    delay_uplink_T_block = zeros(2,T);
    cach_cost_T_block = zeros(2,T);
    for t = 1:T
        %% 当前states (所有方案统一使用的)
        req_type_ind = req_type_T(1,t);
        req_type = zeros(K,1);
        req_type(req_type_ind) = 1;
        data_t = data_se(1,req_type_ind);
        workld_t = workld_se(1,req_type_ind);
        S_bs = S_bs_T(:,t);
        C_bs = C_bs_T(:,t);
        % 计算任务在边缘处理需要的时间
        delay_edge = data_t * workld_t / comp_se(req_type_ind);
        delay_bkb = 10 * delay_edge;
        % 生成当前时隙信道系数，大小为MA*U，服从瑞利分布
        ch = sum(ch_T(:,:,t),3);
        % generate the clustering state of other users (independent from the typical user)
        clus_all = clus_all_T(:,:,t);
        %% proposed algorithm
        a_t = xi' * X_t * (C_t');
        costq_t = max(costq_last + a_t - cost_th, 0.1);
        % 当前时隙变量求解
        [C_t, X_t, feasible_flag] = benders()
        while feasible_flag == 0  % 问题不可解，放宽约束
            [C_t, X_t, feasible_flag] = benders()
        end
        disp('req_type_ind = ')
        disp(req_type_ind);
        % 计算上行链路时延
        delay_uplink = data_t / log2(1 + calcu_sinr(clus_all, C_t, ch));
        delay_uplink_T(1,t) = delay_uplink;
        delay_uplink_T(2,t) = sum(delay_uplink_T(1,:),2)/t;

        costq_last = costq_t;

        % 计算当前时隙用户的任务完成时延以及基站的缓存开销（用更新后的缓存状态计算）
        delay_pro = max((req_type' * X_t).*C_t) * delay_edge + (1 - max((req_type' * X_t).*C_t)) * delay_bkb;
        delay_t = delay_pro + delay_uplink;
        delay_T(1,t) = delay_t;
        delay_T(2,t) = sum(delay_T(1,:),2)/t;
        cach_cost_T(1,t) = sum(xi'*X_t,2);
        cach_cost_T(2,t) = sum(cach_cost_T(1,:),2)/t;
        %% comparison: opti
        [C_opti_t, sinr_u] = opti_cluster(ch);
        a_t_opti = xi' * X_t_opti * (C_opti_t');
        costq_t_opti = max(costq_last_opti + a_t_opti - cost_th, 0.1);
        % 当前时隙变量求解
        X_t_opti = uplink_optimal(C_opti_t)
        % 计算上行链路时延
        delay_uplink_opti = data_t / log2(1 + sinr_u);
        delay_uplink_T_opti(1,t) = delay_uplink_opti;
        delay_uplink_T_opti(2,t) = sum(delay_uplink_T_opti(1,:),2)/t;

        costq_last_opti = costq_t_opti;

        % 计算当前时隙用户的任务完成时延以及基站的缓存开销（用更新后的缓存状态计算）
        delay_pro_opti = max((req_type' * X_t_opti).*C_opti_t) * delay_edge + (1 - max((req_type' * X_t_opti).*C_opti_t)) * delay_bkb;
        delay_t_opti = delay_pro_opti + delay_uplink_opti;
        delay_T_opti(1,t) = delay_t_opti;
        delay_T_opti(2,t) = sum(delay_T_opti(1,:),2)/t;
        cach_cost_T_opti(1,t) = sum(xi'*X_t_opti,2);
        cach_cost_T_opti(2,t) = sum(cach_cost_T_opti(1,:),2)/t;
        %% 对比算法：Lyapunov-based 瞬时最优
        a_t_instant = xi' * X_t_instant * (C_t_instant');
        costq_t_instant = max(costq_last_instant + a_t_instant - cost_th, 0.1);
        % 当前时隙变量求解
        [C_t_instant, X_t_instant] = instant_optimal()
        % 计算上行链路时延
        delay_uplink_instant = data_t / log2(1 + calcu_sinr(clus_all, C_t_instant, ch));
        delay_uplink_T_instant(1,t) = delay_uplink_instant;
        delay_uplink_T_instant(2,t) = sum(delay_uplink_T_instant(1,:),2)/t;

        costq_last_instant = costq_t_instant;

        % 计算当前时隙用户的任务完成时延以及基站的缓存开销（用更新后的缓存状态计算）
        delay_pro_instant = max((req_type' * X_t_instant).*C_t_instant) * delay_edge + (1 - max((req_type' * X_t_instant).*C_t_instant)) * delay_bkb;
        delay_t_instant = delay_pro_instant + delay_uplink_instant;
        delay_T_instant(1,t) = delay_t_instant;
        delay_T_instant(2,t) = sum(delay_T_instant(1,:),2)/t;
        cach_cost_T_instant(1,t) = sum(xi'*X_t_instant,2);
        cach_cost_T_instant(2,t) = sum(cach_cost_T_instant(1,:),2)/t;
        %% 对比算法：block descent
        a_t_block = xi' * X_t_block * (C_t_block');
        costq_t_block = max(costq_last_block + a_t_block - cost_th, 0.1);
        % 当前时隙变量求解
        [C_t_block, X_t_block] = block_descent()
        % 计算上行链路时延
        delay_uplink_block = data_t / log2(1 + calcu_sinr(clus_all, C_t_block, ch));
        delay_uplink_T_block(1,t) = delay_uplink_block;
        delay_uplink_T_block(2,t) = sum(delay_uplink_T_block(1,:),2)/t;

        costq_last_block = costq_t_block;

        % 计算当前时隙用户的任务完成时延以及基站的缓存开销（用更新后的缓存状态计算）
        delay_pro_block = max((req_type' * X_t_block).*C_t_block) * delay_edge + (1 - max((req_type' * X_t_block).*C_t_block)) * delay_bkb;
        delay_t_block = delay_pro_block + delay_uplink_block;
        delay_T_block(1,t) = delay_t_block;
        delay_T_block(2,t) = sum(delay_T_block(1,:),2)/t;
        cach_cost_T_block(1,t) = sum(xi'*X_t_block,2);
        cach_cost_T_block(2,t) = sum(cach_cost_T_block(1,:),2)/t;
        
        theta = 1;
    end
    delay_B(1,B) = delay_T(2,T);
    delay_uplink_B(1,B) = delay_uplink_T(2,T);
    cach_cost_B(1,B) = cach_cost_T(2,T);
    delay_B_opti(1,B) = delay_T_opti(2,T);
    delay_uplink_B_opti(1,B) = delay_uplink_T_opti(2,T);
    cach_cost_B_opti(1,B) = cach_cost_T_opti(2,T);
    delay_B_instant(1,B) = delay_T_instant(2,T);
    delay_uplink_B_instant(1,B) = delay_uplink_T_instant(2,T);
    cach_cost_B_instant(1,B) = cach_cost_T_instant(2,T);
    delay_B_block(1,B) = delay_T_block(2,T);
    delay_uplink_B_block(1,B) = delay_uplink_T_block(2,T);
    cach_cost_B_block(1,B) = cach_cost_T_block(2,T);
% delay_B(1,B) = delay_T(2,T);
% cach_cost_B(1,B) = cach_cost(2,T);
% delay_uplink_B(1,B) = delay_uplink_T(2,T);
%     th_ind  = th_ind + 1;
%% 保存性能数据
% path = sprintf('%s', ['output\' 'diffclus_size' version]);
% save(path, 'delay_B','cach_cost_B','delay_uplink_B');
disp(time);
time = datetime;
disp(time);
end
% hold off;
% legend;

figure();
plot(delay_B, 'r-','Marker','o','LineWidth',1,'DisplayName','proposed(total)');
hold on;
plot(delay_B_instant,'b-','Marker', '*','LineWidth',1,'DisplayName','instant(total)');
hold on;
plot(delay_B_block,'g-','Marker', 'v','LineWidth',1,'DisplayName','block descent(total)');
hold on;
plot(delay_B_opti,'m-','Marker', 'd','LineWidth',1,'DisplayName','uplink optimal(total)');
hold on;
plot(delay_uplink_B,'r--','Marker', 'o','LineWidth',1,'DisplayName','proposed(uplink)');
hold on;
plot(delay_uplink_B_instant,'b--','Marker', '*','LineWidth',1,'DisplayName','instant(uplink)');
hold on;
plot(delay_uplink_B_block,'g--','Marker', 'v','LineWidth',1,'DisplayName','block descent(uplink)');
hold on;
plot(delay_uplink_B_opti,'m--','Marker', 'd','LineWidth',1,'DisplayName','uplink optimal(uplink)');
grid on;
xlabel('size of cluster');
ylabel('30 slots averaged delay');
legend;

figure();
plot(cach_cost_B,'r-','Marker', 'o','LineWidth',1,'DisplayName','proposed');
hold on;
plot(cost_th*ones(1,B_max),'k-','Marker', '.','LineWidth',1,'DisplayName','threshold');
hold on;
plot(cach_cost_B_instant,'b--','Marker', '*','LineWidth',1,'DisplayName','instant');
hold on;
plot(cach_cost_B_block,'g--','Marker', 'v','LineWidth',1,'DisplayName','block descent');
hold on;
plot(cach_cost_B_opti,'m--','Marker', 'd','LineWidth',1,'DisplayName','uplink optimal');
title('caching cost');
grid on;
xlabel('size of cluster');
ylabel('30 slots averaged cost');
legend;
