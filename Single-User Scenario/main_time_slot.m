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
% B_max = 5;
A = 3; % 每个基站天线个数
P = 1; % user's transmit power
N_agents = 10; %number of BSs, i.e., M
% C_bs = randi([2,5],N_agents,1);
% S_bs = randi([2,5],N_agents,1);%computing and storage of BSs
C_bs = 3*ones(N_agents,1); %单位：GB，变化范围[3-10]
S_bs = 3*ones(N_agents,1);%computing and storage of BSs，单位：GHz
U = 4; %number of users
typical_user = randi(U); % 指定typical用户
% sinr_th = 1;
distance = kron(0.1 * randi([1,10],N_agents,U), ones(A,1)); %基站到用户的距离，单位：Km
beta = 128.1 + 37.6 * log(distance);
pathLoss = 10 ^ (beta / 10);
%% parameters of task
K = 6; %number of service types
% size_se = randi([2,5],K,1); %size vector of services
% comp_se = randi([2,5],K,1); %computing requirement of services
% cost_coef = rand(K,1); %cost coefficient of caching services
size_se = 3 * ones(K,1); %单位：GB
comp_se = 0.3 * ones(K,1); %单位：GHz
data_se = linspace(1,K,K); %单位：GB
workld_se = linspace(0.1,0.1*K,K); %处理每单位数据需要的CPU cycles，单位：CPU cycles/GB？
cost_coef = linspace(0.1,0.1*K,K); %下载每单位数据需要花费的开销，单位：千元/GB
R = 0.05; %data rate of backbone，单位：Gbps
xi = cost_coef' .* size_se;
cost_th = 2;

T = 30;
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
% comparison 2: without Lyapunov dynamic
C_t_perslot = zeros(1,N_agents);
X_t_perslot = zeros(K, N_agents);
delay_T_perslot = zeros(2,T);
delay_uplink_T_perslot = zeros(2,T);
cach_cost_T_perslot = zeros(2,T);

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
    ch_T(:,:,t) = ch .* sqrt(pathLoss);
    S_bs_T(:,t) = round(rand(N_agents,1)) * 3;
    C_bs_T(:,t) = round(rand(N_agents,1)) * 3;
end
for t = 1:T
    [clus_all_T(:,:,t),~] = BScluster(sum(ch_T(:,:,t),3));
end
% begin long term optimization 
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
    C_opti_t = clus_all(typical_user,:);
    %% comparison: opti
    a_t_opti = xi' * X_t_opti * (C_opti_t');
    costq_t_opti = max(costq_last_opti + a_t_opti - cost_th, 0.1);
    % 当前时隙变量求解
    X_t_opti = uplink_optimal(C_opti_t)
    % 计算上行链路时延
    delay_uplink_opti = data_t / log2(1 + calcu_sinr(clus_all, C_opti_t, ch));
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
    
    %% proposed algorithm
    a_t = xi' * X_t * (C_t');
    costq_t = max(costq_last + a_t - cost_th, 0.1);
    % 当前时隙变量求解
    [C_t, X_t, feasible_flag] = benders()
    while feasible_flag == 0
        clus_all = BScluster(ch);
        clus_all_T(:,:,t) = clus_all;
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
    
    %% 对比算法2：每时隙最优
    % 当前时隙变量求解
    [C_t_perslot, X_t_perslot] = perslot_optimal()
    % 计算上行链路时延
    delay_uplink_perslot = data_t / log2(1 + calcu_sinr(clus_all, C_t_perslot, ch));
    delay_uplink_T_perslot(1,t) = delay_uplink_perslot;
    delay_uplink_T_perslot(2,t) = sum(delay_uplink_T_perslot(1,:),2)/t;

    % 计算当前时隙用户的任务完成时延以及基站的缓存开销（用更新后的缓存状态计算）
    delay_pro_perslot = max((req_type' * X_t_perslot).*C_t_perslot) * delay_edge + (1 - max((req_type' * X_t_perslot).*C_t_perslot)) * delay_bkb;
    delay_t_perslot = delay_pro_perslot + delay_uplink_perslot;
    delay_T_perslot(1,t) = delay_t_perslot;
    delay_T_perslot(2,t) = sum(delay_T_perslot(1,:),2)/t;
    cach_cost_T_perslot(1,t) = sum(xi'*X_t_perslot,2);
    cach_cost_T_perslot(2,t) = sum(cach_cost_T_perslot(1,:),2)/t;
    
    theta = 1;
end
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
%% 画图
figure;
plot(delay_T(2,:),'r-o','Linewidth',1);
hold on;
plot(delay_uplink_T(2,:),'r--o','Linewidth',1);
hold on;
plot(delay_T_instant(2,:),'b-*','Linewidth',1);
hold on;
plot(delay_uplink_T_instant(2,:),'b--*','Linewidth',1);
hold on;
plot(delay_T_perslot(2,:),'g-d','Linewidth',1);
hold on;
plot(delay_uplink_T_perslot(2,:),'g--d','Linewidth',1);
hold on;
plot(delay_T_opti(2,:),'m-+','Linewidth',1);
hold on;
plot(delay_uplink_T_opti(2,:),'m--+','Linewidth',1);
hold on;
plot(delay_T_block(2,:),'y-x','Linewidth',1);
hold on;
plot(delay_uplink_T_block(2,:),'y--x','Linewidth',1);
grid on;
xlabel('time slots');
ylabel('averaged delay');
legend('proposed(total)', 'porposed(uplink)','instant(total)', 'instant(uplink)', 'perslot(total)', 'perslot(uplink)', 'uplink optimal(total)', 'uplink optimal(uplink)', 'block descent(total)', 'block descent(uplink)');

figure;
plot(cach_cost_T(2,:),'r-o','Linewidth',1);
hold on;
plot(cost_th*ones(1,T),'k-','Linewidth',1);
hold on;
plot(cach_cost_T_instant(2,:),'b--d','Linewidth',1);
hold on;
plot(cach_cost_T_perslot(2,:),'g--+','Linewidth',1);
hold on;
plot(cach_cost_T_opti(2,:),'m--v','Linewidth',1);
hold on;
plot(cach_cost_T_block(2,:),'y--x','Linewidth',1);
grid on;
xlabel('time slots');
ylabel('averaged cost');
legend('proposed', 'threshold', 'instant', 'perslot','uplink optimal','block descent');

