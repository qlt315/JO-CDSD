%% 对目标用户基站簇的划分（同时也会更新其他用户的基站簇划分方试，通过随机方式产生）
% 输入;信道系数
% 输出：所有用户的基站簇划分，以及目标用户的信干噪比
function [clus_all, sinr_u] = BScluster(ch)
    %% 声明全局变量
    global N_agents;
    global B;
    global A;
    global typical_user;
    %% b&b
    clus = zeros(1,N_agents); % 初始化解向量

    % 生成其他用户的分簇
    % 计算每个基站的负荷量，如果所有基站全部满载，需要重新生成
    overload = true;
    while overload
        clus_all = clus_state();
        overload = false;
        clus_all(typical_user,:) = clus;
        load = sum(clus_all);
        cnt = 0;
        for bs=1:N_agents
            if load(bs)>A
                overload = true;
                continue;
            end
            if load(bs)==A
                cnt = cnt+1;
            end
        end
        if cnt==N_agents
            overload = true;
        end
    end
    
    % 遍历搜索初始化
    max_sinr = 0;
    clus = zeros(1,N_agents);
    for ind = 1:1023
        node_solu = dec2bin(ind)-'0';
        if size(node_solu,2) < N_agents
            less_len = N_agents - size(node_solu,2);
            prex = zeros(1,less_len);
            node_solu = [prex, node_solu];
        end
        clus_all(typical_user,:) = node_solu;
        if sum(node_solu) > B || max(sum(clus_all,1)) > A
            continue;
        end
        cur_sinr = calcu_sinr(clus_all, node_solu, ch);
        if max_sinr < cur_sinr
            max_sinr = cur_sinr;
            clus = node_solu;
        end
    end
    sinr_u = max_sinr;
    clus_all(typical_user,:) = clus;
end
