function clus_all=clus_state()
%% generate a random clustering state of other users that meets the constraints
%% output: the clustering state of the other users
    global N_agents;
    global B;
    global A;
    global U;
    global typical_user;

    % 其他用户当前时隙的分簇情况,用户u的分簇初始化为0
    clus_all = zeros(U,N_agents);
%     clus_all(typical_user,:) = clus_u;
    for i=1:U
        if i==typical_user
            continue;
        end
        for j=1:N_agents
            attenta = sum(clus_all);
            clus_size = sum(clus_all,2);
            if (attenta(j)<A) && (clus_size(i)<B)
                seed = rand();
                if seed>0.5
                    clus_all(i,j) = 1;
                else
                    clus_all(i,j) = 0;
                end
            else 
                clus_all(i,j) = 0;
            end
        end
    end
end
    