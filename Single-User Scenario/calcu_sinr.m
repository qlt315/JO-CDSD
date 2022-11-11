%% 根据分簇情况以及当前信道计算tipical用户的信干噪比
function sinr_u=calcu_sinr(clus_all, clus_u, ch)
global N_agents;
global A;
global P;
global U;
global typical_user;

clus_all(typical_user,:) = clus_u;
if clus_all(typical_user,:)==zeros(1,N_agents)
    sinr_u = 0;
    return;
end
clus_size = sum(clus_all(typical_user,:));
I = eye(A*clus_size, A*clus_size);
% 求出G(phi_u,omega_u)矩阵，维度为|phi_u|*A,|omega_u|-1，不包括用户u
G_u = ch;
for bs=N_agents:-1:1
    if clus_all(typical_user,bs)==0
        G_u((bs-1)*A+1:bs*A,:) = [];
    end
end
for user=U:-1:1
    if sum(clus_all(user,:).*clus_all(typical_user,:))==0 || user==typical_user
        G_u(:,user) = [];
    end
end
% 求g_u^u=g_u(typical_user)，维度为|phi_u|*A,U
g_u = ch;
for bs=N_agents:-1:1
    if clus_all(typical_user,bs)==0
        g_u((bs-1)*A+1:bs*A,:) = [];
    end
end
% 求用户u的波束成形矢量
temp = (I-G_u*pinv(G_u))*g_u(:,typical_user);
beamform_u = temp./norm(temp);
% test
test = beamform_u' * G_u; %结果应该十分接近0，表示intra干扰消除   
% 求用户u的信干噪比
signal = P*abs(beamform_u'*g_u(:,typical_user))^2;
interf = 0;
for user=1:U
    if sum(clus_all(user,:).*clus_all(typical_user,:))==0
        interf = interf + P*abs(beamform_u'*g_u(:,user))^2;
    end
end
noise = 10^(-3.5);
sinr_u = signal/(interf+noise^2);
end