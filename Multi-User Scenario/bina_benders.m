function [C_t, X_t] = bina_benders()
global delay_edge;
global delay_bkb;
global B0_size;
global K;
global theta;
global clus_all;
global ch;
global data_t;
global req_type;

theta_l = 0;
theta_r = 0;
for kk=1:K
    theta_l = theta_l + delay_edge(kk);
    theta_r = theta_r + delay_bkb(kk);
end
theta_l = theta_l / B0_size;
theta_r = theta_r / B0_size;
theta = theta_r;
[C_t, X_t, ~] = benders(); % 先用最差的计算一个初始值
delay_uplink = 0;
for kk=1:K
    delay_uplink = delay_uplink + data_t(kk) / log2(1 + calcu_sinr(clus_all, C_t, ch));
end
delay_uplink = delay_uplink / B0_size;
delay_pro = 0;
for kk=1:K
    delay_pro = delay_pro + max(C_t .* (req_type(kk) * X_t(kk,:))) * (delay_edge(kk) - delay_bkb(kk)) + delay_bkb(kk);
end
delay_pro = delay_pro / B0_size;
min_delay = delay_pro + delay_uplink;
% intialization
theta_mid = (theta_l + theta_r) / 2;
epsilon = 0.1;
error = inf;
iter_max = 12;
iter = 1;
while error > epsilon && iter <= iter_max
    theta = theta_mid;
    [C_t_cur, X_t_cur, feasible_flag] = benders();
    if feasible_flag == 0
        theta_l = theta;
        theta_mid = (theta_l + theta_r) / 2;
    else
        % 计算当前解的总时延
        delay_uplink = 0;
        for kk=1:K
            delay_uplink = delay_uplink + data_t(kk) / log2(1 + calcu_sinr(clus_all, C_t_cur, ch));
        end
        delay_uplink = delay_uplink / B0_size;
        delay_pro = 0;
        for kk=1:K
            delay_pro = delay_pro + max(C_t_cur .* (req_type(kk) * X_t_cur(kk,:))) * (delay_edge(kk) - delay_bkb(kk)) + delay_bkb(kk);
        end
        delay_pro = delay_pro / B0_size;
        delay_t = delay_pro + delay_uplink;
        if delay_t < min_delay
            C_t = C_t_cur;
            X_t = X_t_cur;
        end
        theta_r = theta;
        theta_mid = (theta_l + theta_r) / 2;
        error = abs(theta - theta_mid);
    end
    iter = iter + 1;
end
end