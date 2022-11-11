function [primal_solu, lambda_infea] = solveInfeaPrimal(clus_u)
global K;
global N_agents;

VLB(1:K,:) = zeros(K, N_agents);
VUB(1:K,:) = ones(K, N_agents);
X_init(1:K,:) = ones(K, N_agents);
% problem = createOptimProblem('fmincon','objective',@(X) funcInfeaPrimal(X,clus_u),'x0',X_init,'lb',VLB,'ub',VUB,'nonlcon',@(X) consInfeaPrimal(X, clus_u));
% gs = GlobalSearch;
% [primal_solu, ~, ~, ~, allmins] = run(gs, problem);
options = optimoptions('fmincon','MaxIter',100000,'MaxFunEvals',100000);
[primal_solu, ~, ~, ~, lambda] = fmincon(@(X) funcInfeaPrimal(X,clus_u),...
    X_init,[],[],[],[],VLB,VUB,@(X) consInfeaPrimal(X, clus_u), options);
lambda_infea = [lambda.ineqnonlin(2:N_agents*2+1); lambda.eqnonlin];
end