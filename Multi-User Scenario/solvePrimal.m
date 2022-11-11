function [primal_solu, lambda_fea, feasible, opti_val] = solvePrimal(clus_u)
global K;
global N_agents;

feasible = 1;
VLB = zeros(K, N_agents);
VUB = ones(K, N_agents);
X_init = ones(K, N_agents);
% problem = createOptimProblem('fmincon','objective',@(X) funcPrimal(X,clus_u),'x0',X_init,'lb',VLB,'ub',VUB,'nonlcon',@(X) consPrimal(X, clus_u));
% gs = GlobalSearch;
% [primal_solu, opti_val, exitflag] = run(gs, problem);
% if exitflag == -2
%     feasible = 0;
%     lambda_fea = [];
% else  % get the multipiers
%     problem_lambda = createOptimProblem('fmincon','objective',@(miu) funcLambda(miu,primal_solu,clus_u),'x0',zeros(N_agents*2+1, 1),'lb',zeros(N_agents*2+1, 1));
%     gs_lambda = GlobalSearch;
%     [lambda_fea, fval_lambda] = run(gs_lambda, problem_lambda);
% end

options = optimoptions('fmincon','MaxIter',100000,'MaxFunEvals',100000);
[primal_solu, opti_val, exitflag, ~, lambda] = fmincon(@(X) funcPrimal(X,clus_u),...
    X_init,[],[],[],[],VLB,VUB,@(X) consPrimal(X, clus_u), options);
if exitflag == -2
    feasible = 0;
    lambda_fea = [];
else 
    lambda_fea = [lambda.ineqnonlin; lambda.eqnonlin];
end
end