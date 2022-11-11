function req = gener_req(T, K, rho)
zipfpdf = @(x) x.^(-(rho + 1)) ./ zeta(rho + 1);
req = zeros(1,T);
types = linspace(1,K,K);
probas = zipfpdf(types);
sum_proba = trapz(types, probas);
% probas = probas / sum_proba;
t = 1;
while t<=T
    ran_num = randi(K);
    pro_ran = zipfpdf(ran_num) / sum_proba;
    proba = rand(1);
    if proba <= pro_ran
        req(t) = ran_num;
        t = t+1;
    end
end
end