clear all;

lambda = 1;
alpha = 3;
n = 1000;

eosappr = zeros(n,1);
eos = zeros(n,1);


for k=1:n
    eos(k) = ParetoOS(lambda,alpha,1,k,n);
    a = 1/alpha;
    betapr = (n-k)/n;
    % eosappr(k) = lambda*exp(1/alpha-1-(1/alpha-1)/(betapr))*sqrt((betapr*n+1)/(betapr*n+1-1/alpha));
    % eosappr(k) =
    % lambda*sterappr(n)/sterappr(n-k)*sterappr(n-k-1/alpha)/sterappr(n-1/alpha);
    % sterling approximation
    
    eosappr(k) =lambda*sqrt((betapr*n-a)/(betapr*n))*betapr^(-a); 
end

figure();
hold on;
plot(1:n,eos);
plot(1:n,eosappr);
axis([0 n 0 50]);
legend('eos','eosappr');


function [nster] = sterappr(n)
    nster = sqrt(2*pi*n)*(n/exp(1))^n;
end