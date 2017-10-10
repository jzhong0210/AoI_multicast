%% pareto order statistics

function [os] = ParetoOS(lambda,alpha,m,k,n)
    intlist = n-k+1:1:n;
    offsetlist = n-k+1-m/alpha:1:n-m/alpha;
    fraclist = intlist./offsetlist;
    os = lambda^m*prod(fraclist);
    % os = lambda^m*gamma(n+1)/gamma(n-k+1)*gamma(n-k+1-m/alpha)/gamma(n+1-m/alpha);
end
