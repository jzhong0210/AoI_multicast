function [n, nmax, success] = ext_update(k,threshold,m,pool)
    %% hard deadline, restart after time=threshold
    % n: service time of an update with size k
    % nmax: service time of m updates with size k, the max over all m
        
    % FIR: wait until threshold 
    
    nmind = randi(length(pool),m,1);
    nm = pool(nmind);
    n = nm(1);
    if (n<=threshold)
        success = 1;
    else
        success = 0;
        n = threshold;
    end
    nmax = max(nm);
    if (nmax>=threshold)
        nmax = threshold;
    end 
end