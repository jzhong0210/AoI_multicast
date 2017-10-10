function [delay, delay_max, success] = ext_UpdOrderFix(write,read,n,pool)
    %% soft deadline, commit after write-th receiver finishes
    % return a specific update system time regardless it's in the write quorum or not
    % delay: service time of an update
    % delay_max: max service time of the write quorum

    % obtain n different delay values from the pool
    delayind = randi(length(pool),n,1);
    delayset = pool(delayind);
    
    % obtain the index of all updates in the read quorum
    readind = randi(n,read,1);
    % minimum of the read quorum
    delay = min(delayset(readind));
    delay_max = max(delayset(randi(n,write,1)));
    
    if (delay<=delay_max)
        success = 1;
    else
        success = 0;
        delay = delay_max;
    end

end