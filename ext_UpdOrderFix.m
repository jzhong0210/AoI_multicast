function [delay, delay_max, success] = ext_UpdOrderFix(write,read,delayset)
    %% soft deadline, commit after write-th receiver finishes
    % return a specific update system time regardless it's in the write quorum or not
    % delay: service time of an update
    % delay_max: max service time of the write quorum

    n = length(delayset);
    
    % obtain the index of all updates in the read quorum
    % readind = randperm(n,read,1);
    readind = 1; % save time for read==1
    % minimum of the read quorum
    delay = min(delayset(readind));
    delay_max = max(delayset(randperm(n,write))); % pick write numbers without replacement
    
    if (delay<=delay_max)
        success = 1;
    else
        success = 0;
        delay = delay_max;
    end

end