function [delay, delay_max, success] = ext_UpdOrderFixFree(write,read,delayset)
    %% soft deadline, commit after write-th receiver finishes
    % return a specific update system time regardless it's in the write quorum or not
    % premium group vs free group
    % delay: service time of an update
    % delay_max: max service time of the write quorum
    n = length(delayset);
    
    % obtain the index of all updates in the read quorum
    if write<n
        readind = write+1;
    else 
        readind = write;
    end
    % minimum of the read quorum
    delay = delayset(readind);
    
    delay_max = max(delayset(1:write));
    
    if (delay<=delay_max)
        success = 1;
    else
        success = 0;
        delay = delay_max;
    end

end