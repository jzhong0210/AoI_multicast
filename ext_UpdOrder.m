function [delay, delay_max, success] = ext_UpdOrder(write,read,delayset)
    %% soft deadline, commit after write-th receiver finishes, return the min among read number of delay
    % delay: service time of an update
    % delay_max: max service time of the write quorum

    n = length(delayset);
    
    % obtain the index of all updates in the read quorum
    % readind = randperm(n,read,1);
    readind = 1; % save time for read==1
    % minimum of the read quorum
    delay = min(delayset(readind));
    % sort all n delays 
    [delayset_sort, ~] = sort(delayset);
    % the order of a particular node i
    order = find(delayset_sort==delay,1);
    
    if (order<=write)
        success = 1;
    else
        success = 0;
        delay = delayset_sort(write);
    end
    
    % max delay over all write quorum
    delay_max = delayset_sort(write);
end