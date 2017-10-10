n = 100; write = 100; read = 1;
pool = exprnd(1/2,10000,1)+1;

delayind = randi(length(pool),n,1);
delayset = pool(delayind);

% obtain the index of all updates in the read quorum
readind = randi(n,read,1);
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


delay1 = delay;
delay_max1 = max(delayset(randi(n,write,1)));

if (delay1<=delay_max)
    success = 1;
else
    success = 0;
    delay1 = delay_max;
end