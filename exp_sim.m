%% Simulation for m multi user case, exponential distribution
% vary the number of users

clear all;
close all;

% erasrue rate
lambda = 5;

% number of symbols in an update 
k = 1;
% k = 1;
% number of users
marray = 1:50;
% m = 500;


% length of update messages
msglen = 10000;


%% IIR (simulation of analysis)
fprintf('IR pool \n');
% pool_iir = k + nbinrnd(k,1-delta,length(marray)*msglen,1);
pool_ir = exprnd(1/lambda,length(marray)*msglen,1);

normage_ir = zeros(length(marray),1);
fprintf('IR \n');
for m_ind = 1:length(marray)
    m = marray(m_ind);
    n1 = zeros(msglen,1);
    nmax = zeros(msglen,1);
    for msgind = 1:msglen
        nmind = randi(length(pool_ir),m,1);
        nm = pool_ir(nmind);
        n1(msgind) = nm(1);
        nmax(msgind) = max(nm);
    end
    age_ir = mean(n1) + mean(nmax.^2)/(2*mean(nmax));
    normage_ir(m_ind) = age_ir/k;
end


%% FIR with feedback (true experiments in a sample pool)

fprintf('FR pool \n');

poolsize = msglen*m/100;
% pool = gen_update(delta,poolsize,poollen);
% pool = k + nbinrnd(k,1-delta,poolsize,1);
pool_fr = exprnd(1/lambda,length(marray)*msglen,1);


normage_fr = zeros(length(marray),1);
optth_fr = zeros(length(marray),1);

fprintf('FR \n');
for m_ind = 1:length(marray)
    m = marray(m_ind);
    
    % start the search for optimal n_th
    age_pending = inf;
    
    
    % n out of m (FR) threshold
    narray_th = 1:m;
    
    for n_th = narray_th
        endpt = 0;
        polygons = 0;
        response = 0;
        for msgind = 1:msglen
            [n_fir, nmax_fir, suc_fir] = ext_update(k,n_th,m,pool_fr);
            if suc_fir==0 % fail
                response = response + n_fir;
            else
                polygons = polygons + 1/2*((response+n_fir)^2-n_fir^2);
                % endpt = endpt + nmax_fir;
                endpt = endpt + response;
                % response = nmax_fir-n_fir;
                response = nmax_fir;
            end        
        end
        age_fir = polygons/endpt;
        if age_fir<age_pending
            optth_fr(m_ind) = n_th;
        end
        age_pending = min(age_pending,age_fir);
    end
    normage_fr(m_ind) = age_pending/k;
end


%% FR


%% plot

figure(1)
set(gcf,'units','pixels','position',[10,10,400,200]);
hold on;
plot(marray,normage_ir,'.-','linewidth',1.5);
% plot(marray,normage_fr,'.-','linewidth',1.5);
plot(marray,normage_fr,'.-','linewidth',1.5);
xlabel('number of users: \it{m}','Fontsize',14,'FontName','Times');
ylabel('normalized age','Fontsize',14,'FontName','Times');
title('\delta = 0.2','Fontsize',14,'FontName','Times');
leg = legend('IR','FR','location','Northwest');
set(leg,'Fontsize',12,'FontName','Times');
% axis([0 2000 1.5 1.54]);
grid on;