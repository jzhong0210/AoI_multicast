%% Simulation for m multi user case, geometric distribution

clear all;
close all;

% erasrue rate
delta = 0.2;

% number of symbols in an update 
k = 1;
% k = 1;
% number of users
marray = 1:20;
% m = 500;

% FIR & FR thresholde
narray_th = 1:20;

% length of update messages
msglen = 10000;


%% IIR (simulation of analysis)
fprintf('IR pool \n');
pool_iir = k + nbinrnd(k,1-delta,length(marray)*msglen,1);

normage_iir = zeros(length(marray),1);
fprintf('IR \n');
for m_ind = 1:length(marray)
    m = marray(m_ind);
    n1 = zeros(msglen,1);
    nmax = zeros(msglen,1);
    for msgind = 1:msglen
        nmind = randi(length(pool_iir),m,1);
        nm = pool_iir(nmind);
        n1(msgind) = nm(1);
        nmax(msgind) = max(nm);
    end
    age_iir = mean(n1) + mean(nmax.^2)/(2*mean(nmax));
    normage_iir(m_ind) = age_iir/k*(1-delta);
end


%% FIR with feedback (true experiments in a sample pool)

fprintf('FR pool \n');

poolsize = msglen*m/100;
% pool = gen_update(delta,poolsize,poollen);
pool = k + nbinrnd(k,1-delta,poolsize,1);

normage_fir = zeros(length(marray),1);
optth_fir = zeros(length(marray),1);

fprintf('FR \n');
for m_ind = 1:length(marray)
    m = marray(m_ind);
    
    % start the search for optimal n_th
    age_pending = inf;
    for n_th = narray_th
        endpt = 0;
        polygons = 0;
        response = 0;
        for msgind = 1:msglen
            [n_fir, nmax_fir, suc_fir] = ext_update(k,n_th,m,pool);
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
            optth_fir(m_ind) = n_th;
        end
        age_pending = min(age_pending,age_fir);
    end
    normage_fir(m_ind) = age_pending/k*(1-delta);
end


%% FR


%% plot

figure(1)
set(gcf,'units','pixels','position',[10,10,400,200]);
hold on;
plot(marray,normage_iir,'.-','linewidth',1.5);
% plot(marray,normage_fr,'.-','linewidth',1.5);
plot(marray,normage_fir,'.-','linewidth',1.5);
xlabel('number of users: \it{m}','Fontsize',14,'FontName','Times');
ylabel('normalized age','Fontsize',14,'FontName','Times');
title('\delta = 0.2','Fontsize',14,'FontName','Times');
leg = legend('IR','FR','location','Northwest');
set(leg,'Fontsize',12,'FontName','Times');
% axis([0 2000 1.5 1.54]);
grid on;