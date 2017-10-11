%% Simulation for n multiple user case, comparing multiple distributions
% fix the number of users, vary the stopping threshold k

clear all;
close all;

% exponential rate
lambda1 = 1;
lambda2 = 6;
p = 0.4;
lambda = 1/(p/lambda1+(1-p)/lambda2);

% number of users
n = 100;
karray = 1:n;
% constant exponential shift (for shifted exponential)
C = 0.0001;

% length of update messages
msglen = 500;

%% generate pool for experiments 
% k out of n with feedback (true experiments in a sample pool)

pool_hyper = [exprnd(1/lambda1,p*msglen,n) ; exprnd(1/lambda2,(1-p)*msglen,n)];
pool_hyper = pool_hyper(randperm(msglen),:);

pool_exp = exprnd(1/lambda,msglen,n);

% pool for pareto distribution with pscale \gamma and pshape \alpha
pscale = 1;
pshape = 1.5;
% parameters for generating pareto
Genpar_k = 1/pshape;
Genpar_sigma = pscale*Genpar_k;
Genpar_theta = pscale;
pool_pareto = gprnd(Genpar_k,Genpar_sigma,Genpar_theta,msglen,n);

% pool for log cauchy distribution with 
% CauchyMu = 1;
% CauchySigma = 1;
% pool_logcau = exp(CauchyMu+CauchySigma*tan(pi*(rand(length(narray)*msglen,1)-1/2)));
% normage = zeros(length(narray),3);

%% hyper exponential

age_hyper = zeros(n,1);   
for k = karray
    endpt = 0;
    polygons = 0;
    response = 0;
    for msgind = 1:msglen
        delayset = pool_hyper(msgind,:);
        [n_fir, nmax_fir, suc_fir] = ext_UpdOrder(k,1,delayset);
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
    age_hyper(k) = polygons/endpt;
end


%% exponential

age_exp = zeros(n,1);  
for k = karray
    endpt = 0;
    polygons = 0;
    response = 0;
    for msgind = 1:msglen
        delayset = pool_exp(msgind,:);
        [n_fir, nmax_fir, suc_fir] = ext_UpdOrder(k,1,delayset);
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
    age_exp(k) = polygons/endpt;
end

%% pareto

age_pareto = zeros(n,1);
for k = karray
    endpt = 0;
    polygons = 0;
    response = 0;
    for msgind = 1:msglen
        delayset = pool_pareto(msgind,:);
        [n_fir, nmax_fir, suc_fir] = ext_UpdOrder(k,1,delayset);
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
    age_pareto(k) = polygons/endpt;
end

%% age approximation for pareto

apprage_pareto = zeros(n,1);
for k = karray
    beta = k/n;
    a = 1/pshape; 
    osappr = pscale*(1-beta)^(-a); %*sqrt(((1-beta)*m-a)/((1-beta)*m))
    apprage_pareto(k) = pshape*pscale*n/k/(pshape-1) + (pshape*k-2*n+k)/(2*k*(pshape-1))*osappr;
    
end

%% exponential approximation by log function

apprage_exp = zeros(length(karray),1);
for k = karray
    apprage_exp(k) = 1/lambda - log(1-k/n)/2/lambda + C/(k/n) - C/2;
end

alpha_opt = sqrt(4*lambda^2*C^2 + 8*lambda*C)/2 - lambda*C;

%% plot the exponential class distribution

blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
purple = [0.4940    0.1840    0.5560];
color = {blue, red, purple};

figure(1)
set(gcf,'units','pixels','position',[10,10,400,250]);
hold on;
plot(karray(2:2:n),age_exp(2:2:n),'-','Color',color{2},'linewidth',1.5);
plot(karray(4:4:n),apprage_exp(4:4:n),'x','Color',color{2},'linewidth',1.5);
plot(karray(2:2:n),age_hyper(2:2:n),'-','Color',color{1},'linewidth',1.5);

%     plot(alpha_opt(i)*n,lambda_list(i)*apprage(floor(alpha_opt(i)*n),i),'Color',color{i},'Marker','o','MarkerSize',10,'linewidth',2);

xlabel('k','Fontsize',14,'FontName','Times');
% ylabel('average age \Delta_{(k,r=1)}','Fontsize',14,'FontName','Times');
ylabel('average age \Delta_{(k)}','Fontsize',14,'FontName','Times');
% title(['preemption at k out of n=' num2str(m)],'Fontsize',14,'FontName','Times');
leg = legend('exponential','log approximation of exp.','hyper-exponential','location','North');
set(leg,'Fontsize',14,'FontName','Times');
axis([0 n 0 2]); 
grid on; box on;

%% plot the Pareto distribution

figure(2)
set(gcf,'units','pixels','position',[10,10,400,250]);
hold on;
plot(karray(1:1:n),age_pareto(1:1:n),'-','Color',color{3},'linewidth',1.5);
plot(karray(4:4:n),apprage_pareto(4:4:n),'x','Color',color{3},'linewidth',1.5);
xlabel('k','Fontsize',14,'FontName','Times');
ylabel('average age \Delta_{(k)}','Fontsize',14,'FontName','Times');
% ylabel('average age \Delta_{(k,r=1)}','Fontsize',14,'FontName','Times');
% title(['preemption at k out of n=' num2str(m)],'Fontsize',14,'FontName','Times');
leg = legend('Pareto \gamma=1, \sigma=5','approximation','location','North');
set(leg,'Fontsize',14,'FontName','Times');
axis([0 n 0 20]); 
grid on; box on;


%% find the optimal point for pareto

% syms beta;
% eqn = pshape/(pshape-1)/beta^2 == ...
%     ( 1/(pshape-1)/beta^2 + ((pshape+1)*beta-2)/(2*pshape*beta*(1-beta)^2) ) * exp(-(1-1/pshape-(1-1/pshape)/(1-beta)));
% solbeta = solve(eqn,beta);
% vpa(solbeta)

% for i=1:100
%     beta = i/100;
%     testbeta(i) = pshape/(pshape-1)/beta^2 - ...
%     ( 1/(pshape-1)/beta^2 + ((pshape+1)*beta-2)/(2*pshape*beta*(1-beta)^2) ) * exp(-(1-1/pshape-(1-1/pshape)/(1-beta)));
% end