%% Simulation for n multi user case, Pareto distribution
% fix the number of users, vary the threshold

clear all;
close all;

% exponential rate
pshape_list = [1.2 1.4 2];

% number of users
n = 100;
% m = 500;
write_list = 1:n;
% constant shift 
C = 0.5;

% number of update messages
msglen = 1000;

% read quorum
read = 1;


%% examine different lambda

avgage = zeros(length(write_list),3);
apprage = zeros(length(write_list),3);
beta_opt = zeros(1,3);


for pshape_ind = 1:length(pshape_list)

    pshape = pshape_list(pshape_ind);
    
    %% preemptive multicast: true experiments obtained by a sample pool

    % generate the pool for transmission time
    poolsize = msglen*n/100;% pool for pareto distribution with pscale \gamma and pshape \alpha
    pscale = 1;
    % parameters for generating pareto
    Genpar_k = 1/pshape;
    Genpar_sigma = pscale*Genpar_k;
    Genpar_theta = pscale;
    pool = gprnd(Genpar_k,Genpar_sigma,Genpar_theta,msglen,n);


    for write = write_list
        endpt = 0;
        polygons = 0;
        response = 0;
        for msgind = 1:msglen
            delayset = pool(msgind,:);
            [delay, delay_max, suc] = ext_UpdOrder(write,read,delayset);
            if suc==0 % fail
                response = response + delay;
            else
                polygons = polygons + 1/2*((response+delay)^2-delay^2);
                % endpt = endpt + nmax_fir;
                endpt = endpt + response;
                % response = nmax_fir-n_fir;
                response = delay_max;
            end        
        end
        avgage(write,pshape_ind) = polygons/endpt;
    end

    
    %% wait for fixed k out of n: preemptive multicast experiment
    
    %{
    for write = write_list
        endpt = 0;
        polygons = 0;
        response = 0;
        for msgind = 1:msglen
            [delay, delay_max, suc] = ext_UpdOrderFix(write,read,n,pool);
            if suc==0 % fail
                response = response + delay;
            else
                polygons = polygons + 1/2*((response+delay)^2-delay^2);
                % endpt = endpt + nmax_fir;
                endpt = endpt + response;
                % response = nmax_fir-n_fir;
                response = delay_max;
            end        
        end
        avgage_fix(write,lambda_ind) = polygons/endpt;
    end
    %}
    
    %% age approximation for pareto

    for k = write_list
        beta = k/n;
        a = 1/pshape; 
        osappr = pscale*(1-beta)^(-a); %*sqrt(((1-beta)*m-a)/((1-beta)*m))
        apprage(k,pshape_ind) = pshape*pscale*n/k/(pshape-1) + (pshape*k-2*n+k)/(2*k*(pshape-1))*osappr;
    end
    
    %% solve the optimal stopping threshold
    syms x
    eqn = (1+1/pshape)*x^2-x/pshape-(1-x)^(1/pshape) == 0;
    solx = vpasolve(eqn,x,[0 1]);
    beta_opt(pshape_ind) = vpa(solx);

end

%% plot the Pareto distribution

blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
purple = [0.4940 0.1840 0.5560];
color = {blue, red, purple};

figure(1)
set(gcf,'units','pixels','position',[10,10,400,250]);
hold on;
for i = 1:length(pshape_list)
    l(i) = plot(write_list(2:2:n),avgage(2:2:n,i),'-','Color',color{i},'linewidth',1.5);
    plot(write_list(2:2:n),apprage(2:2:n,i),'x','Color',color{i},'linewidth',1.5);
    plot(beta_opt(i)*n,apprage(floor(beta_opt(i)*write),i),'Color',color{i},'Marker','o','MarkerSize',10,'linewidth',2);
end

xlabel('k','Fontsize',14,'FontName','Times');
% ylabel('average age \Delta_{(w,r=50)}','Fontsize',14,'FontName','Times');
ylabel('average age \Delta_{(k)}','Fontsize',14,'FontName','Times');
% title(['preemption at k out of n=' num2str(m)],'Fontsize',14,'FontName','Times');
leg = legend(l(1:3),{['\sigma=' num2str(pshape_list(1))],['\sigma=' num2str(pshape_list(2))],['\sigma=' num2str(pshape_list(3))]},'location','Northeast');
set(leg,'Fontsize',14,'FontName','Times');
axis([0 100 0 8]); 
grid on; box on;