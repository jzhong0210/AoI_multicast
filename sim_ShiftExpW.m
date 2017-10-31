%% Simulation for n multi user case, shifted exponential distribution
% fix the number of users, vary the threshold

clear all;
close all;
rng(31415926);

% exponential rate
lambda_list = [0.5 1 2];

% number of users
n = 100;
% m = 500;
write_list = 1:n;
% write_list = n;
% constant shift 
C = 1;

% number of update messages
msglen = 20000;

% read quorum
read = 1;

% initialization
avgage = zeros(length(write_list),3);
avgage_fix = zeros(length(write_list),3);
avgage_fixfree = zeros(length(write_list),3);
avgage_fixprem = zeros(length(write_list),3);
apprage = zeros(length(write_list),3);
apprage_fix = zeros(length(write_list),3);
alpha_opt = zeros(1,3);
avgage_fr = zeros(length(write_list),3);
apprage_fr = zeros(length(write_list),3);


%% examine different lambda, true experiments obtained by a sample pool

for lambda_ind = 1:length(lambda_list)

    lambda = lambda_list(lambda_ind);
    
    %% earliest k scheme: wait for earliest k out of n

    % generate the pool for transmission time
    pool = exprnd(1/lambda,msglen,n)+C;

    for write = write_list
        endpt = 0;
        polygons = 0;
        response = 0;
        
        endpt_s = 0;
        polygons_s = 0;
        response_s = 0;
        
        for msgind = 1:msglen
            delayset = pool(msgind,:);
            [delay, delay_max, suc] = ext_UpdOrder(write,read,delayset);
%             [delay1, delay_max1, suc1] = ext_UpdOrderFix(write,read,delayset)
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
        avgage(write,lambda_ind) = polygons/endpt;
    end
    
    % approximation by log function and obtain optimal point    
    for write = write_list
        apprage(write,lambda_ind) = 1/lambda - log(1-write/n)/2/lambda + C/(write/n) + C/2;
    end
    alpha_opt(lambda_ind) = sqrt(4*lambda^2*C^2 + 8*lambda*C)/2 - lambda*C;

    
    %% selected k scheme: wait for fixed k out of n

    for write = write_list
        endpt = 0;
        polygons = 0;
        response = 0;
        for msgind = 1:msglen
            delayset = pool(msgind,:);
            [delay, delay_max, suc] = ext_UpdOrderFix(write,read,delayset);
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
    
    % evaluate selected k age formula
    for write = write_list
        apprage_fix(write,lambda_ind) = write/n*(C+1/lambda) + (n-write)/n*(C + 1/lambda + 1/lambda/write*(harmonic(write+1)-1)) ...
            + (n+2^(-write)*(n-write))/(n-2^(-write)*(n-write))/2 *(C+harmonic(write)/lambda+harmonic2(write)/(lambda^2*C+lambda*harmonic(write)));
    end
    
    %% selected k scheme with biase: premium group and free group

    for write = write_list
        endpt = 0;
        polygons = 0;
        response = 0;
        for msgind = 1:msglen
            delayset = pool(msgind,:);
            [delay, delay_max, suc] = ext_UpdOrderFixFree(write,read,delayset);
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
        avgage_fixfree(write,lambda_ind) = polygons/endpt;
    end
    
    for write = write_list
        endpt = 0;
        polygons = 0;
        response = 0;
        for msgind = 1:msglen
            delayset = pool(msgind,:);
            [delay, delay_max, suc] = ext_UpdOrderFixPrem(write,read,delayset);
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
        avgage_fixprem(write,lambda_ind) = polygons/endpt;
    end
   
    %% fixed redundancy scheme
    
    for write = write_list
        due = C+1/lambda*(harmonic(n)-harmonic(n-write));
        endpt = 0;
        polygons = 0;
        response = 0;
        for msgind = 1:msglen
            delay = pool(msgind,1);
            if delay > due % fail
                response = response + due;
            else
                polygons = polygons + 1/2*((response+delay)^2-delay^2);
                endpt = endpt + response;
                response = due;
            end        
        end
        avgage_fr(write,lambda_ind) = polygons/endpt;
    end
    
    for write = write_list
        due = C+1/lambda*(harmonic(n)-harmonic(n-write));
        apprage_fr(write,lambda_ind) = 1/(1-exp(-lambda*(due-C)))* (C-due*exp(-lambda*(due-C))) + 1/lambda + ...
            due/2* (1+exp(-lambda*(due-C)))/(1-exp(-lambda*(due-C)));
    
    end
    %% selective multicast (select k in advance and send to k nodes)
    %{
    for write = write_list
        sel_age(write,find(lambda==lambda_list)) = C+1/lambda+(2*n-write)/(2*write)*(C+harmonic(write)/lambda);
    end
    %}
end


%% generate the figure 1

figure(1)
set(gcf,'units','pixels','position',[10,10,400,250]);
hold on;
% color schemes
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
purple = [0.4940 0.1840 0.5560];

color = {blue, red, purple};
for i = 1:3
    %{
    % earliest k
    plot(write_list(2:2:100),avgage(2:2:100,i),'-','Color',color{i},'linewidth',1.5);    
    plot(write_list(4:4:100),apprage(4:4:100,i),'x','Color',color{i},'linewidth',1.5);
    plot(alpha_opt(i)*write,apprage(floor(alpha_opt(i)*write),i),'Color',color{i},'Marker','o','MarkerSize',10,'linewidth',2);
    %}
    % selected k
    l(i) = plot(write_list(2:2:100),avgage_fix(2:2:100,i),'-.','Color',color{i},'linewidth',1.5);
    % plot(write_list(4:4:100),apprage_fix(4:4:100,i),'s','Color',color{i},'linewidth',1.5);
    
    % selected k with free and premium group
    plot(write_list(2:2:100),avgage_fixfree(2:2:100,i),'-^','Color',color{i},'linewidth',1.5);
    plot(write_list(2:2:100),avgage_fixprem(2:2:100,i),'-v','Color',color{i},'linewidth',1.5);
end


xlabel('k','Fontsize',14,'FontName','Times');
% ylabel('average age \Delta_{(w,r=50)}','Fontsize',14,'FontName','Times');
ylabel('average age \Delta_{(k)}','Fontsize',14,'FontName','Times');
% title(['preemption at k out of n=' num2str(m)],'Fontsize',14,'FontName','Times');
leg = legend(l(1:3),{['\lambdac = ' num2str(lambda_list(1)*C)],['\lambdac = ' num2str(lambda_list(2)*C)],['\lambdac = ' num2str(lambda_list(3)*C)]},'location','North');
set(leg,'Fontsize',14,'FontName','Times');
axis([0 100 2 12]); 
grid on; box on;

%% generate the figure 2

figure(2);
set(gcf,'units','pixels','position',[10,10,400,250]);
hold on;
% color schemes
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
purple = [0.4940 0.1840 0.5560];
color = {blue, red, purple};
for i = 1:3
    % earliest k
    l(i) = plot(write_list(2:2:100),avgage(2:2:100,i),'-','Color',color{i},'linewidth',1.5);    
    % plot(write_list(4:4:100),apprage(4:4:100,i),'x','Color',color{i},'linewidth',1.5);
    plot(alpha_opt(i)*write,apprage(floor(alpha_opt(i)*write),i),'Color',color{i},'Marker','o','MarkerSize',10,'linewidth',2);
    % fixed redundancy
    plot(write_list(2:2:100),avgage_fr(2:2:100,i),'-.','Color',color{i},'linewidth',1.5);
    plot(write_list(2:2:100),apprage_fr(2:2:100,i),'+','Color',color{i},'linewidth',1.5);
    
    % selected k
    % plot(write_list(2:2:100),avgage_fix(2:2:100,i),'-.','Color',color{i},'linewidth',1.5);
    % plot(write_list(4:4:100),apprage_fix(4:4:100,i),'s','Color',color{i},'linewidth',1.5);   
    % selected k with free and premium group
    % plot(write_list(2:2:100),avgage_fixfree(2:2:100,i),'-^','Color',color{i},'linewidth',1.5);
    % plot(write_list(2:2:100),avgage_fixprem(2:2:100,i),'-v','Color',color{i},'linewidth',1.5);
    
end

xlabel('k','Fontsize',14,'FontName','Times');
% ylabel('average age \Delta_{(w,r=50)}','Fontsize',14,'FontName','Times');
ylabel('average age \Delta_{(k)}','Fontsize',14,'FontName','Times');
% title(['preemption at k out of n=' num2str(m)],'Fontsize',14,'FontName','Times');
leg = legend(l(1:3),{['\lambda = ' num2str(lambda_list(1))],['\lambda = ' num2str(lambda_list(2))],['\lambda = ' num2str(lambda_list(3))]},'location','North');
set(leg,'Fontsize',14,'FontName','Times');
axis([0 100 2 10]); 
grid on; box on;