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
msglen = 50000;

% read quorum
read = 1;


%% examine different lambda

avgage = zeros(length(write_list),3);
avgage_fix = zeros(length(write_list),3);
apprage = zeros(length(write_list),3);
apprage_fix = zeros(length(write_list),3);
alpha_opt = zeros(1,3);


for lambda_ind = 1:length(lambda_list)

    lambda = lambda_list(lambda_ind);
    
    %% preemptive multicast: true experiments obtained by a sample pool

    % generate the pool for transmission time
    poolsize = n*msglen;
    % pool = gen_update(delta,poolsize,poollen);
    % pool = k + nbinrnd(k,1-delta,poolsize,1);
    pool = exprnd(1/lambda,poolsize,1)+C;


    for write = write_list
        endpt = 0;
        polygons = 0;
        response = 0;
        for msgind = 1:msglen
            [delay, delay_max, suc] = ext_UpdOrder(write,read,n,pool);
%             [delay1, delay_max1, suc1] = ext_UpdOrderFix(write,read,n,pool);
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

    
    %% selected k scheme: wait for fixed k out of n

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
    
    % evaluate selected k age formula
    for write = write_list
        apprage_fix(write,lambda_ind) = write/n*(C+1/lambda) + (n-write)/n*(C + 1/lambda + 1/lambda/write*(harmonic(write+1)-1)) ...
            + (n+2^(-write)*(n-write))/(n-2^(-write)*(n-write))/2 *(C+harmonic(write)/lambda+harmonic2(write)/(lambda^2*C+lambda*harmonic(write)));
%         apprage_fix(write,lambda_ind) = C + 1/lambda + 1/lambda/write*(log(write+1)+double(eulergamma)-1) ...
%             + (n+2^(-write)*(n-write))/(n-2^(-write)*(n-write))/2 *(C+(log(write)+double(eulergamma))/lambda);

    end
    
    
    %% approximation by log function
    
    for write = write_list
        apprage(write,lambda_ind) = 1/lambda - log(1-write/n)/2/lambda + C/(write/n) + C/2;
    end

    alpha_opt(lambda_ind) = sqrt(4*lambda^2*C^2 + 8*lambda*C)/2 - lambda*C;
    
    %% evaluation of general age equation

    %{
    for write = write_list
        exos1(write) = C + (harmonic(n)-harmonic(n-write))/lambda;
        exos2(write) = C^2 + 2*C*(harmonic(n)-harmonic(n-write))/lambda + ((harmonic(n)-harmonic(n-write))^2+harmonic2(n)-harmonic2(n-write))/(lambda^2); % emit second order
    end    
    
    for write = write_list
        if n-write <= read
            pfail = 0;
        else
            pfail = nchoosek(n-write,read) / nchoosek(n,read);
        end
        condsum = 0;
        if n-write <= read
            condsum = mean(exos1(1:write));
        else
            for i=1:write
                condsum = condsum + exos1(i) * nchoosek(n-i,read-1)/(nchoosek(n,read)-nchoosek(n-write,read));
            end            
        end
        eqnage(write,find(lambda==lambda_list)) = condsum + (1+pfail)/2/(1-pfail)*exos2(write)/exos1(write);
    end
    %}
    
    %% approximation by stering approximation

    %{
    for write = write_list
        condsum = 0;
        if n-write <= read
            for i=1:n-read+1
                condsum = condsum + exos1(i) * (n-i)^(read-1)*read/(n^read);
            end  
        else
            for i=1:write
                condsum = condsum + exos1(i) * (n-i)^(read-1)*read/(n^read-(n-write)^read);
            end            
        end
        ster_age(write,find(lambda==lambda_list)) = condsum + (n^read+(n-write)^read)/2/(n^read-(n-write)^read)*exos2(write)/exos1(write);
    end
    
    [~, ster_opt(find(lambda==lambda_list))] = min(ster_age(:,find(lambda==lambda_list)));
    %}
    
   
    %% selective multicast (select k in advance and send to k nodes)
    
%     for write = write_list
%         sel_age(write,find(lambda==lambda_list)) = C+1/lambda+(2*n-write)/(2*write)*(C+harmonic(write)/lambda);
%     end
    
end


%% generate the figure

figure(1)
set(gcf,'units','pixels','position',[10,10,400,250]);
hold on;
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
purple = [0.4940    0.1840    0.5560];

color = {blue, red, purple};
for i = 1:3
    l(i) = plot(write_list(2:2:100),avgage(2:2:100,i),'-','Color',color{i},'linewidth',1.5);    
%     plot(narray(2:2:100),lambda_list(i)*eqnage(2:2:100,i),'o-','Color',color{i},'linewidth',1.5);
%     plot(write_list(2:2:100),lambda_list(i)*deepappr_age(2:2:100,i),'x','Color',color{i},'linewidth',1.5);
    plot(write_list(4:4:100),apprage(4:4:100,i),'x','Color',color{i},'linewidth',1.5);
    plot(alpha_opt(i)*write,apprage(floor(alpha_opt(i)*write),i),'Color',color{i},'Marker','o','MarkerSize',10,'linewidth',2);
%     plot(ster_opt(i),lambda_list(i)*ster_age(floor(ster_opt(i)),i),'Color',color{i},'Marker','o','MarkerSize',10,'linewidth',2);
    plot(write_list(2:2:100),avgage_fix(2:2:100,i),'-.','Color',color{i},'linewidth',1.5);
    plot(write_list(4:4:100),apprage_fix(4:4:100,i),'s','Color',color{i},'linewidth',1.5);
end


xlabel('k','Fontsize',14,'FontName','Times');
% ylabel('average age \Delta_{(w,r=50)}','Fontsize',14,'FontName','Times');
ylabel('average age \Delta_{(k)}','Fontsize',14,'FontName','Times');
% title(['preemption at k out of n=' num2str(m)],'Fontsize',14,'FontName','Times');
leg = legend(l(1:3),{['\lambdac = ' num2str(lambda_list(1)*C)],['\lambdac = ' num2str(lambda_list(2)*C)],['\lambdac = ' num2str(lambda_list(3)*C)]},'location','North');
set(leg,'Fontsize',14,'FontName','Times');
axis([0 100 2 12]); 
grid on; box on;