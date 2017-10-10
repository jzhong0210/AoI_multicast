%% Simulation for m multi user case, shifted exponential distribution
% fix the number of users, vary the threshold

clear all;
close all;

% erasrue rate
lambda_list = [0.1 0.5 2];

% number of symbols in an update 
k = 1;
% k = 1;
% number of users
m = 100;
% m = 500;
% write_list = 1:m;
write_list = 10;
write = write_list;
read_list = 1:m;
% constant shift 
C = 2;
% C = 0.001;


% number of update messages
msglen = 2000;

% read quorum
% read = 10;


%% examine different lambda

normage_fr = zeros(length(write_list),3);
apprage = zeros(length(write_list),3);
alpha_opt = zeros(1,3);


for lambda = lambda_list

    %% n out of m (FR) with feedback (true experiments in a sample pool)

    % generate the pool for transmission time
    poolsize = msglen*m/100;
    % pool = gen_update(delta,poolsize,poollen);
    % pool = k + nbinrnd(k,1-delta,poolsize,1);
    pool_fr = exprnd(1/lambda,length(write_list)*msglen,1)+C;


    for read_th = read_list
        endpt = 0;
        polygons = 0;
        response = 0;
        for msgind = 1:msglen
            [n_fir, nmax_fir, suc_fir] = ext_UpdOrder(write,read_th,m,pool_fr);
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

        normage_fr(read_th,find(lambda==lambda_list)) = age_fir/k;
    end
 
     %% approximation by log function
% 
%     
%     for write = write_list
%         apprage(write,find(lambda==lambda_list)) = 1/lambda - log(1-write/m)/2/lambda + C/(write/m) + C/2;
%     end
% 
%     alpha_opt(find(lambda==lambda_list)) = sqrt(4*lambda^2*C^2 + 8*lambda*C)/2 - lambda*C;
    
    %% evaluation of general age equation
    
%     for write = write_list
%         exos1(write) = C + (harmonic(m)-harmonic(m-write))/lambda;
%         exos2(write) = C^2 + 2*C*(harmonic(m)-harmonic(m-write))/lambda + ((harmonic(m)-harmonic(m-write))^2+harmonic2(m)-harmonic2(m-write))/(lambda^2); % emit second order
%     end    
%     
%     for write = write_list
%         if m-write <= read
%             pfail = 0;
%         else
%             pfail = nchoosek(m-write,read) / nchoosek(m,read);
%         end
%         condsum = 0;
%         if m-write <= read
%             condsum = mean(exos1(1:write));
%         else
%             for i=1:write
%                 condsum = condsum + exos1(i) * nchoosek(m-i,read-1)/(nchoosek(m,read)-nchoosek(m-write,read));
%             end            
%         end
%         eqnage(write,find(lambda==lambda_list)) = condsum + (1+pfail)/2/(1-pfail)*exos2(write)/exos1(write);
%     end

    %% approximation by stering approximation

    for i = 1:m
        exos1(i) = C + (harmonic(m)-harmonic(m-i))/lambda;
        exos2(i) = C^2 + 2*C*(harmonic(m)-harmonic(m-i))/lambda + ((harmonic(m)-harmonic(m-i))^2+harmonic2(m)-harmonic2(m-i))/(lambda^2); % emit second order
    end  
    
    for read = read_list
        condsum = 0;
        if m-write <= read
            for i=1:m-read+1
                condsum = condsum + exos1(i) * (m-i)^(read-1)*read/(m^read);
            end  
        else
            for i=1:write
                condsum = condsum + exos1(i) * (m-i)^(read-1)*read/(m^read-(m-write)^read);
            end            
        end
        ster_age(read,find(lambda==lambda_list)) = condsum + (m^read+(m-write)^read)/2/(m^read-(m-write)^read)*exos2(write)/exos1(write);
    end
    
%     [~, ster_opt(find(lambda==lambda_list))] = min(ster_age(:,find(lambda==lambda_list)));

    %% further approximation

    
    for read = read_list
        alpha = write/m;
%         if m-write < read
%             deepappr_age(read,find(lambda==lambda_list)) =  - (2*(1-alpha)+1)*log(1-alpha)/(2*lambda) + ((C-1/(lambda*read))*(1+(1-alpha)^read)) + C/2;
%         else
%             deepappr_age(read,find(lambda==lambda_list)) = -1/(lambda*read) - (1+3*(1-alpha)^read)*log(1-alpha)/(2*lambda*(1-(1-alpha)^read)) + C + (C*(1+(1-alpha)^read))/(2*(1-(1-alpha)^read));         
%         end
        deepappr_age(read,find(lambda==lambda_list)) = -1/(lambda*read) - (1+3*(1-alpha)^read)*log(1-alpha)/(2*lambda*(1-(1-alpha)^read)) + C + (C*(1+(1-alpha)^read))/(2*(1-(1-alpha)^read));
    end
   
end


%% generate the figure

figure(1)
set(gcf,'units','pixels','position',[10,10,400,200]);
hold on;
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
purple = [0.4940    0.1840    0.5560];

color = {blue, red, purple};
for i = 1:3
    l(i) = plot(read_list,lambda_list(i)*normage_fr(:,i),'-','Color',color{i},'linewidth',1.5);    
%     plot(narray(2:2:100),lambda_list(i)*eqnage(2:2:100,i),'o-','Color',color{i},'linewidth',1.5);
    plot(read_list(2:2:100),lambda_list(i)*deepappr_age(2:2:100,i),'x','Color',color{i},'linewidth',1.5);
%     plot(narray(2:2:100),lambda_list(i)*apprage(2:2:100,i),'x-','Color',color{i},'linewidth',1.5);
%     plot(alpha_opt(i)*write,lambda_list(i)*apprage(floor(alpha_opt(i)*write),i),'Color',color{i},'Marker','o','MarkerSize',10,'linewidth',2);
%     plot(ster_opt(i),lambda_list(i)*ster_age(floor(ster_opt(i)),i),'Color',color{i},'Marker','o','MarkerSize',10,'linewidth',2);
end


xlabel('r','Fontsize',14,'FontName','Times');
ylabel('normalized age \lambda\Delta_{(w=10,r)}','Fontsize',14,'FontName','Times');
% title(['preemption at k out of n=' num2str(m)],'Fontsize',14,'FontName','Times');
leg = legend(l(3:-1:1),{['\lambda=' num2str(lambda_list(3))],['\lambda=' num2str(lambda_list(2))],['\lambda=' num2str(lambda_list(1))]},'location','Northeast');
set(leg,'Fontsize',14,'FontName','Times');
axis([0 100 0 12]); 
grid on; box on;