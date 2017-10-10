m = 100;
la = 0.5;
c = 2;
alpha = sqrt(2*la*c+la^2*c^2)-la*c;
age = zeros(m,1);
for n=1:m
    k = floor(n*alpha);
    h=zeros(k,1);
    for i=1:k
        h(i)=harmonic(i);
    age(n) = 1/k*sum(h)+(2*n-k)/(2*k)*harmonic(k);
    end
end
%%
figure(1)
set(gcf,'units','pixels','position',[10,10,400,250]);
hold on;

plot(1:m,age,'-','linewidth',1.5);    

xlabel('n','Fontsize',14,'FontName','Times');
ylabel('Minimum AoI','Fontsize',14,'FontName','Times');
% title(['preemption at k out of n=' num2str(m)],'Fontsize',14,'FontName','Times');
% set(leg,'Fontsize',14,'FontName','Times');
% axis([0 100 0 20]); 
grid on; box on;