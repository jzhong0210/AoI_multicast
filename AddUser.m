clear all;
n = 100;
lambda = 1;
c = 1;
alpha = sqrt(2*lambda*c+lambda^2*c^2)-lambda*c;

%% main

age = zeros(n,1);
for j = 1:n
    k = round(j*alpha);
    h = zeros(k,1);
    for i = 1:k
        h(i) = c+1/lambda*(harmonic(j)-harmonic(j-i));
    end
    age(j) = 1/k*sum(h)+(2*j-k)/(2*k)*h(k);
end

%% plot
figure(1)
set(gcf,'units','pixels','position',[10,10,400,250]);
hold on;
plot(1:n,age,'-','linewidth',1.5);
% plot(1:n,)
xlabel('n','Fontsize',14,'FontName','Times');
ylabel('Minimum AoI','Fontsize',14,'FontName','Times');
% set(leg,'Fontsize',14,'FontName','Times');
axis([1 100 2.8 4]); 
grid on; box on;