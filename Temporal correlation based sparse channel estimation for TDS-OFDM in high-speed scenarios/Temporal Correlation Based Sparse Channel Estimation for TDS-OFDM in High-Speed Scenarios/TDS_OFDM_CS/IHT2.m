clear; clc

% 信号长度n， 测量次数m, 稀疏度s
n  = 500;m  = 300;s  = 30;

x = 2*(rand(1,n)-0.5); 
x_sparse = rand(1,n); 
x(x_sparse<0.95)=0; 

%初始化
X = randn(m,n); 
Y=  X*x';
snr_db = 80;
Y = Y + 10^(-snr_db/10)*1*(randn(length(Y),1)+j*randn(length(Y),1))/sqrt(2);

% Normalization and parameters.
colnorm=mean((sum(X.^2)).^(.5));
X=X./colnorm;
Y=Y./colnorm;

delta=m/n;
nsweep=m;
tol=1e-6;
xinitial=zeros(n,1);
far=.04; 
relaxation_p=.65;   % Optimal value of relaxation.
it = 0;             % 计算迭代次数

x1 = xinitial;
lambda =  norminv(1-far/2,0,1);

for sweep=1:nsweep,
    % estimate noise
    r        = Y - X*x1;
    c        = abs(X'*r);
    sigmaest = relaxation_p*median(c)/.6745;
    %
    x2 = x1 + relaxation_p*X'*r;
    x1 = (abs(x2) >= lambda.*sigmaest) .* x2;
    it = it +1;
    if norm(Y-X*x1)/norm(Y)< tol
        break
    end
end

disp(it)
stem(x1,'LineStyle','none')  %重构信号，篮圈
hold on
stem(x,'r','.','LineStyle','none')   %原始的信号，红点
title('recover')
legend('Recovered signal', 'Original signal')