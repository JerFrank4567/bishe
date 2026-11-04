function s = IHT(y,Phi,K)

% clear; clc
% 
% % 信号长度n， 测量次数m, 稀疏度s
% n  = 500;m  = 300;s  = 30;
% 
% x = 2*(rand(1,n)-0.5); 
% x_sparse = rand(1,n); 
% x(x_sparse<0.95)=0; 
% 
% %初始化
% X = randn(m,n); 
% Y=  X*x';
% 
% % Normalization and parameters.
% colnorm=mean((sum(X.^2)).^(.5));
% X=X./colnorm;
% Y=Y./colnorm;



% Normalization and parameters.
Y = y;
X = Phi;
colnorm=mean((sum(X.^2)).^(.5));
X=X./colnorm;
Y=Y./colnorm;

m = size(Phi,1);
n = size(Phi,2);

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

s = x1;
