function [C,indice]=pa_somp(X,signal,T, T_0)

% M_SOMP for sparse signal recovery with partial prioir information T_0  (partial support), 2012/4/13, Linglong Dai
% T_0: known partial support
nbsignal=size(signal,2);
nbbase=size(X,2);
% lambdareg=1e-8;
% lambdareg=1e-4;
lambdareg= 0;  % standard SOMP

% R=signal;

% Phi=[];
% indice=[];

% Initialize the recovery signal
indice=T_0;
Phi=X(:,T_0);
R = signal - Phi * (Phi\signal);

for i=size(T_0,2)+1:T  
    C=zeros(i,nbsignal);
    [maxi,ind]=max(sum(abs(X'*R),2));
%     [maxi,ind]=max(sum(abs(X'*R),1));
    Phi=[Phi X(:,ind)];
    indice=[indice ind];
    regI=lambdareg*eye(i);
    for j=1:nbsignal
        C(:,j)=(Phi'*Phi + regI)\(Phi'*signal(:,j));
    end;
    A=Phi*C;
    R=signal-A;
    sum(sum(R.^2));
end;

Caux=zeros(size(X,2),size(signal,2));
Caux(indice,:)=C;
C=Caux;