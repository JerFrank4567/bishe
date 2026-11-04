function [C,indice]=somp(X,signal,T);

nbsignal=size(signal,2);
nbbase=size(X,2);
lambdareg=1e-8;

R=signal;

Phi=[];
indice=[];
for i=1:T  
    C=zeros(i,nbsignal);
    [maxi,ind]=max(sum(abs(X'*R),2));
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