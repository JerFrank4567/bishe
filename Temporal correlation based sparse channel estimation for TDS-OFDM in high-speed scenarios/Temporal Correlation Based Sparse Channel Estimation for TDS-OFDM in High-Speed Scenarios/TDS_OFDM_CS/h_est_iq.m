function [ha,h2]=h_est_iq(h0,h1,mod,mod_qpsk)
y_abs=zeros(1,length(h1));
y_re=abs(real(h1));y_im=abs(imag(h1));
abs1=y_re+y_im/2;
abs2=y_im+y_re/2;
index=abs1-abs2;
id1=find(index<0);
id2=find(index>=0);
y_abs(id1)=abs2(id1);
y_abs(id2)=abs1(id2);
hmax=max(y_abs);
    
% hmax=max(abs(h1));

if(mod_qpsk==1) 
   thr_h=0.1*hmax;   %qpsk
else
   thr_h=0.02*hmax;
end

index=find(y_abs>thr_h);
dif=zeros(1,length(h1));
dif(index)=h1(index)-h0(index);
    
if(mod==0)   %first iteration
    h2=h1+0.02*dif;
    ha=h1+0.04*dif;
else
    h2=h1+0.02*dif;
    ha=h1+0.04*dif;
end

if(mod_qpsk==1) 
    h2_thr = 0.1*max(abs(h2));
    index = find(abs(h2)<h2_thr);
    h2(index) = 0;
    index = find(abs(h2)>1);
    h2(index) = 1;
    ha_thr = 0.1*max(abs(ha));
    index = find(abs(ha)<ha_thr);
    ha(index) = 0;  %qpsk
    index = find(abs(ha)>1);
    ha(index) = 1;
else
    h2_thr = 0.05*max(abs(h2));
    index = find(abs(h2)<h2_thr);
    h2(index) = 0;
    index = find(abs(h2)>1);
    h2(index) = 1;
    ha_thr = 0.05*max(abs(ha));
    index = find(abs(ha)<ha_thr);
    ha(index) = 0;
    index = find(abs(ha)>1);
    ha(index) = 1;
end
    
        
    