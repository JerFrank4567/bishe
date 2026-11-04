function out=mmse(in,weigh)
%in=randint(1,420,[0 2048])+j*randint(1,420,[0 2048]);weigh=1;
nh=length(in);
% weigh=1;
if(weigh==0)alfa=0.9;
elseif(weigh==1) alfa=0.99;
elseif(weigh==2) alfa=0.998;
else alfa=0.9995;
end

real_in=real(in);
imag_in=imag(in);
hmax=max([abs(real_in); abs(imag_in)]);
real_index=find(abs(real_in)<hmax*0.1);
imag_index=find(abs(imag_in)<hmax*0.1);
real_fa=ones(nh,1);
imag_fa=ones(nh,1);
real_fa(real_index)=1./(alfa+4*(1-alfa)*(hmax/2./real_in(real_index)).^2);
imag_fa(imag_index)=1./(alfa+4*(1-alfa)*(hmax/2./imag_in(imag_index)).^2);
real_in=real_in.*real_fa;
imag_in=imag_in.*imag_fa;
real_index=find(abs(real_in)<hmax*0.05);
imag_index=find(abs(imag_in)<hmax*0.05);
real_in(real_index)=0;
imag_in(imag_index)=0;%%Ê±ÓòÂË²¨
out=(real_in+j*imag_in);
out = out(1:420);
% out = out./norm(out);
% out = out*1.1;
