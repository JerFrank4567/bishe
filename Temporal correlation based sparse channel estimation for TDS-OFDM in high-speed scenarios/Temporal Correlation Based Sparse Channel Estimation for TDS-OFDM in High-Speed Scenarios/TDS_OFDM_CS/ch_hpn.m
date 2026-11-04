function y=ch_hpn(in0,in1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nh=2048;
%%test fft_peng
% in0=in0(1:nh)+in0(nh+1:2*nh);
in0=[in0 zeros(1,2048-length(in0))];
in1=[in1 zeros(1,2048-length(in1))];
fft0=fft(in0);
fft1=fft(in1);
% %phase codic
h_ph=angle(fft1);
h_mag=abs(fft1);
cos_sin=cos(h_ph)-j*sin(h_ph);
fft0=fft0.*cos_sin;
h_mag(find(h_mag<0.001))=0.001;
div=fft0./h_mag;
div_real=real(div);
div_imag=imag(div);
Nend=20;
% subplot(2,1,1);plot(div_real)
% subplot(2,1,2);plot(div_imag)

div_real(find(div_real>Nend))=0;
div_real(find(div_real<-Nend))=0;
div_imag(find(div_imag>Nend))=0;
div_imag(find(div_imag<-Nend))=0;
ifft_in=div_real+j*div_imag;
y=ifft(ifft_in,nh);
y=y(1:420);
