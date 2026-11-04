function y=ch_tdfilt(fft0,fft1)

h_ph=angle(fft1);
h_mag=abs(fft1);
cos_sin=cos(h_ph)-j*sin(h_ph);
fft0=fft0.*cos_sin;
mult=fft0;
mult_real=real(mult);
mult_imag=imag(mult);
% subplot(2,1,1);plot(mult_real)
% subplot(2,1,2);plot(mult_imag)

Nend=2e6;
mult_real(find(real(mult)>Nend))=Nend;
mult_real(find(real(mult)<-Nend))=-Nend;
mult_imag(find(imag(mult)>Nend))=Nend;
mult_imag(find(imag(mult)<-Nend))=-Nend;
fft0=mult_real+j*mult_imag;

h_mag(find(h_mag<0.001))=0.001;
div=fft0./h_mag;
div_real=real(div);
div_imag=imag(div);
% subplot(2,1,1);plot(div_real)
% subplot(2,1,2);plot(div_imag)

Nend=2e5;
div_real(find(div_real>Nend))=0;
div_real(find(div_real<-Nend))=0;
div_imag(find(div_imag>Nend))=0;
div_imag(find(div_imag<-Nend))=0;

ifft_in=div_real+j*div_imag;
y=ifft(ifft_in,6144);



