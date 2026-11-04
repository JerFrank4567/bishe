function y=ch_dlt(in0,in1);

fft0=fft(in0,2048);
fft1=fft(in1,2048);

mult=fft0.*fft1;
mult_real=real(mult);
mult_imag=imag(mult);
% subplot(2,1,1);plot(mult_real)
% subplot(2,1,2);plot(mult_imag)

Nend=1.5e6;
mult_real(find(real(mult)>Nend))=Nend;
mult_real(find(real(mult)<-Nend))=-Nend;
mult_imag(find(imag(mult)>Nend))=Nend;
mult_imag(find(imag(mult)<-Nend))=-Nend;
ifft_in=mult_real+j*mult_imag;

y=ifft(ifft_in,2048);
