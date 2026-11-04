%For Simulation, Increase T and repeat the channel. This is just for an
%example, not for complete curve. For HB-Kalman, when the dynamic channge
%is little, you can obtain good performance.
clc;

%% System Parameter
T = 150;
N = 150;
M = 100;
P1 = 0.05;
K = P1*N;
P10 = 0.01;
P01 = K*P10/(N-K);
sigma_h = 1;
sigma_n = 0.05;
sigma_a = 0.01;
Pth1 = 5*sigma_n/sqrt(M);
Pth2 = 5*sigma_n/sqrt(M);
OMP_Kmax = 20;
D = 1;
SNR = linspace(5, 20, 6);

%% Generate Phi
load pn_all.mat pn420;
PN_len = N;
PN255 = transpose(pn420(1,420-PN_len+1:420));
PN255 = sqrt(PN_len)*ifft(PN255);
PN255_matrix = zeros(PN_len, PN_len);
for i = 1:PN_len
    PN255_matrix(i:PN_len, i) = PN255(1:PN_len-i+1);
end
Phi_l = zeros(2*PN_len, PN_len);
for i = 1:PN_len
    Phi_l(i+PN_len:2*PN_len, i) = PN255(1:PN_len-i+1);
    Phi_l(i:i+PN_len-1,i) = PN255(1:PN_len);
end
Phi = Phi_l(N*2-M+1:N*2, :);
%random('Normal', 0, 1, M, N);
for i = 1:N
    Phi(:,i) = Phi(:,i)/norm(Phi(:,i));
end

%% Generate Channel
    s = zeros(N, T);
    a = zeros(N, T);
    h = zeros(N, T);
    w = zeros(N, T);
    s(:,1) = rand(1,N)<P1;
    a(:,1) = random('Normal', 0, sigma_h, N, 1);
    h(:,1) = s(:,1).*a(:,1);
    %figure;
    %stem(h(:,1));
    for i = 2:T
        K_pre = sum(s(:,i-1));

        for j = 1:N
            if (s(j,i-1)==1)
                if (K_pre>5)
                    s(j, i) = rand(1,1)>P10;
                else
                    s(j, i) = 1;
                end
            else
                if (K_pre<13) 
                    s(j, i) = rand(1,1)<P01;
                else
                    s(j, i) = 0;
                end
            end
        end

        %if (i==T/2)
        %    s(1:N/2,i) = rand(1,N/2)<P1;
        %end
        w(:, i) = random('Normal', 0, sigma_a, N, 1);
        a(:, i) = a(:, i-1) + w(:, i);
        h(:, i) = s(:, i).*a(:, i);
        %stem(h(:,i));
    end
    
mse_lce = zeros(length(SNR), T);
mse_omp = zeros(length(SNR), T);   
mse_asomp = zeros(length(SNR), T);
mse_Kalman = zeros(length(SNR), T);      
mse_somp = zeros(length(SNR), T);

for i_snr = 1:length(SNR)
    
    display(SNR(i_snr));
    %% Generate Y
    y = Phi*h;
    %y = y + random('Normal', 0, sigma_n, M, T);
    y = awgn(y, SNR(i_snr), 'measured');

    %% Linear Channel Estimation
    y_l = Phi_l*h;
    y_l = awgn(y_l, SNR(i_snr), 'measured');

    for i = 1:T
        y_t = y_l(PN_len+1:2*PN_len,i);
        h_lce = ifft(fft(y_t,PN_len)./fft(PN255,PN_len),PN_len);
        mse_lce(i_snr,i) = norm(h_lce-h(:,i))^2;
    end
    %figure;
    %stem(h_lce);

    %% Standard OMP

    for i = 1:T
        h_omp = OMP(N, Phi, y(:,i),17);
        mse_omp(i_snr,i) = norm(h_omp-h(:,i))^2;
    end

    %% A-SOMP
    R = 5;

    for i = 1: T
        if (i<R)
            y_i = mean(y(:,1:R), 2);
        else
            y_i = mean(y(:,i-R+1:i), 2);
        end
        h_asomp = OMP(N, Phi, y_i, 20);
        d_set = find(abs(h_asomp)>0);
        %h_asomp = zeros(N, 1);
        h_asomp(d_set) = pinv(Phi(:,d_set))*y(:,i);
        mse_asomp(i_snr,i) = norm(h_asomp-h(:,i))^2;    
    end

    %% Kalman sparse filter
    h_est_pre = zeros(N, T);
    h_est = zeros(N, T);
    mse_pre = zeros(N, N, T);
    mse = zeros(N, N, T);
    h_est(:, 1) = OMP(N, Phi, y(:,1),15);
    temp = zeros(1,N);
    d_set = find(s(:,1)>0);
    temp(d_set) = 1;
    mse(:, :, 1) = diag(temp)*sigma_h;
    Kmse = zeros(N, M, T);
    for i = 2:T
        h_est_pre(:,i) = h_est(:,i-1);
        temp = zeros(1,N);
        d_set = find(s(:,i)>0);
        temp(d_set) = 1;
        mse_pre(:,:,i) = mse(:,:,i-1) + diag(temp)*0.001;
        Kmse(:,:,i) = mse_pre(:,:,i)*Phi'*(eye(M)*0.0001+Phi*mse_pre(:,:,i)*Phi')^(-1);
        h_est(:,i) = h_est_pre(:,i) + Kmse(:,:,i)*(y(:,i)-Phi*h_est_pre(:,i));
        mse(:,:,i) = (eye(N)-Kmse(:,:,i)*Phi)*mse_pre(:,:,i);
    end

    for i = 1:T
        %h_mse(i) = norm(h(T(i,1:K(i)),i)-h_est(T(i,1:K(i)),i))^2/norm(h(T(i,1:K(i)),i))^2;
        mse_Kalman(i_snr,i) = norm(h(:,i)-h_est(:,i))^2;
    end


    %% S-OMP

    h_somp = zeros(N, 1);
    h_somp_pre = zeros(N, 1);
    h_add = zeros(N, 1);

    %figure;
    %stem(h(:,1));
    K_somp = zeros(1, T);
    for i = 1:T
        %if (i ==1)
            y_d = y(:, i) - Phi*h_somp_pre;
        %else
        %    h_t = zeros(N, 1);
        %    h_t(lambda) = pinv(Phi(:,lambda))*y(:,i);
        %    y_d = y(:,i) - Phi*h_t;
        %end
        [val, pos] = sort(abs(Phi'*y_d),'descend');
        lambda_add = pos(1:D);
        [val, pos] = sort(abs(h_somp_pre), 'descend');
        temp = sum(abs(h_somp_pre)>0);
        lambda = union(lambda_add, pos(1:temp));
        h_somp(lambda) = pinv(Phi(:, lambda))*y(:, i);
        h_somp(abs(h_somp)<Pth2) = 0;
        if (i == 1) % || norm(h_somp_pre-h_somp)>norm(h_somp))
            h_somp = OMP(N, Phi, y(:,i), OMP_Kmax);
        end
        lambda = find(abs(h_somp)>0);
        mse_somp(i_snr,i) = norm(h_somp-h(:,i))^2; 
        h_somp_pre = h_somp;
        K_somp(i) = sum(abs(h_somp)>0);
    end

end

%% Plot
%{
figure;
Ks = sum(s);
plot(1:T, Ks, 'b');
hold on;
plot(1:T, K_somp, 'g');
axis([1 T 0 15]);
legend('Real', 'S-OMP');
%}
figure;
semilogy(SNR, mean(mse_lce,2), 'k^-', 'LineWidth', 1.5);
hold on;
semilogy(SNR, mean(mse_omp,2), 'rv-', 'LineWidth', 1.5);
hold on;
semilogy(SNR, mean(mse_asomp,2), 'm<-', 'LineWidth', 1.5);
hold on;
semilogy(SNR, mean(mse_Kalman,2), 'c>-', 'LineWidth', 1.5);
hold on;
semilogy(SNR, mean(mse_somp,2), 'gs-', 'LineWidth', 1.5);
hold on;
grid on;
%axis([1 T 0.5e-3 0.5e-0]);
legend('Linear Method','Standard OMP', 'A-SOMP', 'HB-Kalman', 'Proposed S-OMP');
xlabel('Time');
ylabel('MSE');

%subplot(3,1,3);
%[val pos] = sort(abs(Ks(20:T)-K_somp(20:T)),'descend');
%stem(h(:,pos(1)+19));
%display(pos(1));





















