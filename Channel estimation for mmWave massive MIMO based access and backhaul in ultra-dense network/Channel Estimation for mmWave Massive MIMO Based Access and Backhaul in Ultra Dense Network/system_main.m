clear;
% 仿真图表 1 10°（7）散射体，snr 128天线 vesus snr/ T时间开销√
% 128天线  T 固定 vesus snr/角度变大                 散射角度
% 
% 128天线  T 固定 10°vesus snr 和导频个数           导频个数
% 未来 反馈，导频污染方面                                              
%                                                   信道容量*
%                                                   收敛*
%                                                   压缩感知恢复概率*
% 是否可以写一个函数训练门限呢？
%% 基本参数
% s = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(s)
M = 128; % 信道状态信息矢量
S = 8;   % 稀疏度
%T = 30;  % 时间开销
P =1;  % 导频个数
sel_pilot = 3;% 1: ±1;2: 高斯 3; [0~2pi]均匀分布
sel_DOA_distribution = 2;% 1: DOA 180°均匀分布,support;2:DOA内集中分布; 3 DOA拉普拉斯分布
awgn_en = 1;
sel_CS = 3;% 1: G-SAMP;2 SAMP;3: LS;
step_size=1;
step_size_max = 20;
sim_N = 1e3;
snr_vector = [30];

T_set = [18];
MSE_record = zeros(numel(snr_vector),numel(T_set));
for snr_ii = 1:numel(snr_vector)
    for T_ii = 1:numel(T_set)
    T =    T_set(T_ii);
    mse_add = 0;
    snr = snr_vector(snr_ii);
    clearvars H_measurement;
    clearvars y_noise;clearvars y;
    
     H = zeros(M,P); % 设置的多维稀疏信道状态信息矢量
    
    H_T = dftmtx(M)/sqrt(M);
    if sel_pilot ==1
        H_P = randi(2,[T,M,P])*2-3;
    elseif sel_pilot ==2
        H_P = (randn(T,M,P)+randn(T,M,P)*1i)/sqrt(2);
    elseif sel_pilot ==3
        H_P = exp(1i*2*pi*rand(T,M,P));
    end
    
    for sim_ii = 1:sim_N
    %% 观测矩阵 傅里叶阵和导频矩阵
   
    %% CSI 稀疏信号生成
    h = zeros(M,P);
    support = [];
    if sel_DOA_distribution==1 %DOA 180°均匀分布,support
        while numel(support)~=S
            support = randi(M,[S,1]);% DOA 角度集合
            support = unique(support);
        end
        h(support,:) = (randn(S,P)+randn(S,P)*1i)/sqrt(2);
    elseif sel_DOA_distribution==2
         support1 = randi(M,[1,1]);
         support = [support1:support1+S-1];
         support(find(support>M)) = support(find(support>M))-M;
         h(support,:) = (randn(S,P)+randn(S,P)*1i)/sqrt(2);
    elseif sel_DOA_distribution==3
    end
    for  pp=1:P
        h(:,pp) = h(:,pp)/norm(h(:,pp));
    end
    %% 信号接收
    for pp = 1:P
        H_measurement(:,:,pp) = H_P(:,:,pp)*H_T;
%         y(:,pp) = H_measurement(:,:,pp)*h(:,pp);
     if awgn_en==1
           % y_noise(:,pp) = y(:,pp)+randn(T,P)/sqrt(10^(snr/10));%awgn(y,snr,'measured');
            y_noise(:,pp) =  H_measurement(:,:,pp)*h(:,pp)+randn(T,1)/sqrt(10^(snr/10));
%             y_noise(:,pp) =  H_measurement(:,:,pp)*(h(:,pp)+randn(M,1)/sqrt(10^(snr/10)));
        else
            y_noise(:,pp) =  H_measurement(:,:,pp)*h(:,pp);
    end    
    end
  
    if sel_CS ==1
         if snr==10  
             threshold = 0.1;
         elseif snr==20
              threshold = 0.075;
         elseif snr==30
              threshold = 0.035;
         end
         [x_est support_est support_num] = SAMP_distributed(y_noise, H_measurement, step_size, threshold,20,0);
        
    elseif sel_CS==2
    elseif sel_CS==3
        x_est = zeros(size(h));
        for pp = 1:P
            x_est(support,pp) =  pinv(H_measurement(:,support,pp))*y_noise(:,pp);
        end
    end

%      figure(1);stem(abs(x_est(:,1)));hold on;stem(abs(h(:,1)),'rx');hold off;legend('Est','True')
%     figure(2);stem(abs(x_est(:,2)));hold on;stem(abs(h(:,2)),'rx');hold off;legend('Est','True');
    
    mse_add = mse_add + norm(x_est(1:end)-h(1:end))^2;
   fprintf('sim_num = %d, snr=%d, T_overhead = %d, mse = %7.9f,mse_s = %7.9f \n',sim_ii,snr,T,mse_add/sim_ii/P,norm(x_est(1:end)-h(1:end))^2/P);
    
    end
    MSE_record(snr_ii,T_ii) = mse_add/sim_ii/P;
    end
end

save LSvsT1_1_30P64snr10dB.mat ;