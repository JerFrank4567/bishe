% 1) TFT-OFDM without channel coding: Preamble + Sub_frame * (OFDM (3780) + Single PN (420))
%    confige ZP or PN by the parameter PN_ZP, and increase the pilot group number G>L to solve the ill-posed problem
% 2) Compressive Sensing for channel estimation based on IBI-free region
% 3) CS recovery algorithm with paritial priori information of the channel
% 4) When 2L-1<=M, using last >=L samples without IBI to recovery channel via LS

clear;clc;close all;
warning off all;
Global_Var_sys;
PN_ZP = 1;   % 1: PN; 0:ZP
PN_TF = 1;   % 1: PN in time domain (T); 0: PN in frequency domain (F)
MultipathType = 2;
Sim_frame_num = 60;
Sub_frame = 2;
num_obs = 5:10:60;
Mod_order =[2];    % 8: 256QAM; 6:64QQM; 4: 16QAM; 2:QPSK; 1: BPSK
SNR = linspace(20,20,1);  %dB
weigh_acq = 3; % 0 1 2 3 for alpha = 0.9 0.99 0.998 0.9995
Chan_Len = PN_len - 150;

load pn_all.mat pn420;
% PN256=transpose(pn420(1,:));
PN256=transpose(pn420(1,166:420));
if(PN_TF == 1)     % PN_time(+1,-1)
%     PN_offset =(22.6274-1)/PN_len;  % for PN511
%     PN_offset =(16-1)/PN_len;       % for PN255
    PN_offset =0/PN_len; 
    PN256 = PN256 + PN_offset;
elseif(PN_TF == 0) % PN-MC
    PN256=sqrt(length(PN256)) * ifft(PN256);
end

% % CAZAC
% PN256 = zeros(PN_len,1);
% for k = 1: PN_len
%     PN256(k) = exp(sqrt(-1)*(2*pi/PN_len)*(PN_len-1)*(k^2/2));
% end

% Sensing matrix
PN256_matrix = zeros(PN_len,PN_len);
for i = 1: PN_len
    PN256_matrix(i:PN_len,i) = PN256(1:PN_len-i+1);
end

if PN_ZP ==1
   ZP256 = PN256;
else 
   ZP256 = zeros(PN_len,1); 
end

awgnSignalPower = 10*log10((PN_total_Len*PN_ZP+Idft_len*1)/(PN_total_Len+Idft_len));

Stored_channel_gains = zeros(length(Rayleigh_path),Sim_frame_num * Sub_frame);  % Restore block fading channel
for k = 1: Sim_frame_num * Sub_frame   
  Stored_channel_gains(:,k)= RayleighChannel.PathGains;
  tx = randint(Frame_len, 1, 2);       % SuperBlock Fading
  dpskSig = dpskmod(tx, 2);      
  y = filter(RayleighChannel, dpskSig);          % Run signal through channel
end

State_srrc_trans = [];
State_srrc_receive = [];

CE = zeros(Sim_frame_num, PN_len);
storage = zeros(Sim_frame_num, 1);

Send_data_aft_srrc_buf = [];
Send_data_aft_srrc_down_buf = [];
Receive_data_aft_srrc_buf = [];
Receive_data_aft_srrc_down_buf = zeros(2*Frame_len, 1);
Fnum_demod_buf = [];
PN_Selected_d1 = zeros(PN_len,1);

State_channel = [];
BER = zeros(1,length(SNR));
MSE_channel = zeros(Sim_frame_num * Sub_frame,length(num_obs));
MSE_channel_new = zeros(Sim_frame_num * Sub_frame,length(num_obs));
MSE_channel_new_ML = zeros(Sim_frame_num * Sub_frame,length(num_obs));

Probability_correct_recovery = zeros(Sim_frame_num*Sub_frame,length(num_obs));
Probability_correct_recovery_new = zeros(Sim_frame_num*Sub_frame,length(num_obs));
Probability_correct_recovery_new_ML = zeros(Sim_frame_num*Sub_frame,length(num_obs));
BER_store = zeros(Sub_frame,Sim_frame_num);

% Assume perfect delay estimation at first
if MultipathType == 2    
  Tao_l = transpose(Rayleigh_path)-1;   % To ensure perfect delay estimation, the power delay of the channel model can be adjusted to be just intergeral times of the sampling rate
elseif MultipathType == 1             
  Tao_l = transpose(Static_path)-1;  
elseif MultipathType == 0             
  Tao_l = 1;  
end 

c = 0;                                  % (2*Pilot_span+1) pilots in a small group;
% Pilot_group = length(Static_path)+4;   
Pilot_group = 0;                                
Pilot_number = (2 * c + 1) * Pilot_group;
Pilot_amptitude = 1;                   % 1: NO Power boosting --> has no impact on the joint CE performance; 10: power boosting

% Pilot_distance = round(Idft_len/Pilot_group) - 1;     % pilot distance
% Pilot_offset = floor(Pilot_distance/2);               % pilot offset
% Pilot_group_center =  Pilot_offset + (0:1:Pilot_group - 1).* Pilot_distance;

Pilot_group_center = DIG(Idft_len,Pilot_group); %  Generate Pilot_group different integers in the range of [1:Idft_len]
% Pilot_group_center = 60 + 1:round(Idft_len/Pilot_group):Idft_len;
% Pilot_group_center = [(1:18) (Idft_len-17:Idft_len)];
%  Pilot_group_center = [88 115 320  546 793  831 1080  1174 1201 1396 1417 1625 1732 1756  1846  2501  2635 2921 2969  3459];  % 20 pilots
% Pilot_group_center = [71  693  935  1018  1257 1316 1608 1856  2239  2327  2393 2508  3304 3372 3397]; % 15 pilots

IDFT_block = zeros(Idft_len,1);
for center_index = 1:length(Pilot_group_center)
     IDFT_block(Pilot_group_center(center_index) - c: Pilot_group_center(center_index) + c) = 1;
end
Pilot_index = find(IDFT_block==1);
Data_index = find(IDFT_block~=1);
IDFT_block_freq = zeros(Idft_len,1);
IDFT_block_freq(Pilot_index) = 0;  % configure 2c zeros pilots around the central pilot
IDFT_block_freq(Pilot_group_center) = Pilot_amptitude; 
Central_Pilot_freq = IDFT_block_freq(Pilot_group_center);

W = exp(-sqrt(-1)*2*pi/Idft_len);
W_matrix =  zeros(length(Pilot_group_center),length(Tao_l));
for row = 1: length(Pilot_group_center)
    for column = 1:length(Tao_l);
        W_matrix(row,column) = W.^((Pilot_group_center(row)-1).*Tao_l(column));
    end
end
cond(W_matrix)
Moore_Penrose_matrix = pinv(W_matrix);  % reverse matrix 

PN256_matrix_direct = PN256_matrix(end-90+1:end,Tao_l + 1);
cond(PN256_matrix_direct)
Moore_Penrose_PN256_matrix = pinv(PN256_matrix_direct);  % reverse matrix 

for m_index = 1: length(num_obs)

 for ii = 1:length(Mod_order)
    
    if (Mod_order(ii) == 1)
        Mod_power = 1;
    elseif (Mod_order(ii) == 2)
        Mod_power = 1*2;
    elseif (Mod_order(ii) == 4)
        Mod_power = 5*2;
    elseif (Mod_order(ii) == 6)
        Mod_power = 21*2;
    elseif (Mod_order(ii) == 8)
        Mod_power = 75*2;
    end
    mod_QPSK = (Mod_order(ii) == 2);
    
  for i_snr = 1:length(SNR)
      
        reset(RayleighChannel,[10; 20]);  % reset Rayleigh channels to ensure the same channel under different SNR. 2010.3.11
        
        Send_data_buf = zeros(Sub_frame*(Idft_len-Pilot_number)*Mod_order(ii), Sim_frame_num);
        Receive_data_buf = zeros(Sub_frame*(Idft_len-Pilot_number)*Mod_order(ii), Sim_frame_num);

   for i_frame = 1:Sim_frame_num
%% Transmitter    
    Payload_integer = randint((Idft_len-Pilot_number), Sub_frame,2^Mod_order(ii));
    IDFT_block_time = zeros(Idft_len, Sub_frame);
    for i_subframe = 1:Sub_frame
        Payload_bit = dec2bin(Payload_integer(:,i_subframe), Mod_order(ii));
        Payload_bit = Payload_bit(:);
        Send_data_buf((i_subframe-1)*(Idft_len-Pilot_number)*Mod_order(ii)+1: i_subframe*(Idft_len-Pilot_number)*Mod_order(ii),i_frame) = Payload_bit; % storage for BER 
        Mod_symbol_freq = qammod(Payload_integer(:,i_subframe),2^Mod_order(ii),0,'gray')/sqrt(Mod_power);  % normalize the power
        IDFT_block_freq(Data_index) = Mod_symbol_freq;
        IDFT_block_time(:,i_subframe) = sqrt(Idft_len) * ifft(IDFT_block_freq, Idft_len);
    end
    Preamble = [PN256(2:end); PN256];     % cyclic extension of length M-1, and the preamble length is M
    TS_matrx = repmat(ZP256, 1, Sub_frame);
    Tranmistted_frames = [IDFT_block_time; TS_matrx];
    Send_data_bef_srrc_down = [Preamble; ZP256; Tranmistted_frames(:)];  % Preamble + ZP + Sub_frame * (OFDM (3780) + Single PN or ZP (420)) 
%% Channel   
   if MultipathType == 2    
%       Ideal_channel = RayleighChannel.PathGains;
%       Send_data_channel = filter(RayleighChannel, Send_data_bef_srrc_down);  
         Send_data_channel = zeros(length(Send_data_bef_srrc_down),1);
         for i_subframe = 0:Sub_frame
              Block_Fading_multipath_channel = zeros(PN_len,1);
              if i_subframe == 0
                 Block_Fading_multipath_channel(Rayleigh_path) = Stored_channel_gains(:,(i_frame-1)*Sub_frame+1); % Load block fading channels
                 Send_data_subframe = Send_data_bef_srrc_down(1:3 * PN_len - 1);
                 [Send_data_channel(1:3 * PN_len - 1) State_channel] = filter(Block_Fading_multipath_channel,1,Send_data_subframe,State_channel); 
              else
                 Block_Fading_multipath_channel(Rayleigh_path) = Stored_channel_gains(:,(i_frame-1)*Sub_frame + i_subframe); % Load block fading channels
                 Send_data_subframe = Send_data_bef_srrc_down((i_subframe -1) * Frame_len + 3 * PN_len: i_subframe * Frame_len + 3 * PN_len -1);
                 [Send_data_channel((i_subframe -1) * Frame_len + 3 * PN_len: i_subframe * Frame_len + 3 * PN_len -1)  State_channel] = filter(Block_Fading_multipath_channel,1,Send_data_subframe,State_channel); 
              end
         end
    elseif MultipathType == 1             
      [Send_data_channel State_channel] = filter(Static_multipath_channel,1,Send_data_bef_srrc_down,State_channel);
    elseif MultipathType == 0
      Send_data_channel = Send_data_bef_srrc_down;
   end 
   Send_data_channel = awgn(Send_data_channel,SNR(i_snr),awgnSignalPower);  % don't use 'Measured' to guarantee the identical noise power
   
%% Receiver        
    Receive_data = Send_data_channel;
    
    % MMV
    Received_frames = Receive_data(3 * PN_len:end);
    Received_signal = reshape(Received_frames, Frame_len, Sub_frame);  
    M = num_obs (m_index); % number of measurements (much less observation number is required for modified CoSaMP: 20 VS. 60)
    Index_sel = PN_len-M+1:1:PN_len;  % ony the last M symbols for random selection;
    X = PN256_matrix(Index_sel,:);                    %  sensing matrix
    signal = Received_signal(Frame_len - M +1:Frame_len,:);
    T = 6;
    [C,indiceOMP]=somp(X,signal,T);
    
    % Averaging among several received TSs to improve path delay estimation accuracy
    PN_in_current = PN256;
    Averaged_Ch_est_direct_correlation = zeros(PN_len,1);
    Averaged_subframe = Sub_frame;  % number of subframe for averaging
    temp = Receive_data(end-Averaged_subframe*Frame_len+1-PN_len:end-PN_len);
    Received_matrix = reshape(temp, Frame_len, Averaged_subframe);
    if Sub_frame > 1
        Received_matrix_mean = transpose(mean(transpose(Received_matrix)));
        Averaged_Received_PN_with_Interference = Received_matrix_mean(1:PN_len);
    else
        Averaged_Received_PN_with_Interference = temp(1:PN_len);
    end
    for aa =1:PN_len
        Averaged_Ch_est_direct_correlation(aa) = (1/PN_len)*sum(Averaged_Received_PN_with_Interference.*conj([PN_in_current(PN_len-aa+2:PN_len); PN_in_current(1:PN_len-aa+1)]));
    end
    Th = 0.1; % The threshold should be optimized
    T_0 = find(abs(Averaged_Ch_est_direct_correlation)> Th);  % Most significant taps
    T_0 = T_0';
    
    % PA-SOMP
    [C_new,indiceOMP_new]=pa_somp(X,signal,T,T_0);
    
    % Step 3: ML path coefficient estimation
    C_new_ML = zeros(PN_len, Sub_frame);
    support = find(abs(C_new(:,1)) > 0);
     for i_subframe = 1: Sub_frame       
        Timing_offset = 0;
        Frame_current = Receive_data((i_subframe -1) * Frame_len + 2 * PN_len - Timing_offset: i_subframe * Frame_len + 3 * PN_len -1 - Timing_offset);
        y =  Frame_current(Index_sel + Frame_len);
%         Num = size (support,1) + 10;
        Num = M;
        Row_sel = randperm(M);
%         rows = sort(Row_sel(M-Num+1:M),'ascend');
        rows = M -Num + 1: M;
        Phi_for_ML = X(rows,support);
        Path_coefficient = pinv(Phi_for_ML) * y(rows);
        C_new_ML(support,i_subframe) = Path_coefficient;
     end
    
     % MSE calculation
     Estimated_CIR_Group = C;
     Actual_Channel = zeros(PN_len,Sub_frame);
     Actual_Channel(Rayleigh_path,:) = Stored_channel_gains(:,(i_frame-1)*Sub_frame +1:i_frame*Sub_frame);
     Th_mse = 10^(-2);
     for i_subframe = 1: Sub_frame
         if MultipathType == 2  
            error = sum(abs(Estimated_CIR_Group(:,i_subframe) -  Actual_Channel(:,i_subframe)).^2)/sum(abs(Actual_Channel(:,i_subframe)).^2);
            MSE_channel((i_frame-1)*Sub_frame + i_subframe,m_index) = error;
            if error < Th_mse
                Probability_correct_recovery((i_frame-1)*Sub_frame + i_subframe,m_index) = 1;
            end
         else
            error = sum(abs(Estimated_CIR_Group(:,i_subframe) -  Static_multipath_channel).^2)/sum(abs(Static_multipath_channel).^2);
            MSE_channel((i_frame-1)*Sub_frame + i_subframe,m_index) = error;
            if error < Th_mse
                Probability_correct_recovery((i_frame-1)*Sub_frame + i_subframe,m_index) = 1;
            end
         end
     end  
     
     % MSE_new calculation
     Estimated_CIR_Group = C_new;
     Th_mse = 10^(-2);
     for i_subframe = 1: Sub_frame
         if MultipathType == 2  
            error = sum(abs(Estimated_CIR_Group(:,i_subframe) -  Actual_Channel(:,i_subframe)).^2)/sum(abs(Actual_Channel(:,i_subframe)).^2);
            MSE_channel_new((i_frame-1)*Sub_frame + i_subframe,m_index) = error;
            if error < Th_mse
                Probability_correct_recovery_new((i_frame-1)*Sub_frame + i_subframe,m_index) = 1;
            end
         else
            error = sum(abs(Estimated_CIR_Group(:,i_subframe) -  Static_multipath_channel).^2)/sum(abs(Static_multipath_channel).^2);
            MSE_channel_new((i_frame-1)*Sub_frame + i_subframe,m_index) = error;
            if error < Th_mse
                Probability_correct_recovery_new((i_frame-1)*Sub_frame + i_subframe,m_index) = 1;
            end
         end
     end
     
     % MSE_new_ML calculation
     Estimated_CIR_Group = C_new_ML;
     Th_mse = 10^(-2);
     for i_subframe = 1: Sub_frame
         if MultipathType == 2  
            error = sum(abs(Estimated_CIR_Group(:,i_subframe) -  Actual_Channel(:,i_subframe)).^2)/sum(abs(Actual_Channel(:,i_subframe)).^2);
            MSE_channel_new_ML((i_frame-1)*Sub_frame + i_subframe,m_index) = error;
            if error < Th_mse
                Probability_correct_recovery_new_ML((i_frame-1)*Sub_frame + i_subframe,m_index) = 1;
            end
         else
            error = sum(abs(Estimated_CIR_Group(:,i_subframe) -  Static_multipath_channel).^2)/sum(abs(Static_multipath_channel).^2);
            MSE_channel_new_ML((i_frame-1)*Sub_frame + i_subframe,m_index) = error;
            if error < Th_mse
                Probability_correct_recovery_new_ML((i_frame-1)*Sub_frame + i_subframe,m_index) = 1;
            end
         end
     end
     
     
     
     
     
     
%      
%      
%      
%      
%    
%     for i_subframe = 1: Sub_frame
%         
%         fprintf('Observation number  = %d, i_frame = %d\n', num_obs (m_index),i_frame);
%         
%         Timing_offset = 0;
%         Frame_current = Receive_data((i_subframe -1) * Frame_len + 2 * PN_len - Timing_offset: i_subframe * Frame_len + 3 * PN_len -1 - Timing_offset);
%         
%         % Initial channel estimation only at the begining of each superframe
%         if i_subframe == 1
%             PN_out_current = Receive_data(PN_len: 2 * PN_len -1);
%             PN_in_current = PN256;
%             Ch_est_preamble = ifft(fft(PN_out_current)./fft(PN_in_current));
% %             Ch_est_unnoised_preamble = mmse(Ch_est_preamble,weigh_acq);
%             Ch_est_unnoised_preamble = Ch_est_preamble;   
%         end
%         
%         % Rough channel estimation to acquire path delays
%         if i_subframe < Sub_frame-1
%             Ch_est_direct_correlation = zeros(PN_len,1);
%             Data_for_correlation = Receive_data((i_subframe -1) * Frame_len + 2 * PN_len: (i_subframe +1) * Frame_len + 3 * PN_len -1);
%     %         Received_PN_with_Interference = Data_for_correlation(end-PN_len+1:end);
%             Received_PN_with_Interference = (Data_for_correlation(end-PN_len-Frame_len+1:end-Frame_len) + Data_for_correlation(end-PN_len+1:end))/2; % two received TSs
%             for aa =1:PN_len
%                 Ch_est_direct_correlation(aa) = (1/PN_len)*sum(Received_PN_with_Interference.*conj([PN_in_current(PN_len-aa+2:PN_len); PN_in_current(1:PN_len-aa+1)]));
%             end
%         end
%         
%         % Averaging among several received TSs to improve path delay estimation accuracy
%         if i_subframe < Sub_frame-1   
%             Averaged_Ch_est_direct_correlation = zeros(PN_len,1);
%             Averaged_subframe = 5;  % number of subframe for averaging
%             temp = Receive_data(end-Averaged_subframe*Frame_len+1-PN_len:end-PN_len);
%             Received_matrix = reshape(temp, Frame_len, Averaged_subframe);
%             Received_matrix_mean = transpose(mean(transpose(Received_matrix)));
%             Averaged_Received_PN_with_Interference = Received_matrix_mean(1:PN_len);
%             for aa =1:PN_len
%                 Averaged_Ch_est_direct_correlation(aa) = (1/PN_len)*sum(Averaged_Received_PN_with_Interference.*conj([PN_in_current(PN_len-aa+2:PN_len); PN_in_current(1:PN_len-aa+1)]));
%             end
%         end
%         
%         % LS estimation when 2L-1<=M (using last >=L samples without IBI to recovery channel via LS); For simplicity, use the second half
%         
% %        figure
% %        stem(abs(Ch_est_direct_correlation),'bo','LineWidth',1.0);
% %        hold on;
% %        stem(abs(Averaged_Ch_est_direct_correlation),'ro','LineWidth',1.0);
% %        stem(abs(Static_multipath_channel),'k^','LineWidth',1.0);
% % %        axis([0 max(Tao_l)+10 0 1]);
% %        axis([0 max(Static_path)+10 0 1]);
% % %        legend('Averaged correlation','Ideal Channel')
% %        legend('Direct Correlation','Averaged correlation', 'Ideal Channel')
% %        hold off
%         
% %         % Channel predition
% %         hm1 = Ch_est_unnoised_1_current;
% %         if (i_subframe > 1)
% %             h0 = Ch_est_prv;
% %         else
% %             h0 = Ch_est_unnoised_1_current;
% %         end
% %         [Ch_est_middle Ch_est_end] = h_est_iq(hm1,h0,0,mod_QPSK);  % 
% 
%         %%% Compressive sensing based channel estimation
%         M = num_obs (m_index); % number of measurements (much less observation number is required for modified CoSaMP: 20 VS. 60)
% %         Index_sel = DIG(PN_len - max(Tao_l),M) + max(Tao_l);  % maximum random selection;
%         Index_sel = PN_len-M+1:1:PN_len;  % ony the last M symbols for random selection;
%         y = (Frame_current(Index_sel) + Frame_current(Index_sel + Frame_len))/2;
% 
%         %  OMP for signal recovery
% 
%         % Initialization
%         N = PN_len;
%         K = length(Tao_l);
%         Phi = PN256_matrix(Index_sel,:);                    %  sensing matrix
%         hat_x = zeros(N,1);                                 %  estimated signal                     
%         Aug_t = [];                                         %  argment matrix (empty at first)
%         r_n = y;                                            %  residual signal
%         
% %         % Implementation 1: Stop criterion is iteration for K+? times
% %         for k=1:K+0;                                       
% %             for col=1:N;                                  
% %                 product(col)=abs(Phi(:,col)'*r_n);        %  Step 1): calculate the correlation between thecol-th column and the residual signal 
% %             end
% %             [val,pos]=max(product);                       %  Step 2): find the value and position of the most matching column
% %             Aug_t=[Aug_t,Phi(:,pos)]; pos_array(k)=pos;   %  Step 3): update the argment matrix and record the position,
% %             Phi(:,pos)=zeros(M,1);                        %           remove the column just selected
% %             aug_x=(Aug_t'*Aug_t)^(-1)*Aug_t'*y;           %  Step 4): solve the LS problem to obatain a new signal estimate
% %             r_n=y-Aug_t*aug_x;                            %  Step 5): Calculate the new residual
% %         end                                               %  Step 6): return to Step 2) if k < K 
% %         hat_x(pos_array)= aug_x;                          %  Step 7): get the recovered signal
% %         Ch_est_CS_TS = hat_x;
%         
% %         % Implementation 2:  Adaptive stop criterion: iteration for K times or the residual error is small enough    
% %         r_n = y; index_set = []; Aug_t = []; k = 1;                       % Step 1): initialize the residual, the index set, the chosen atom, and the iteration counter
% %         err = 1e-5;          
% %         while (k < K)  && (norm(r_n)>err)                 
% %             [val,pos]=max(Phi'* r_n);                                     %  Step 2): find the value and position of the most matching column
% %             index_set = [index_set; pos]; Aug_t=[Aug_t,Phi(:,pos)];       %  Step 3): argment the index set and the matrix of chosen atoms
% %             Phi(:,pos)=zeros(M,1); aug_x = (Aug_t'*Aug_t)^(-1)*Aug_t'*y;  %  Step 4): solve a least square problem to obtain a new signal estimate
% %             r_n=y-Aug_t*aug_x;                                            %  Step 5): calculate the new apporoxiation of the data and the new residual
% %             k = k + 1;                                                    %  Step 6): inrement k, and return to Step 2) if k< 2*K (in the paper, k=K, normally, m>=K is required)
% %         end
% %         hat_x(index_set)=aug_x;                                           %  Step 7): get the recovered signal
% %         Ch_est_CS_TS = hat_x;
% 
%         % CoSaMP for signal recovery
%         max_iter = K;
%         [Ch_est_CS_TS used_iter]= CoSaMP(y,Phi,K,max_iter);
%         
%         % Modified-CoSaMP wiht partial priori information
%         max_iter = K;
%         Th = 0.1; % The threshold should be optimized
%         T_0 = find(abs(Ch_est_direct_correlation)> Th);  % Most significant taps
%         Known_support = size(T_0,1); % Known_support 
%         [Ch_est_CS_TS_1 used_iter_1]= M_CoSaMP(y,Phi,K,T_0,max_iter-Known_support);
%         
%         % Step 3: ML path coefficient estimation
%         support = find(Ch_est_CS_TS_1 > 0);
% %         Num = size (support,1) + 10;
%         Num = M;
%         Row_sel = randperm(M);
% %         rows = sort(Row_sel(M-Num+1:M),'ascend');
%         rows = M -Num + 1: M;
%         Phi_for_ML = Phi(rows,support);
%         Path_coefficient = pinv(Phi_for_ML) * y(rows);
%         Ch_est_CS_TS_ML = zeros(PN_len,1);
%         Ch_est_CS_TS_ML(support) = pinv(Phi_for_ML) * y(rows);
% 
% %         % Rugular_OMP for signal recovery
% %         [Ch_est_CS_TS, numIts] = Rugular_OMP(K, Phi, y);
% 
% %         % SAMP (Sparse Adaptive MP) for signal recovery: unique feature is that K is not required
% %         step_size = 1; % should be smaller than K
% %         sigma = 10^(-SNR(i_snr)/10); % noise power
% %         [Ch_est_CS_TS iter_num] = SAMP(y, Phi, step_size, sigma);
%         
% %         % Subspace_Pursuit (+ lp-Reweighted Least Squares Problem Solver) for signal recovery
% %         [Ch_est_CS_TS,That]=Subspace_Pursuit(K, Phi, y);
% % %         p=0.9; [xfinal,That]=Subspace_Pursuit(K, Phi, y); [Ch_est_CS_TS] = lp_re1(Phi,y,p,xfinal); % Subspace_Pursuit + lp-Reweighted
% 
% %         % AIHT for signal recovery  [s, err_mse, iter_time]=AIHT(x,P,N,K)
% %         thresh=1e-4;
% % %         Ch_est_CS_TS = AIHT(y, Phi, N, K,thresh);
% % %         Ch_est_CS_TS = IHT(y,Phi,K);
% %         Ch_est_CS_TS = hard_l0_Mterm_mod(y,Phi,N,K,thresh);
%         
% %             figure;
% %             stem(abs(Ch_est_direct_correlation),'bo');
% %             hold on;
% %             stem(abs(Ch_est_unnoised_preamble),'m+');
% %             stem(abs(Ch_est_CS_TS),'k^');
% %             stem(abs(Ch_est_CS_TS_1),'ms');
% %             stem(abs(Ch_est_CS_TS_ML),'cp');
% %             stem(abs(Static_multipath_channel),'r.');
% %             legend('Correlation-Based','CoSaMP','Modified-CoSaMP','ML','Actual Channel');
% %             axis([0  max(Tao_l)+5 0 1]);
% 
% %        % When 2L-1<=M, using last >=L samples without IBI to recovery channel via LS
% %        if  (max(Tao_l))< PN_len/2
% % %             rows_for_LS = max(Tao_l) + 1: PN_len;
% %             rows_for_LS = ceil(PN_len/2): PN_len;
% %             Phi_for_LS = PN256_matrix(rows_for_LS,:);
% % %             Phi_for_LS = PN256_matrix(rows_for_LS,1:ceil(PN_len/2)-1);
% %             cond(Phi_for_LS)     
% %             Second_part_TS = (Frame_current(rows_for_LS) + Frame_current(rows_for_LS + Frame_len))/2;
% %             Ch_est = pinv(Phi_for_LS) * Second_part_TS;
% %             Ch_est_Direct_LS = zeros(PN_len,1);
% %             Ch_est_Direct_LS(1:ceil(PN_len/2)-1) = Ch_est(1:ceil(PN_len/2)-1);
% % %             Ch_est_Direct_LS = zeros(ceil(PN_len/2),1);
% % %             Ch_est_CS_TS_ML(support) = pinv(Phi_for_LS) * y(rows_for_LS);
% %        end
%        
%        % Final CIR selection
% %        CIR_estimation_final = Ch_est_unnoised_preamble;
%        CIR_estimation_final = Ch_est_CS_TS;    % based on CS 
% %        CIR_estimation_final = Ch_est_CS_TS_1;  % based on Modified CS
% %        CIR_estimation_final = Ch_est_CS_TS_ML;   % Modified CS + ML
% %        CIR_estimation_final = Ch_est_Direct_LS;  % Direct_LS when 2L-1<=M
% %        CIR_estimation_final = Static_multipath_channel; % ideal channel
% 
% %        figure
% %        stem(abs(Ch_est_CS_TS));
% %        hold on;
% %        stem(abs(Static_multipath_channel),'r.');
% %        stem(abs(CIR_estimation_final),'m^');
% %        legend('CoSaMP','Ideal Channel', 'Modified CoSaMP')
% 
%         if i_subframe == 1
%             Ch_est_unnoised_1_current = Ch_est_unnoised_preamble;
%         else
%             Ch_est_unnoised_1_current = CIR_estimation_final_prv;
%         end
%         Ch_est_unnoised_1_next = CIR_estimation_final;
%         
%         Frame_data_time = Frame_current(PN_len + 1:Frame_len);
%         if PN_ZP == 1 % TS is the guard interval
%             PN_current_tail = conv(ZP256,Ch_est_unnoised_1_current);
%             PN_next_tail = conv(ZP256,Ch_est_unnoised_1_next);
%     % %         Frame_data_time(1:PN_len-1) = Frame_data_time(1:PN_len-1) - PN_current_tail(PN_len + 1:end) + (Frame_current(Frame_len + 1:Frame_len+PN_len-1)-PN_next_tail(1:PN_len-1));
%             Frame_data_time(1:Chan_Len) = Frame_data_time(1:Chan_Len) - PN_current_tail(PN_len + 1:PN_len+Chan_Len) + (Frame_current(Frame_len + 1:Frame_len+Chan_Len)-PN_next_tail(1:Chan_Len));
%         elseif PN_ZP == 0  % Zero padding
%             Frame_data_time(1:Chan_Len) = Frame_data_time(1:Chan_Len) + Frame_current(Frame_len + 1:Frame_len + Chan_Len); % OLA
%         end
% %        % CPR
% %        Frame_data_time(1:Chan_Len) = Frame_data_time(1:Chan_Len) + Frame_current(Frame_len + 1:Frame_len+Chan_Len) - Frame_current(PN_len + 1: PN_len + Chan_Len);
%                 
%        Frame_data_freq = fft(Frame_data_time, Idft_len) ./ sqrt(Idft_len);
% 
% %        % CE based on grouped pilots, direct method (maybe ill-posed)
% %        Received_Central_pilot = Frame_data_freq(Pilot_group_center);   % Subtract the frequency-domain pilot  (before equalization)
% %        CIR_alpha_direct = Moore_Penrose_matrix * Received_Central_pilot./Pilot_amptitude;
% %        CIR_direct_estimate = zeros(PN_len,1);
% %        CIR_direct_estimate(Tao_l + 1) = CIR_alpha_direct;
% 
% % %             Direct iverse problem using the last part of the TS
% %             Ch_est_direct = Moore_Penrose_PN256_matrix * Frame_current(PN_len - 90 +1:PN_len);                      % Phi will be change in OMP
% %             CIR_direct_estimate = zeros(PN_len,1);
% %             CIR_direct_estimate(Tao_l + 1) = Ch_est_direct;
%             
% %             figure;
% %             stem(abs(CIR_direct_estimate),'r.');
% %             hold on;
% %             stem(abs(Static_multipath_channel),'ks');
% %             legend('Channel Estimation','Ideal Channel');
% 
%        % MSE calculation
%         error = sum(abs(CIR_estimation_final -  Static_multipath_channel).^2)/sum(abs(Static_multipath_channel)).^2;
%         Th_mse = 10^(-2);
%         if error < Th_mse
%             Probability_correct_recovery((i_frame-1)*Sub_frame + i_subframe,m_index) = 1;
%         end
%         MSE_channel((i_frame-1)*Sub_frame + i_subframe,m_index) = error;
% 
%      
%        % 2009.5.27: average H when equalization
%        Ch_est_unnoised_fre = fft(CIR_estimation_final,Idft_len);  % for normalized FFT/IFFT, Y = sqrt(N) * H * X = [sqrt(N) * H] * X, where sqrt(N) * H = fft(h)
%        csi = abs(Ch_est_unnoised_fre);
%        snr = 10.^(SNR(i_snr)/10);  
%        Frame_equalized = Frame_data_freq.*conj(Ch_est_unnoised_fre) ./ (csi).^2;  % ZF Equalizer
% %        Frame_equalized = Frame_data_freq.*conj(Ch_est_unnoised_fre) ./ ((csi).^2 + 1/snr);   % MMSE Equalizer
%        Frame_equalized_interger = Frame_equalized(Data_index) * sqrt(Mod_power);   % restore the power of modulated data
%                    
% %         figure(1);
% %         stem(abs(Static_multipath_channel),'r.');
% %         hold on;
% %         stem(abs(Ch_est_unnoised_preamble),'gs');
% %         stem(abs(CIR_direct_estimate),'m>');
% %         stem(abs(h_lsqi_final),'bp');
% %         stem(abs(h_discrep_final),'k+');
% %         axis([1 PN_len 0 1]);
% %         legend('Ideal Channel', 'Preamble', 'Direct estimate','LSQI','Discrep')
% % %         pause
%         
% %         Mod_symbol_est = qamdemod(Frame_equalized_interger,2^Mod_order(ii),0,'gray');
% %         Payload_bit_est = dec2bin(Mod_symbol_est, Mod_order(ii));
% %         Payload_bit_est = Payload_bit_est(:);
% %         ber = length(find(Payload_bit_est~= Send_data_buf((i_subframe-1)*(Idft_len-Pilot_number)*Mod_order(ii)+1:i_subframe*(Idft_len-Pilot_number)*Mod_order(ii),i_frame)))/length(Payload_bit_est);
% %         BER_store(i_subframe,i_frame) = ber;
% %         Receive_data_buf((i_subframe-1)*(Idft_len-Pilot_number)*Mod_order(ii)+1:i_subframe*(Idft_len-Pilot_number)*Mod_order(ii),i_frame) = Payload_bit_est; % added for longtime SER
%         
%         CIR_estimation_final_prv = CIR_estimation_final; % CIR of the previous subframe
    end
    
   end
%    BER(i_snr) = length(find(Send_data_buf(:)~= Receive_data_buf(:)))/length(Send_data_buf(:)) % get the BER
  end
end



%%  Plot
% figure(1);
% semilogy(SNR,BER,'ro-','LineWidth',1.5);
% axis([5 30 1e-3 1e0]);
% xlabel('SNR (dB)');
% ylabel('BER');
% grid on;
% hold on;
% legend('Compressive Sensing');

% figure(2);
% plot(num_obs,sum(Probability_correct_recovery)/(Sim_frame_num*Sub_frame),'bo-','LineWidth',1.5);
% hold on;
% plot(num_obs,sum(Probability_correct_recovery_new)/(Sim_frame_num*Sub_frame),'rs-','LineWidth',1.5);
% plot(num_obs,sum(Probability_correct_recovery_new_ML)/(Sim_frame_num*Sub_frame),'mp-','LineWidth',1.5);
% axis([min(num_obs) max(num_obs) 0 1]);
% grid on;
% xlabel('Number of Observations');
% ylabel('Recovery Probability');
% legend('SOMP','PA-SOMP','PA-SOMP-ML');

figure(3);
% semilogy(num_obs,mean(MSE_channel),'bo-','LineWidth',1.5);
hold on;
% semilogy(num_obs,mean(MSE_channel_new),'rs-','LineWidth',1.5);
semilogy(num_obs,mean(MSE_channel_new_ML),'mp-','LineWidth',1.5);
axis([min(num_obs) max(num_obs) 1e-4 1e0]);
xlabel('M');
ylabel('MSE');
grid on;
hold on;
% legend('SOMP','PA-SOMP','PA-SOMP-ML');