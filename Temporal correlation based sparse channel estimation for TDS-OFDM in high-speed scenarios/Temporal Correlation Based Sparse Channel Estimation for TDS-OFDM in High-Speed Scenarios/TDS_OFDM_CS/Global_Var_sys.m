
SC =0;                % 1: Single-carrier; 0: Multi-carrier(OFDM)
Virtual_frame = 0;    % 1: Virtual_frame; 0: tail_subtract
% Sim_frame_num = 10; % Number of frames for simulation
% MultipathType = 2;    % 0: AWGN; 1: Static Multi-path; 2: Rayleigh
Frame_number_total = 225;
% PN_len = 420;      % for CAZAC: *2; for PN-MC: *1;
PN_len = 255; 
PN_cyclic_Len = 0; % Dual PN
PN_total_Len = PN_len + PN_cyclic_Len;
PN_power = 3;        % the power of PN is 3dB higher than DFT data
Idft_len = 3780;   
Data_len = Idft_len-36;  % exclude TPS
Frame_len = PN_total_Len + Idft_len;  
Ca_prc_start_lshift = 0;  % 同步模块给出的信号，告诉主径在哪一类的
mmseAlpha = 0.9997; 
% mmseAlpha = 0.99; 

Noise_Threshold = 0.05;

Srrc_oversample = 1; 


Ca_prc_start_index = 0;
Pre_len_shift = 0; %(-3) * Srrc_oversample;
Trans_scale = 24000/10^(PN_power/20)/2/512;
Receive_scale = 0.5;

Symbol_rate = 7.56e6;
Sampling_rate = Symbol_rate * Srrc_oversample;

Fft_delay = 3813;
State_channel = [];
mldpc=8;
% Static_multipath_channel = zeros(255*Srrc_oversample,1);

% Srrc_factor = 0.05; Srrc_trans_delay = 128; 
% Srrc_trans_filter = rcosine(1, Srrc_oversample, 'fir/sqrt', Srrc_factor, Srrc_trans_delay);  % 得到的数据为时域中的，sa2波形
% Srrc_receive_delay = 64; 
% Srrc_receive_filter = rcosine(1, Srrc_oversample, 'fir/sqrt', Srrc_factor, Srrc_receive_delay);
% 
% Ca_prc_start_index = (Srrc_trans_delay + Srrc_receive_delay - Ca_prc_start_lshift) * Srrc_oversample + 1;  % 接收端延迟
% Pre_len_shift = 0; 

Symbol_rate = 7.56e6;
% Symbol_rate = 10e6;
Sampling_rate = Symbol_rate * Srrc_oversample;

%% Channel
Static_multipath_channel = zeros(PN_len*Srrc_oversample,1);  
Static_multipath_channel_1 = zeros(PN_len*Srrc_oversample,1);  
Static_multipath_channel_2 = zeros(PN_len*Srrc_oversample,1);  

% AWGN
% Static_path = round([0] * 1e-6 * Symbol_rate)*Srrc_oversample + 1;   
% Static_multipath_channel(round([0] * 1e-6 * Symbol_rate)*Srrc_oversample + 1) = 10.^([0]/20);

% Static_path = round([0.0 0.31  0.71  1.09  1.73  2.51] * 1e-6 * Symbol_rate)*Srrc_oversample + 1;     % DVB-T channel (ITU-R Vehicular A)
% Static_multipath_channel(round([0.0 0.31  0.71  1.09  1.73  2.51] * 1e-6 * Symbol_rate)*Srrc_oversample + 1) = 10.^([0.0 -1.0 -9.0 -10.0 -15.0 -20]/20);
%  
% Static_path = round([0.0 0.30  8.90 12.90 17.10 20.00] * 1e-6 * Symbol_rate)*Srrc_oversample + 1;     % ITU-R Vehicular B
% Static_multipath_channel(round([0.0 0.30  8.90 12.90 17.10 20.00] * 1e-6 * Symbol_rate)*Srrc_oversample + 1) = 10.^([-2.5 0.0 -12.8 -10.0 -25.2 -16.0]/20);

% Static_path = round([0.0 0.30  8.90 12.90 17.10 20.00 33] * 1e-6 * Symbol_rate)*Srrc_oversample + 1;     % ITU-R Vehicular B with one artifical path
% Static_multipath_channel(round([0.0 0.30  8.90 12.90 17.10 20.00 33] * 1e-6 * Symbol_rate)*Srrc_oversample + 1) = 10.^([-2.5 0.0 -12.8 -10.0 -25.2 -16.0 -16.0]/20);


% Static_path = round([0.0 0.2 0.5] * 1e-6 * Symbol_rate)*Srrc_oversample + 1;     % self-defined for test 1
% Static_multipath_channel(round([0.0 0.2 0.5] * 1e-6 * Symbol_rate)*Srrc_oversample + 1) = 10.^([0 0 0]/20);

% Static_path = round([0.0 0.4 0.8] * 1e-6 * Symbol_rate)*Srrc_oversample + 1;     % self-defined for test 2
% Static_multipath_channel(round([0.0 0.4 0.8] * 1e-6 * Symbol_rate)*Srrc_oversample + 1) = 10.^([0 0 0]/20);
% 

% Static_path = round([0.0 0.15 2.22 3.05 5.86 5.93] * 1e-6 * Symbol_rate)*Srrc_oversample + 1;     % Brazil A channel
% Static_multipath_channel(round([0.0 0.15 2.22 3.05 5.86 5.93] * 1e-6 * Symbol_rate)*Srrc_oversample + 1) = 10.^([0.0 -13.8 -16.2 -14.9 -13.6 -16.4]/20);

% Static_path = round([0.0 0.30 3.50 4.40 9.50 12.7] * 1e-6 * Symbol_rate)*Srrc_oversample + 1;   % Brazil B channel
% Static_multipath_channel(round([0.0 0.30 3.50 4.40 9.50 12.7] * 1e-6 * Symbol_rate)*Srrc_oversample + 1) = 10.^([0.0 -12.0 -4.0 -7.0 -15.0 -22.0]/20); 

Static_path = round([0.15 0.63  2.22  3.05 5.86 5.93] * 1e-6 * Symbol_rate)*Srrc_oversample + 1;     % Brazil D channel
Static_multipath_channel(round([0.15 0.63  2.22  3.05 5.86 5.93] * 1e-6 * Symbol_rate)*Srrc_oversample + 1) = 10.^([-0.1  -3.8  -2.6  -1.3  0.0  -2.8]/20);

% Static_path = round([0    -1.8+18  0+18      0.15+18  1.8+18    5.7+18] * 1e-6 * Symbol_rate)*Srrc_oversample + 1;     % GD6 channel
% Static_multipath_channel(round([0    -1.8+18  0+18      0.15+18  1.8+18    5.7+18] * 1e-6 * Symbol_rate)*Srrc_oversample + 1) = 10.^([-10  -20    0  -20  -10  -14]/20);
 

% Static_path = round([0    1.8      0.15+1.8  1.8+1.8  5.7+1.8   18+1.8] * 1e-6 * Symbol_rate) *Srrc_oversample+ 1;   % Guangdian 7 channel
% Static_multipath_channel(round([0    1.8      0.15+1.8  1.8+1.8  5.7+1.8   18+1.8] * 1e-6 * Symbol_rate)*Srrc_oversample + 1) = 10.^([-20    0  -20  -10    0  -18]/20);

% Static_path = round([0    1.8      0.15+1.8  1.8+1.8  5.7+1.8   18+1.8] * 1e-6 * Symbol_rate) *Srrc_oversample+ 1;   % Modified Guangdian 7 channel
% Static_multipath_channel(round([0    1.8      0.15+1.8  1.8+1.8  5.7+1.8   18+1.8] * 1e-6 * Symbol_rate)*Srrc_oversample + 1) = 10.^([-20    0  -20  -10     -18 0 ]/20);

% Static_path = round([0.0 1.80 1.95 3.60 7.50 31.8] * 1e-6 * Symbol_rate) *Srrc_oversample+ 1;   % Guangdian 8 channel
% Static_multipath_channel(round([0.0 1.80 1.95 3.60 7.50 31.8] * 1e-6 * Symbol_rate)*Srrc_oversample + 1) = 10.^([-18.0 0.0 -20.0  -20.0  -10.0  0.0]/20);

Static_multipath_channel = Static_multipath_channel./norm(Static_multipath_channel);
% Static_multipath_channel_1 = Static_multipath_channel_1./norm(Static_multipath_channel_1);
% Static_multipath_channel_2 = Static_multipath_channel_2./norm(Static_multipath_channel_2);
Static_multipath_channel_2 = Static_multipath_channel_1; % just test

% Static_multipath_channel_1(Static_path) = [-0.1330 + 0.2197i   0.0820 - 0.1236i  -0.0515 + 0.0906i -0.1123 + 0.1824i   0.0511 + 0.0194i  -0.0323 + 0.1669i];
% Static_multipath_channel_2(Static_path) = [-0.1179 - 0.3576i  -0.0109 + 0.0391i  -0.0755 - 0.0206i  0.2009 + 0.2645i  -0.0550 + 0.2810i   0.1013 + 0.1962i];

% load path_power.mat;   % PLC channel
% load path_delay_us.mat;
% path_power = aa;
% path_delay_us = delay_us;
% Static_path = round([path_delay_us] * 1e-6 * Symbol_rate*Srrc_oversample) + 1;    
% Static_multipath_channel(round([path_delay_us] * 1e-6 * Symbol_rate*Srrc_oversample) + 1) = path_power;

Ts = 1/Symbol_rate/Srrc_oversample;

DopplerShift = 85;
% DopplerShift = 20*(256*2 + 3780)/Frame_len; % for CAZAC;

% Rayleigh_path = round([0.0 0.31  0.71  1.09  1.73  2.51] * 1e-6 * Symbol_rate)*Srrc_oversample+1;   % DVB-T  (ITU-R Vehicular A)
% Delay_vec = round([0.0 0.31  0.71  1.09  1.73  2.51] * 1e-6 * Symbol_rate) / Symbol_rate; 
% Power_vec = [0.0 -1.0 -9.0 -10.0 -15.0 -20];  
% Rayleigh_path = round([0.0 0.20 0.80 1.20 2.30 3.70] * 1e-6 *
% Symbol_rate)*Srrc_oversample+1;   % ITU-R Pedestrian B
% Delay_vec = round([0.0 0.20 0.80 1.20 2.30 3.70] * 1e-6 * Symbol_rate) / Symbol_rate; 
% Power_vec = [0.0 -0.9 -4.9 -8.0 -7.8 -23.9]; 

% Rayleigh_path = round([0.0 0.31  0.71  1.09  1.73  2.51] * 1e-6 * Symbol_rate)*Srrc_oversample+1;     % ITU-R Vehicular A
% Delay_vec = [0.0 0.31  0.71  1.09  1.73  2.51] * 1e-6; 
% Power_vec = [0.0 -1.0 -9.0 -10.0 -15.0 -20];  

Rayleigh_path = round([0.0 0.30  8.90 12.90 17.10 20.00] * 1e-6 * Symbol_rate)*Srrc_oversample+1;   % ITU-R Vehicular B
% Delay_vec = [0.0 0.30  8.90 12.90 17.10 20.00] * 1e-6; 
Delay_vec = (Rayleigh_path-1)/Symbol_rate;  % to set the sample-spaced Rayleigh channel
Power_vec = [-2.5 0.0 -12.8 -10.0 -25.2 -16.0];  

% Rayleigh_path = round([0.0 0.15 2.22 3.05 5.86 5.93] * 1e-6 * Symbol_rate)*Srrc_oversample+1;   % Brazil A
% Delay_vec = [0.0 0.15 2.22 3.05 5.86 5.93] * 1e-6 ; 
% Delay_vec = (Rayleigh_path-1)/Symbol_rate;  % to set the sample-spaced Rayleigh channel
% Power_vec = [0.0 -13.8 -16.2 -14.9 -13.6 -16.4];  

% Rayleigh_path = round([0.0 0.30 3.50 4.40 9.50 12.7] * 1e-6 * Symbol_rate)*Srrc_oversample+1;  % Brazil B
% % Delay_vec = [0.0 0.30 3.50 4.40 9.50 12.7] * 1e-6 ; 
% Delay_vec = (Rayleigh_path-1)/Symbol_rate;  % to set the sample-spaced Rayleigh channel
% Power_vec = [0.0 -12.0 -4.0 -7.0 -15.0 -22.0];  
% % % 

% Rayleigh_path = round([0.15 0.63  2.22  3.05 5.86 5.93] * 1e-6 * Symbol_rate)*Srrc_oversample+1;  % Brazil D
% % Delay_vec = [0.15 0.63  2.22  3.05 5.86 5.93] * 1e-6 ; 
% Delay_vec = (Rayleigh_path-1)/Symbol_rate;  % to set the sample-spaced Rayleigh channel
% Power_vec = [-0.1  -3.8  -2.6  -1.3  0.0  -2.8];  

% Rayleigh_path = round([0.0 1.0 2.0] * 1e-6 * Symbol_rate)*Srrc_oversample+1;  % Brazil E
% Delay_vec = [0.0 1.0 2.0] * 1e-6 ; 
% Power_vec = [0.0 0.0 0.0];   

% Rayleigh_path = round([0    -1.8+18  0+18      0.15+18  1.8+18    5.7+18] * 1e-6 * Symbol_rate)*Srrc_oversample+1;   % Guangdian 6 channel
% Delay_vec = [0    -1.8+18  0+18      0.15+18  1.8+18    5.7+18] * 1e-6 ; 
% Power_vec = [-10  -20    0  -20  -10  -14]; 

% Rayleigh_path = round([0.0 1.80 1.95 3.60 7.50 31.8] * 1e-6 * Symbol_rate)*Srrc_oversample+1;  % Guangdian 8 channel
% % Delay_vec = [0.0 1.80 1.95 3.60 7.50 31.8] * 1e-6; 
% Delay_vec = (Rayleigh_path-1)/Symbol_rate;  % to set the sample-spaced Rayleigh channel
% Power_vec = [-18.0 0.0 -20.0  -20.0  -10.0  0.0]; 

RayleighChannel = rayleighchan(Ts,DopplerShift,Delay_vec,Power_vec);
reset(RayleighChannel,[10; 20])    % initialize the seed 
RayleighChannel.ResetBeforeFiltering = 0;

RayleighChannel_1 = rayleighchan(Ts,DopplerShift,Delay_vec,Power_vec);
reset(RayleighChannel_1,[10; 20]);     % initialize the seed 
RayleighChannel_1.ResetBeforeFiltering = 0;

RayleighChannel_2 = rayleighchan(Ts,DopplerShift,Delay_vec,Power_vec);
reset(RayleighChannel_2,[15; 25]);     % initialize the seed
RayleighChannel_2.ResetBeforeFiltering = 0;


