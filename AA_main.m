clear all
close all
clc
tic
%% 滤波器参数
N = 32; % downsampling factor
W = 500e6;
Rool = 0.22;
T_symbol=(1+Rool)/W; % symbol time
Ts=T_symbol/N; % sampling time for the filter
Tc=T_symbol; % sampling time for the output of the receive filter 
h_r_t = generate_filter_h(W,N,Rool);
%%  Parameter for transmitter and receiver planar arrays 
Yt=5; % number of transmit antennas on the y-axis of planar array
Zt=4; % number of transmit antennas on the z-axis of planar array
Yr=5; % number of receiver antennas on the y-axis of planar array
Zr=2; % number of receiver antennas on the z-axis of planar array
Ysurface = 6;
Zsurface = 6;
% 
f=73e09; % carrier frequency
K = 5; % number of serving users 同时服务用户数
% Positions of transmitter and receiver in 3-D plane
TX_pos=[0 0 7];
Surface_pos = [30 30 5];
RX_pos_vector =generate_Kusers_pos(K); 
% scenario: variable that contains information about the use-case scenario,
% it assumes the values:
% - scenario==1  ==> 'Open square'
% - scenario==2  ==> 'Street Canyon'
% - scenario==3  ==> 'Indoor Office'
% - scenario==4  ==> 'Shopping mall'
scenario=2; % UMi-Street Canyon
% 非移动状态
v_RX= 0; % speed of the receiver in km/h
v_TX= 0; % speed of the transmitter in km/h
%% Frequency selective Channel Matrix LTI case 
ifnorm = 1;
Gtemporal =Generate_Channel_frequency_selective_LTI(f,TX_pos,Surface_pos,scenario,Yt,Zt,Ysurface,Zsurface,h_r_t,Ts,Tc);
G = Covert_Htime2Hfrq(Gtemporal,f,Tc,ifnorm);
for kth = 1:K
Htemporal{kth} =Generate_Channel_frequency_selective_LTI(f,TX_pos,RX_pos_vector(kth,:),scenario,Yt,Zt,Yr,Zr,h_r_t,Ts,Tc);
H{kth} = Covert_Htime2Hfrq(Htemporal{kth},f,Tc,ifnorm)
Rtemporal{kth} = Generate_Channel_frequency_selective_LTI(f,Surface_pos,RX_pos_vector(kth,:),scenario,Ysurface,Zsurface,Yr,Zr,h_r_t,Ts,Tc);
R{kth} = Covert_Htime2Hfrq(Rtemporal{kth},f,Tc,ifnorm)
v_vec{kth} = randn*exp(-sqrt(-1)*2*pi*rand(Yr*Zr,1));
end
N_RF = 6;
D = ones(K,1);
W = randn(N_RF,K)+sqrt(-1)*randn(N_RF,K);
F = exp(-sqrt(-1)*2*pi*rand(Yt*Zt,N_RF));
theta_vector = exp(-sqrt(-1)*2*pi*rand(Ysurface*Zsurface,1));
Theta = diag(theta_vector);
X = F*W*D;
D_estimate = zeros(K,1);
Z = randn(Yr*Zr,1)+exp(sqrt(-1))*randn(Yr*Zr,1);
for kth = 1:K
Y{kth} = H{kth}*X+R{kth}*Theta*G*X+Z;
D_estimate(kth) = v_vec{kth}'*Y{kth};
end
SSE = sum(abs(D-D_estimate).^2);%再对D和Z求期望就能得到SMSE
