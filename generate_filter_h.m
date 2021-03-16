function h_r_t = generate_filter_h(W,N,R)
% (Ϊʲôʱ���-4��4�أ�����8����Ԫ�ķ�Χ�������պͷ����˲����������ֵĳ����Ӧ
%  ʱ��ϵͳ����
%  W�źŴ�������ȷ���������Ų��ε�ʱ��
%  N�²���ϵ����һ�����ŵĲ��εĲ�������
%% Parameters of the filter 
%  We consider RRC pulses as transmitt and receive shaping pulses 
R=0.22; % roll off factor
N=32; % downsampling factor
W=500e6;
tt=linspace(-4,4,8*N+1); 
%% Generation of transmit and receive shaping pulse
% We use RRC shaping filters

% RRC transmitter shaping filter
len=length(tt);
rrc_t=zeros(len,1);
for i=1:len
   t=tt(i);      
   if(t==0)
       rrc_t(i)= ( 1-R+4*R/pi ) ;
   elseif(abs(abs(t)-1/4/R)<1e-3)
       rrc_t(i)= (  cos(pi*t*(1-R))*pi*(1-R) + 4*R*cos(pi*t*(1+R)) - 4*R*t*sin(pi*t*(1+R))*(pi*(1+R))  )/(pi)/(1-3*(4*R*t)^2   ) ;   
   else
       rrc_t(i)= ( sin(pi*t*(1-R))+4*R*t*cos(pi*t*(1+R)) ) / (  pi*t*(1- (4*R*t)^2)  );
   end
end

% Normalization of filter as unitary energy filter
E_rrc_t=rrc_t'*rrc_t;
rrc_t=rrc_t/sqrt(E_rrc_t);

% RRC receive shaping filter
rrc_r=rrc_t; 

% Convolution between receiver and transmitter shaping filters and
% normalization
h_r_t=conv(rrc_t,rrc_r);
E_h=h_r_t'*h_r_t;
h_r_t=h_r_t/sqrt(E_h);
end

