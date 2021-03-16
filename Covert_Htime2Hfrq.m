function Hout =Covert_Htime2Hfrq(Ht,f0,tsample)
% 在窄带条件下
% 将时域形式的信道状态转化为某频点f0处的信道矩阵
% 由线性时不变系统的性质可知，即为时域信号在f0处的傅里叶变换
% Ht时域形式的信道冲击响应（Nr*Nt*Nsample)
% tsample Ht的时间采样间隔
% f0 目标频点
% Hout 输出的信道矩阵(Nr*Nt)
[Nr,Nt,Nsample] = size(Ht);
Hout = zeros(Nr,Nt);
for ith_Nr = 1:Nr
    for jth_Nt = 1:Nt
        h = 0;
        for kth_Nsample = 1:Nsample
            h = h+tsample*Ht(ith_Nr,jth_Nt,kth_Nsample)*exp(-j*2*pi*f0*(kth_Nsample-1)*tsample);
        end
        Hout(ith_Nr,jth_Nt) = h;
    end
end
end

