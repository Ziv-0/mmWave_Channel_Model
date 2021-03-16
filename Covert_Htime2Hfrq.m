function Hout =Covert_Htime2Hfrq(Ht,f0,tsample)
% ��խ��������
% ��ʱ����ʽ���ŵ�״̬ת��ΪĳƵ��f0�����ŵ�����
% ������ʱ����ϵͳ�����ʿ�֪����Ϊʱ���ź���f0���ĸ���Ҷ�任
% Htʱ����ʽ���ŵ������Ӧ��Nr*Nt*Nsample)
% tsample Ht��ʱ��������
% f0 Ŀ��Ƶ��
% Hout ������ŵ�����(Nr*Nt)
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

