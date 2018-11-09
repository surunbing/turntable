function [Design] = Designplastic(data, advance, data_cur)
%   �������� ʹ�ó�ǰ�ٺ󻷽ںͷǵ��͵Ķ��׻���   ��״�ж��������ԸĽ�
%   Լ����С��λ�� ����cost function
global PHI_MIN BANDWIDTH

%% ����Ƿ���Խ������
%% �϶����Խ������
num = find(data.fre == advance.fre);
frequence = data.fre(1 : num);
phi_min = min(data_cur.phi);
num = find(data_cur.phi == phi_min);
num = num(1);
fre_phimin = frequence(num);
if fre_phimin < -90
    %% ��Ƶ�ʺ������λ�����Ƶ��
    num = find(data_cur.fre > 2 * pi);
    num = num(1);
    phi_max = data_cur.phi(num : end);
    num = find(data_cur.phi == phi_max);
    num = num(1);
    fre_phimax = data_cur.fre(num);
    if(fre_phimax < BANDWIDTH)
        bPhiFillSlow = 1;
    else 
        bPhiFillSlow = 0;
    end
    %% �»������ǣ����Ҽ�¼
    if bPhiFillSlow == 1
        num = find(data_cur.fre > BANDWIDTH);
        num = num(1);
        phi_max = data_cur.phi(num : end);
        num = find(data_cur.phi == phi_max);
        num = num(1);
        fre_phimax = data_cur.fre(num);
    end
    phi_reg = fre_phimin - PHI_MIN;
    % ��������;
    bNormal = 1;
else
    % �������������
    bNormal = 0;
end




% ��ñջ�
end

