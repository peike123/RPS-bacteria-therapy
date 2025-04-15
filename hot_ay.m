
function r
global B1 B2 B3 E T total_time
global alpha m  A_c gamma alpha_A alpha_AA gamma_A    s n mu eta1 g1 r k e1     eta2 g2 e2 beta g3 e3      v b omega gamma_C t0  tau

total_time = 600;  % ��ʱ��
E = 1e5; T = 1e7; B1 = 0.25; B2 = 0.3; B3 = 0.1; %��ֵ��

m = 2; A_c = 1; alpha_A = 0.4; alpha_AA = 8; gamma_A = 1;  
s = 1.3e4; n = 3.422e-10; mu = 0.0412; eta1 = 0.1245; g1 = 2.019e7;  r = 0.5; k = 5e8; e1 = 1.101e-7;
eta2 = 0.2; g2 = 2.5e7; e2 = 0.4; beta = 5; g3 = 2e7; e3 = 0.8; 
v = 0.05; omega = 5; gamma_C = 1; t0 = 0;  tau = 0.05; b = 0.1;
%alpha = 1; gamma = 4; 

%alpha gamma�䶯����ͼ 
alpha_values = 0:1:30;
gamma_values = 1:1:10;

% �洢����ֵ�����飨������������
final_tumor9 = zeros(length(alpha_values), length(gamma_values));

% %�洢ЭͬЧӦ����E������
% efficacy9 = zeros(length(b_values), length(v_values));

    for idx = 1:length(alpha_values)
        for idy = 1:length(gamma_values)
            alpha = alpha_values(idx);
            gamma = gamma_values(idy);
            
%             B1 = 0.2225 + (0.275 - 0.225) * rand;
%             B2 = 0.27 + (0.33 - 0.27) * rand;
%             B3 = 0.09 + (0.11 - 0.09) * rand; %0.25 0.3 0.1
%             
%             E = randi([1e5 1.5e5]);  %����[]֮�������
%             T = randi([1e7 1.5e7]);
% %             E = randi([0.95e5 1.05e5]);  %����[]֮������� 1e5 1e6
% %             T = randi([0.95e7 1.05e7]);

            [T9,Y9] = jiujie();
        
            % �������չ�ģ
            final_tumor9(idx,idy) = Y9(end,9);
        
        end
    end
    
%     %����ЭͬЧӦ����E
%     for idx = 1:length(b_values)
%         for idy = 1:length(v_values)
%             efficacy9(idx,idy) = (final_tumor9(1,1) - final_tumor9(idx,idy))./final_tumor9(1,1);
%         end
%     end

    % ������ͼ
    figure (1); 
    % ��������
    [Bb, G] = meshgrid(alpha_values, gamma_values); 
    surf(Bb, G, final_tumor9');
    colormap('jet');  % ������ɫӳ��
    shading interp;   % ʹ�ò�ֵ�����ɫ
    colorbar;         % �����ɫ��
    xlabel({'\alpha'},'FontWeight','bold','FontSize',14); 
    ylabel({'\gamma'},'FontWeight','bold','FontSize',14);
    zlabel({'E'},'FontWeight','bold','FontSize',14);

    figure (2); 
    imagesc(alpha_values, gamma_values, final_tumor9'); %��ά��ͼ
    colormap('jet');  % ������ɫӳ��
    colorbar;         % �����ɫ��
    xlabel({'\alpha'},'FontWeight','bold','FontSize',14); 
    ylabel({'\gamma'},'FontWeight','bold','FontSize',14);

end

function [T_all,Y_all] = jiujie(~)
global B1 B2 B3 E T total_time

y0=[B1 B2 B3 0 0 0 0 E T];
timestep = 30;     % ����30��
 
% �洢�����
T_all = [];
Y_all = [];

% �ֶ����
for t_start = 0:timestep:(total_time-timestep)
    tspan0 = [0,30];
    
    % ���� ode45 ������ֵ���
    [Tt, Y] = ode45(@(t, y)jiuwei(t, y), tspan0, y0);
    
    % ������洢���� 
    Y_all = [Y_all; Y];   % ��ÿ��ʱ���Yƴ�ӵ�Y_all
    T_all = [T_all; Tt+ t_start]; % ��ÿ��ʱ���Tƴ�ӵ�T_all
    
    % ���� y0  v=0.01 
    if ismember(t_start, [0, 90, 180, 270, 360, 450, 540])    
        jia1 = 0.1 + (0.3 - 0.1) * rand; %0.25  0.1
        y0 = Y(end, :) + [0 0 jia1 0 0 0 0 0 0]; 
    elseif ismember(t_start, [30, 120, 210, 300, 390, 480, 570])
        jia2 = 0.1 + (0.3 - 0.1) * rand; %0.2
        y0 = Y(end, :) + [0 jia2 0 0 0 0 0 0 0]; 
    elseif ismember(t_start, [60, 150, 240, 330, 420, 510])
        jia3 = 0.1 + (0.3 - 0.1) * rand; %0.2
        y0 = Y(end, :) + [jia3 0 0 0 0 0 0 0 0];
    end
    
end

end

function dy = jiuwei(t, y)  

global alpha m  A_c gamma alpha_A alpha_AA gamma_A    s n mu eta1 g1 r k e1     eta2 g2 e2 beta g3 e3      v b omega gamma_C t0 tau

    dy = zeros(9,1);
    dy(1) = alpha.*y(1)*exp(-v*(t-t0))-((y(7).^m)./(A_c.^m+y(7).^m)).*gamma.*y(1)-b*y(1)*y(6); 
    dy(2) = alpha.*y(2)*exp(-v*(t-t0))-((y(7).^m)./(A_c.^m+y(7).^m)).*gamma.*y(2)-b*y(2)*y(4);
    dy(3) = alpha.*y(3)*exp(-v*(t-t0))-((y(7).^m)./(A_c.^m+y(7).^m)).*gamma.*y(3)-b*y(3)*y(5);
    dy(4) = ((y(7).^m)./(A_c.^m+y(7).^m)).*gamma.*y(1)*omega - gamma_C*y(4);
    dy(5) = ((y(7).^m)./(A_c.^m+y(7).^m)).*gamma.*y(2)*omega - gamma_C*y(5);
    dy(6) = ((y(7).^m)./(A_c.^m+y(7).^m)).*gamma.*y(3)*omega - gamma_C*y(6);
    dy(7) = (alpha_A+alpha_AA.*((y(7).^m)./(A_c.^m+y(7).^m))).*(y(1)+y(2)+y(3))-gamma_A.*y(7);
    dy(8) = tau * (s-n*y(8)*y(9)-mu*y(8)+(eta1*y(8)*y(9))/(g1+y(9))+(eta2*y(8)*(y(1)+y(2)+y(3)))/(g2+(y(1)+y(2)+y(3))));
    dy(9) = tau * (r*y(9)*(1-y(9)/k)-e1*y(8)*y(9)-(e2*(beta*((y(7)^m)/(A_c^m+y(7)^m))*gamma*(y(1)+y(2)+y(3)))^2*y(9))/((g3)^2+(beta*((y(7)^m)/(A_c^m+y(7)^m))*gamma*(y(1)+y(2)+y(3)))^2)-e3*(y(1)+y(2)+y(3))*y(9));
end

function dy = erwei( t,y )
global s n mu eta1 g1 r k e1 tau

    dy = zeros(2,1);
    dy(1) = tau * (s-n*y(1)*y(2)-mu*y(1)+(eta1*y(1)*y(2))/(g1+y(2)));
    dy(2) = tau * (r*y(2)*(1-y(2)/k)-e1*y(1)*y(2));

end

function dy = siAHL( t,y )

global alpha m  A_c gamma alpha_A alpha_AA gamma_A    s n mu eta1 g1 r k e1     eta2 g2 e2 beta g3 e3      v t0 tau % b omega gamma_C
    
    dy = zeros(4,1);
    dy(1) = alpha.*y(1)*exp(-v*(t-t0))-((y(2).^m)./(A_c.^m+y(2).^m)).*gamma.*y(1);
    dy(2) = (alpha_A+alpha_AA.*((y(2).^m)./(A_c.^m+y(2).^m))).*y(1)-gamma_A.*y(2);
    dy(3) = tau * (s-n*y(3)*y(4)-mu*y(3)+(eta1*y(3)*y(4))/(g1+y(4))+(eta2*y(3)*y(1))/(g2+y(1)));
    dy(4) = tau * (r*y(4)*(1-y(4)/k)-e1*y(3)*y(4)-(e2*(beta*((y(2)^m)/(A_c^m+y(2)^m))*gamma*y(1))^2*y(4))/((g3)^2+(beta*((y(2)^m)/(A_c^m+y(2)^m))*gamma*y(1))^2)-e3*y(1)*y(4));
end
