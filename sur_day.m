function r
global B B1 B2 B3 E T total_time
global alpha m  A_c gamma alpha_A alpha_AA gamma_A    s n mu eta1 g1 r k e1     eta2 g2 e2 beta g3 e3      v b omega gamma_C t0  tau

total_time = 600;  % ��ʱ��
alpha = 1; gamma = 4; alpha_A = 0.4; alpha_AA = 8; 
s = 1.3e4; g1 = 2.019e7; eta2 = 0.2; g2 = 2.5e7; e2 = 0.4; 
beta = 5; g3 = 2e7; b = 0.1; omega = 5; gamma_C = 1; t0 = 0; tau = 0.05;
v = 0.05;  % v 0.01-0.05    % e3 0.1-0.5 

total_simulations = 50;
death_threshold = 3e8;  % ����ϸ�����ﵽ2e8ʱ����
death_times4 = NaN(total_simulations, 1);  % ��¼ÿ�����������ʱ��
death_times9 = NaN(total_simulations, 1);  % ��¼ÿ�����������ʱ��
    
    for i = 1:total_simulations
        B = 0.27 + (0.33 - 0.27) * rand;  %  10%
        B1 = 0.225 + (0.275 - 0.225) * rand;
        B2 = 0.27 + (0.33 - 0.27) * rand;
        B3 = 0.09 + (0.11 - 0.09) * rand; %0.25 0.3 0.1
          
        E = randi([1e5 1e6]);  %����[]֮�������
        T = randi([1e7 2e8]);
      
%         E = randi([0.5e5 9.5e5]);  %����[]֮�������  90%
%         T = randi([0.5e7 9.5e7]);
        
        gamma_A = 0.5 + (1.5 - 0.5) * rand;  % 1   50%
        m = 1 + (3 - 1) * rand;  %2
        A_c = 0.5 + (1.5 - 0.5) * rand;  %1
        e3 = 0.4 + (1.2 - 0.4) * rand; %0.3����0.8����1
        
        k = 4.5e8 + (5.5e8 - 4.5e8) * rand;  %5e8    10%
        r = 0.45 + (0.55 - 0.45) * rand; % 0.18����0.5
        
        n = 2.6e-10 + (4.2e-10 - 2.6e-10) * rand;  % 3.422e-10;    25%
        mu = 0.03 + (0.05 - 0.03) * rand;  %0.0412;
        eta1 = 0.09 + (0.15 - 0.09) * rand;  %0.1245;
        e1 = 0.83e-7 + (1.37e-7 - 0.83e-7) * rand;  %1.101e-7;

%         gamma_A = 0.5 + (1.5 - 0.5) * rand;  %1    
%         m = 1.5 + (2.5 - 1.5) * rand;  %2 
%         A_c = 0.5 + (1.5 - 0.5) * rand;  %1
%         e3 = 0.5 + (1.5 - 0.5) * rand; %1  0.3    50%
%         k = 4.5e8 + (5.5e8 - 4.5e8) * rand;  %5e8;
%         r = 0.4 + (0.7 - 0.4) * rand; %0.5   0.18  10%
%         n = 3.41e-10 + (3.44e-10 - 3.41e-10) * rand;  % 3.422e-10; 
%         mu = 0.03 + (0.05 - 0.03) * rand;  %0.0412; 
%         eta1 = 0.11 + (0.14 - 0.11) * rand;  %0.1245;   
%         e1 = 1.1e-7 + (1.2e-7 - 1.1e-7) * rand;  %1.101e-7; 
        
        y2 = [E T];
        [T2,Y2] = ode45(@erwei,[0 total_time],y2);
        
        y4=[B 0 E T]; 
        [T4,Y4] = ode45(@siAHL,[0 total_time],y4);
        
        [T9,Y9] = jiujie();
        
        % �ҵ�����ϸ���� y(2) ������ֵ��ʱ��
        tumor_cells2 = Y2(:, 2);  % ȡ������ϸ���� y(2)
        death_idx2 = find(tumor_cells2 > death_threshold, 1);   % �ҵ���һ��������ֵ��λ��

        if isempty(death_idx2)
            death_times2(i) = total_time+1;  % �������ϸ����δ������ֵ����Ϊδ����
        else
            death_times2(i) = T2(death_idx2);  % ��¼����ʱ��
        end

        % �ҵ�����ϸ���� y(4) ������ֵ��ʱ��
        tumor_cells4 = Y4(:, 4);  % ȡ������ϸ���� y(4)
        death_idx4 = find(tumor_cells4 > death_threshold, 1);   % �ҵ���һ��������ֵ��λ��

        if isempty(death_idx4)
            death_times4(i) = total_time+1;  % �������ϸ����δ������ֵ����Ϊδ����
        else
            death_times4(i) = T4(death_idx4);  % ��¼����ʱ��
        end
        
        % �ҵ�����ϸ���� y(9) ������ֵ��ʱ��
        tumor_cells9 = Y9(:, 9);  % ȡ������ϸ���� y(9)
        death_idx9 = find(tumor_cells9 > death_threshold, 1);  % �ҵ���һ��������ֵ��λ��

        if isempty(death_idx9)
            death_times9(i) = total_time+1;  % �������ϸ����δ������ֵ����Ϊδ����
        else
            death_times9(i) = T9(death_idx9);  % ��¼����ʱ��
        end
    end
    
    valid_death_times2 = death_times2(~isnan(death_times2)); % ���˵� NaN ֵ��ȷ��ֻ������Ч������ʱ��
    death_counts2 = zeros(1, 6);    % Ԥ���������� ÿ��������
    cumulatively_death_counts2 = zeros(1, total_time);    % Ԥ���������� �ۼ�������
    for t = 1:total_time     % �����ۻ���������
        cumulatively_death_counts2(t) = sum(valid_death_times2 <= t);        % ���㵽��ǰʱ�� t ���ۼ���������  1*75 
    end 
    for j = 1:1:6
        death_counts2(j) = sum(valid_death_times2 > (j*100-100) & valid_death_times2 <= j*100);        % ���� ��ǰʱ�� t ����������  1*75
    end
   
    valid_death_times4 = death_times4(~isnan(death_times4)); % ���˵� NaN ֵ��ȷ��ֻ������Ч������ʱ��
    death_counts4 = zeros(1, 6);    % Ԥ���������� ÿ��������
    cumulatively_death_counts4 = zeros(1, total_time);    % Ԥ���������� �ۼ�������
    for t = 1:total_time     % �����ۻ���������
        cumulatively_death_counts4(t) = sum(valid_death_times4 <= t);        % ���㵽��ǰʱ�� t ���ۼ���������  1*75 
    end
    for j = 1:1:6
        death_counts4(j) = sum(valid_death_times4 > (j*100-100) & valid_death_times4 <= j*100);        % ���� ��ǰʱ�� t ����������  1*75
    end
    
    valid_death_times9 = death_times9(~isnan(death_times9)); % ���˵� NaN ֵ��ȷ��ֻ������Ч������ʱ��
    death_counts9 = zeros(1, 6);    % Ԥ���������� ÿ��������
    cumulatively_death_counts9 = zeros(1, total_time);    % Ԥ���������� �ۼ�������
    for t = 1:total_time     % �����ۻ���������
        cumulatively_death_counts9(t) = sum(valid_death_times9 <= t);        % ���㵽��ǰʱ�� t ���ۼ���������  1*75 
    end
    for j = 1:1:6
        death_counts9(j) = sum(valid_death_times9 > (j*100-100) & valid_death_times9 <= j*100);        % ���� ��ǰʱ�� t ����������  1*75
    end
    
    survival_simulations2 = total_simulations - cumulatively_death_counts2;   %�������   1*75
    survival_rate2 = survival_simulations2./total_simulations;    %������ 1*75
    
    survival_simulations4 = total_simulations - cumulatively_death_counts4;   %�������   1*75
    survival_rate4 = survival_simulations4./total_simulations;    %������ 1*75
    
    survival_simulations9 = total_simulations - cumulatively_death_counts9;   %�������   1*75
    survival_rate9 = survival_simulations9./total_simulations;    %������ 1*75

    % ���ɸ����ܼ���ʱ������
    dense_time2 = linspace(1, total_time, 750);  % ��ԭʼ����С���ʱ���֮������ 20 ����
    dense_survival_rate2 = interp1(1:total_time, survival_rate2, dense_time2, 'linear');  % ʹ�����Բ�ֵ����������ʵĲ�ֵ
    
    dense_time4 = linspace(1, total_time, 750);  % ��ԭʼ����С���ʱ���֮������ 20 ����
    dense_survival_rate4 = interp1(1:total_time, survival_rate4, dense_time4, 'linear');  % ʹ�����Բ�ֵ����������ʵĲ�ֵ
    
    dense_time9 = linspace(1, total_time, 750);  % ��ԭʼ����С���ʱ���֮������ 20 ����
    dense_survival_rate9 = interp1(1:total_time, survival_rate9, dense_time9, 'linear');  % ʹ�����Բ�ֵ����������ʵĲ�ֵ

    %�������� 
    z_value = 1.645;  % 90%���������Zֵ 1.645
    se2 = sqrt((survival_rate2 .* (1 - survival_rate2)) / total_simulations);  % ��׼������
    ci_lower2 = max(survival_rate2 - z_value * se2, 0);  % ���޲���С��0
    ci_upper2 = min(survival_rate2 + z_value * se2, 1);  % ���޲��ܴ���1

    se4 = sqrt((survival_rate4 .* (1 - survival_rate4)) / total_simulations);  % 4άģ�͵ı�׼���
    ci_lower4 = max(survival_rate4 - z_value * se4, 0);
    ci_upper4 = min(survival_rate4 + z_value * se4, 1);

    se9 = sqrt((survival_rate9 .* (1 - survival_rate9)) / total_simulations);  % 9άģ�͵ı�׼���
    ci_lower9 = max(survival_rate9 - z_value * se9, 0);
    ci_upper9 = min(survival_rate9 + z_value * se9, 1);
    
    % ������������
    figure (1); 
    plot([0,dense_time2], [1,dense_survival_rate2], 'LineWidth', 1.5, 'Color', [0, 0, 1, 0.2]);  % ��ֵ�������
    hold on;
    plot([0,dense_time4], [1,dense_survival_rate4], 'LineWidth', 1.5, 'Color', [1, 0, 0, 0.2]);  % ��ֵ�������
    hold on;
    plot([0,dense_time9], [1,dense_survival_rate9], 'LineWidth', 1.5, 'Color', [0, 1, 0, 0.2]);  % ��ֵ�������
    hold on;
    fill([[0:1:total_time] fliplr([0:1:total_time])], [1 ci_upper2 fliplr([1 ci_lower2])], 'b', 'FaceAlpha', 0.05, 'EdgeColor', 'b','EdgeAlpha', 0.4); % �����������
    hold on;
    fill([[0:1:total_time] fliplr([0:1:total_time])], [1 ci_upper4 fliplr([1 ci_lower4])], 'r', 'FaceAlpha', 0.05, 'EdgeColor', 'r','EdgeAlpha', 0.4);  % �����������
    hold on;
    fill([[0:1:total_time] fliplr([0:1:total_time])], [1 ci_upper9 fliplr([1 ci_lower9])], 'g', 'FaceAlpha', 0.05, 'EdgeColor', 'g','EdgeAlpha', 0.4);  % �����������
    legend({'No treatment','Single strain','Three strains'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);
    xlabel('Time (hr)','FontWeight','bold','FontSize',14);
    ylabel('Survival rate (%)','FontWeight','bold','FontSize',14);
    grid on;
    ylim([0 1]);
    xlim([0 total_time]);
    
    % ���ƴ����
    figure (2); 
    x = 1:6;
%     death_counts2 = [30, 20, 0, 0, 0, 0];
%     death_counts4 = [7, 30, 11, 5, 0, 3];
%     death_counts9 = [5, 10, 0, 1, 0, 0];
    fig = bar(x, [death_counts2; death_counts4; death_counts9]');
    
    % ���ò�ͬ��ɫ
    fig(1).FaceColor = [0 0 1];  % ��ɫ
    fig(2).FaceColor = [1 0 0];  % ��ɫ
    fig(3).FaceColor = [0 1 0];  % ��ɫ

    % ����͸����
    fig(1).FaceAlpha = 0.15;  % ���õ�1������͸����
    fig(2).FaceAlpha = 0.15;  % ���õ�2������͸����
    fig(3).FaceAlpha = 0.15;  % ���õ�3������͸����

%     % ��ÿ����״ͼ�Ϸ������ֵ��ǩ
%     for ii = 1:3  % 3�����Ʋ���
%         for jj = 1:6  % ÿ��ʱ���
%             % ��ȡÿ�����ӵ�x�����yֵ
%             x_pos = x(jj) + (ii - 2) * 0.2;  % ����xλ�ã�ȷ����ֵλ�������Ϸ�
%             y_value = [death_counts2(jj), death_counts4(jj), death_counts9(jj)];  % ÿ�����������
%             text(x_pos, y_value(ii) - 1, num2str(y_value(ii)), 'HorizontalAlignment', 'center', ...
%                  'VerticalAlignment', 'bottom', 'FontSize', 14, 'FontWeight', 'bold');
%         end
%     end

    xlabel('Time (hr)','FontWeight','bold','FontSize',14);
    ylabel('Number of dead mouse','FontWeight','bold','FontSize',14);
    legend({'No treatment','Single strain','Three strains'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);
    grid on;
    
    set(gca,'Xtick',x);
    set(gca,'xticklabel',{'(0, 100]','(100, 200]','(200, 300]','(300, 400]','(400, 500]','(500, 600]'});
    
    figure (3);
    f1 = plot(T9,Y9(:,1), 'g','LineWidth',2);
    hold on
    f2 = plot(T9,Y9(:,2), 'b','LineWidth',2);
    hold on
    f3 = plot(T9,Y9(:,3), 'r','LineWidth',2);
    hold on
    f4 = plot(T4,Y4(:,1),'m','LineWidth',2);
    xlabel('Time','FontWeight','bold','FontSize',14);
    ylabel({'Bacterial'},'FontWeight','bold','FontSize',14);
    legend([f1, f2, f3 ,f4],{'B1','B2','B3','Single strain'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);


  

end

function dy = siAHL( t,y )

global alpha m  A_c gamma alpha_A alpha_AA gamma_A    s n mu eta1 g1 r k e1     eta2 g2 e2 beta g3 e3      v t0 tau % b omega gamma_C
    
    dy = zeros(4,1);
    dy(1) = alpha.*y(1)*exp(-v*(t-t0))-((y(2).^m)./(A_c.^m+y(2).^m)).*gamma.*y(1);
    dy(2) = (alpha_A+alpha_AA.*((y(2).^m)./(A_c.^m+y(2).^m))).*y(1)-gamma_A.*y(2);
    dy(3) = tau * (s-n*y(3)*y(4)-mu*y(3)+(eta1*y(3)*y(4))/(g1+y(4))+(eta2*y(3)*y(1))/(g2+y(1)));
    dy(4) = tau * (r*y(4)*(1-y(4)/k)-e1*y(3)*y(4)-(e2*(beta*((y(2)^m)/(A_c^m+y(2)^m))*gamma*y(1))^2*y(4))/((g3)^2+(beta*((y(2)^m)/(A_c^m+y(2)^m))*gamma*y(1))^2)-e3*y(1)*y(4));
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
%     if t_start == 0                            %e=0.3    %e=0.5
%         jia1 = 0.2 + (0.3 - 0.2) * rand; %0.25
%         y0 = Y(end, :) + [0 0 jia1 0 0 0 0 0 0]; %0.13   0.2
%     elseif t_start ==30
%         jia2 = 0.15 + (0.25 - 0.15) * rand; %0.2
%         y0 = Y(end, :) + [0 jia2 0 0 0 0 0 0 0]; %0.07   0.1
%     elseif t_start == 60
%         jia3 = 0.15 + (0.25 - 0.15) * rand; %0.2
%         y0 = Y(end, :) + [jia3 0 0 0 0 0 0 0 0]; %0.025   0.1
%     elseif t_start == 90
%         jia4 = 0.05 + (0.15 - 0.05) * rand; %0.1
%         y0 = Y(end, :) + [0 0 jia4 0 0 0 0 0 0]; %0.007   0.05
%     elseif t_start >= 120
%         y0 = Y(end, :);
%     end

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

