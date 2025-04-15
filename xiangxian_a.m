function r
global B B1 B2 B3 E T total_time
global alpha m  A_c gamma alpha_A alpha_AA gamma_A    s n mu eta1 g1 r k e1     eta2 g2 e2 beta g3 e3      v b omega gamma_C t0  tau

total_time = 600;  % 总时间
% B = 0.3; B1 = 0.25; B2 = 0.3; B3 = 0.1; E = 1e5; T = 1e7;
m = 2; A_c = 1; gamma = 4; alpha_A = 0.4; alpha_AA = 8; gamma_A = 1;  
s = 1.3e4; n = 3.422e-10; mu = 0.0412; eta1 = 0.1245; g1 = 2.019e7;  r = 0.5; k = 5e8; e1 = 1.101e-7; 
eta2 = 0.2; g2 = 2.5e7; e2 = 0.4; beta = 5; g3 = 2e7; e3 = 0.8;
b = 0.1; v = 0.05; omega = 5; gamma_C = 1; t0 = 0; tau = 0.05;
% alpha = 1; 
alpha_value= 0:0.2:2;

total_simulations = 50;
tumor_cells9 = [length(alpha_value), total_simulations];

for j = 1:length(alpha_value)
    for i = 1:total_simulations
        alpha = alpha_value(j);
        
        B = 0.27 + (0.33 - 0.27) * rand;
        B1 = 0.225 + (0.275 - 0.225) * rand;
        B2 = 0.27 + (0.33 - 0.27) * rand;
        B3 = 0.09 + (0.11 - 0.09) * rand; %0.25 0.3 0.1
        
        E = randi([1e5 1.5e5]);  %生成[]之间的整数
        T = randi([1e7 1.5e7]);
%         E = randi([1e5 2e5]);  %生成[]之间的整数
%         T = randi([1e7 2e7]);
        
        [T9,Y9] = jiujie();
        
        % 找到肿瘤细胞数 y(9) 超过阈值的时间
        tumor_cells9(j,i) = Y9(end, 9);  % 取出肿瘤细胞数 y(9)
       
    end
end
    

    % 绘制箱线图
    figure (1); 
    h = boxplot(tumor_cells9'); 
    set(h, 'Linewidth', 1.2);  %设置线宽
    % 设置箱体颜色
    set(findobj(h, 'Tag', 'Box'), 'Color', [1, 0, 0, 0.5]); % 自定义颜色 
    % 设置中位数线颜色
    set(findobj(h, 'Tag', 'Median'), 'Color', [1, 0, 0, 1]); % 中位数为红色
%     % 设置须线颜色
%     set(findobj(h, 'Tag', 'Whisker'), 'Color', [1, 0, 0, 1]); 
    xlabel('alpha','FontWeight','bold','FontSize',14);
    ylabel('Number of tumor cells at 600 hr','FontWeight','bold','FontSize',14);
    grid on;
%     x = 1:5;
%     set(gca,'Xtick',x);
%     set(gca,'xticklabel',{'0.05','0.1','0.15','0.2','0.25'});
    
end

function [T_all,Y_all] = jiujie(~)
global B1 B2 B3 E T total_time

y0=[B1 B2 B3 0 0 0 0 E T];
timestep = 30;     % 步长30秒
 
% 存储求解结果
T_all = [];
Y_all = [];

% 分段求解
for t_start = 0:timestep:(total_time-timestep)
    tspan0 = [0,30];
    
    % 调用 ode45 进行数值求解
    [Tt, Y] = ode45(@(t, y)jiuwei(t, y), tspan0, y0);
    
    % 将结果存储起来 
    Y_all = [Y_all; Y];   % 将每段时间的Y拼接到Y_all
    T_all = [T_all; Tt+ t_start]; % 将每段时间的T拼接到T_all
    
    % 更新 y0  v=0.01 
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

