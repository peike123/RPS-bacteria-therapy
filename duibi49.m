function mc9

y4=[0.2 0 2e5 5e7]; %7e5 7e8
% y4=[0.5 0 1e3 1e4];
y9=[0.2 0.3 0.1 0 0 0 0 2e5 5e7];
total_time = 150;  % 总时间
timestep = 30;     % 步长30秒
 
% 存储求解结果
T_all = [];
Y_all = [];

% 分段求解
for t_start = 0:timestep:(total_time-timestep)
    % 时间段 tspan = [t_start, t_start + timestep];
    tspan0 = [0,30];
    
    % 调用 ode45 进行数值求解
    [T, Y] = ode45(@(t, y)jiuwei(t, y), tspan0, y9);
    Y (size (Y,1),:)
    
    % 将结果存储起来 
    Y_all = [Y_all; Y];   % 将每段时间的Y拼接到Y_all
    T_all = [T_all; T+ t_start]; % 将每段时间的T拼接到T_all
    
    % 更新 y0
    if t_start == 0                            %e=0.3    %e=0.5
        y9 = Y(end, :) + [0 0 0.25 0 0 0 0 0 0]; %0.13   0.2
    elseif t_start ==30
        y9 = Y(end, :) + [0 0.2 0 0 0 0 0 0 0]; %0.07   0.1
    elseif t_start == 60
        y9 = Y(end, :) + [0.2 0 0 0 0 0 0 0 0]; %0.025   0.1
    elseif t_start == 90
        y9 = Y(end, :) + [0 0 0.1 0 0 0 0 0 0]; %0.007   0.05
    end

%         if ismember(t_start, [0, 90, 180, 270, 360, 450, 540])    
%             jia1 = 0.1 + (0.3 - 0.1) * rand; %0.25  0.1
%             y9 = Y(end, :) + [0 0 jia1 0 0 0 0 0 0]; 
%         elseif ismember(t_start, [30, 120, 210, 300, 390, 480, 570])
%             jia2 = 0.1 + (0.3 - 0.1) * rand; %0.2
%             y9 = Y(end, :) + [0 jia2 0 0 0 0 0 0 0]; 
%         elseif ismember(t_start, [60, 150, 240, 330, 420, 510])
%             jia3 = 0.1 + (0.3 - 0.1) * rand; %0.2
%             y9 = Y(end, :) + [jia3 0 0 0 0 0 0 0 0];
%         end

end

[T4,Y4]=ode45(@siwei,[0 total_time],y4);

hold on
figure(1)
line([0 0], ylim, 'Color', [0.5 0.5 0.5], 'LineWidth', 5); % RGB 颜色值表示灰色
hold on
line([30 30], ylim, 'Color', [0.5 0.5 0.5], 'LineWidth', 5);
hold on
line([60 60], ylim, 'Color', [0.5 0.5 0.5], 'LineWidth', 5);
hold on
line([90 90], ylim, 'Color', [0.5 0.5 0.5], 'LineWidth', 5);
hold on
line([120 120], ylim, 'Color', [0.5 0.5 0.5], 'LineWidth', 5);
hold on
f1 = plot(T_all,Y_all(:,1), 'g','LineWidth',2);
hold on
f2 = plot(T_all,Y_all(:,2), 'b','LineWidth',2);
hold on
f3 = plot(T_all,Y_all(:,3), 'r','LineWidth',2);
hold on
f4 = plot(T4,Y4(:,1),'m','LineWidth',2);
legend([f1, f2,f3 ,f4],{'B1','B2','B3','Single strain'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);
xlabel('Time');
ylabel({'Bacterial'});
ylim([0 1]);
xlim([0 total_time]);
hold on

figure(2)
plot(T_all,Y_all(:,7),'k',T4,Y4(:,2),'m','LineWidth',2);
hold on;
legend({'Nine-dimensional','Four-dimensional'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);
xlabel({'Time'},'FontWeight','bold','FontSize',14);
ylabel({'AHL'},'FontWeight','bold','FontSize',14);
xlim([0 total_time]);

figure(3)
plot(T_all,Y_all(:,8),'k',T4,Y4(:,3),'m','LineWidth',2);
hold on;
xlabel({'Time'},'FontWeight','bold','FontSize',14);
ylabel({'Effector cells'},'FontWeight','bold','FontSize',14); 
legend({'Nine-dimensional','Four-dimensional'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);
xlim([0 total_time]);

figure(4)
plot(T_all,Y_all(:,9),'k',T4,Y4(:,4),'m','LineWidth',2); 
hold on;
xlabel({'Time'},'FontWeight','bold','FontSize',14);
ylabel({'Tumor cells'},'FontWeight','bold','FontSize',14); 
legend({'Nine-dimensional','Four-dimensional'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);
xlim([0 total_time]);

figure(5)
plot(Y_all(:,8),Y_all(:,9),'k','LineWidth',2); 
hold on;
plot(Y4(:,3),Y4(:,4),'m','LineWidth',2);
xlabel({'Effector cells'},'FontWeight','bold','FontSize',14);
ylabel({'Tumor cells'},'FontWeight','bold','FontSize',14);
legend({'Nine-dimensional','Four-dimensional'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);
hold on 
end


function dy = jiuwei(t, y )
    alpha = 1; m = 2; A_c = 1; gamma = 4; alpha_A = 0.4; alpha_AA = 8; gamma_A = 1;  
    s = 1.3e4; n = 3.422e-10; mu = 0.0412; eta1 = 0.1245; g1 = 2.019e7;  r = 0.18; k = 5e8; e1 = 1.101e-7; 
    eta2 = 0.2; g2 = 2.5e7; e2 = 0.4; beta = 5; g3 = 2e7; e3 = 0.55; %gamma = 4;极限环 e3 = 0.5;
    v = 0.05; b = 0.05; omega = 5; gamma_C = 1; t0 = 0; tau = 0.05; 
   
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

function dy = siwei( t,y )
    alpha = 1; m = 2; A_c = 1; gamma = 4; alpha_A = 0.4; alpha_AA = 8; gamma_A = 1;  
    s = 1.3e4; n = 3.422e-10; mu = 0.0412; eta1 = 0.1245; g1 = 2.019e7;  r = 0.18; k = 5e8; e1 = 1.101e-7; 
%     s = 1.3; n = 5e-10; mu = 0.0412; eta1 = 0.1245; g1 = 2.019e7; r = 0.18; k = 5e8; e1 = 1.101e-7; 
    eta2 = 0.2; g2 = 2.5e7; e2 = 0.4; beta = 5; g3 = 2e7; e3 = 0.5; 
    v = 0.05; tau = 0.05; %e3 = 0.3

    dy = zeros(4,1);
    dy(1) = alpha.*y(1)*exp(-v*t)-((y(2).^m)./(A_c.^m+y(2).^m)).*gamma.*y(1);
    dy(2) = (alpha_A+alpha_AA.*((y(2).^m)./(A_c.^m+y(2).^m))).*y(1)-gamma_A.*y(2);
    dy(3) = tau * (s-n*y(3)*y(4)-mu*y(3)+(eta1*y(3)*y(4))/(g1+y(4))+(eta2*y(3)*y(1))/(g2+y(1)));
    dy(4) = tau * (r*y(4)*(1-y(4)/k)-e1*y(3)*y(4)-(e2*(beta*((y(2)^m)/(A_c^m+y(2)^m))*gamma*y(1))^2*y(4))/((g3)^2+(beta*((y(2)^m)/(A_c^m+y(2)^m))*gamma*y(1))^2)-e3*y(1)*y(4));
end

