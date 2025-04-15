function mc9

y0=[0.25 0.3 0.1 0 0 0 0 7e5 7e8];
total_time = 90;  % 总时间
timestep = 30;     % 步长30
 
% 存储求解结果
T_all = [];
Y_all = [];

% 分段求解
for t_start = 0:timestep:(total_time-timestep)
    tspan0 = [0,timestep];
    
    % 调用 ode45 进行数值求解
    [T, Y] = ode45(@(t, y)jiuwei(t, y), tspan0, y0);
% %     Y (size (Y,1),:) % format('long')
    
    % 将结果存储起来 
    Y_all = [Y_all; Y];   % 将每段时间的Y拼接到Y_all
    T_all = [T_all; T+ t_start]; % 将每段时间的T拼接到T_all
    
    % 更新 y0  v=0.01
    if t_start == 0 
        y0 = Y(end, :) + [0 0 0.12 0 0 0 0 0 0];
    elseif t_start ==30
        y0 = Y(end, :) + [0 0.17 0 0 0 0 0 0 0];
%     elseif t_start ==60 
%         y0 = Y(end, :) + [0.13 0 0 0 0 0 0 0 0];
    end
end

hold on
figure(1)
line([0 0], ylim, 'Color', [0.5 0.5 0.5], 'LineWidth', 5); % RGB 颜色值表示灰色
hold on
line([30 30], ylim, 'Color', [0.5 0.5 0.5], 'LineWidth', 5);
hold on
line([60 60], ylim, 'Color', [0.5 0.5 0.5], 'LineWidth', 5);
hold on
f1 = plot(T_all,Y_all(:,1), 'g','LineWidth',2);
hold on
f2 = plot(T_all,Y_all(:,2), 'b','LineWidth',2);
hold on
f3 = plot(T_all,Y_all(:,3), 'r','LineWidth',2);
hold on
legend([f1,f2,f3],{'B1','B2','B3'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);
xlabel('Time');
ylabel({'Bacterial'});
ylim([0 1]);
xlim([0 total_time]);

figure(2)
plot(T_all,Y_all(:,4),'g',T_all,Y_all(:,5),'b',T_all,Y_all(:,6),'r','LineWidth',2);
hold on;
legend({'C_{1}','C_{2}','C_{3}'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);
xlabel({'Time'},'FontWeight','bold','FontSize',14);
ylabel({'Colicins'},'FontWeight','bold','FontSize',14);
xlim([0 total_time]);

figure(3)
plot(T_all,Y_all(:,7),'k','LineWidth',2);
hold on;
xlabel({'Time'},'FontWeight','bold','FontSize',14);
ylabel({'AHL'},'FontWeight','bold','FontSize',14);
xlim([0 total_time]);

figure(4)
plot(T_all,Y_all(:,8),'k','LineWidth',2);
hold on;
xlabel({'Time'},'FontWeight','bold','FontSize',14);
ylabel({'Effector cells'},'FontWeight','bold','FontSize',14);
xlim([0 total_time]);

figure(5)
plot(T_all,Y_all(:,9),'k','LineWidth',2);
hold on;
xlabel({'Time'},'FontWeight','bold','FontSize',14);
ylabel({'Tumor cells'},'FontWeight','bold','FontSize',14);
xlim([0 total_time]);

figure(6)
plot(Y_all(:,8),Y_all(:,9),'k','LineWidth',2); 
xlabel({'Effector cells'},'FontWeight','bold','FontSize',14);
ylabel({'Tumor cells'},'FontWeight','bold','FontSize',14);
hold on 

% figure(7)
% plot(Y_all(:,4)+Y_all(:,5)+Y_all(:,6),Y_all(:,9),'k','LineWidth',2); 
% xlabel({'Colicins'},'FontWeight','bold','FontSize',14);
% ylabel({'Tumor cells'},'FontWeight','bold','FontSize',14);
% hold on 

end


function dy = jiuwei(t, y)
%     alpha = 1; m = 2; A_c = 1; gamma = 3; alpha_A = 0.4; alpha_AA = 8; gamma_A = 1;  
%     s = 1.3; n = 3.422e-10; mu = 0.0412; eta1 = 0.1245; g1 = 2.019e7;  r = 0.18; k = 5e8; e1 = 1.101e-7; %1.101e-7
%     eta2 = 0.2; g2 = 2.5e7; e2 = 0.4; beta = 5; g3 = 2e7; e3 = 0.3; %gamma = 4;极限环   s=1.3e4原    s=1.3双稳
%     v = 0.04; b = 0.1; omega = 5; gamma_C = 1; t0 = 0; 

    alpha = 1; m = 2; A_c = 1; gamma = 4; alpha_A = 0.4; alpha_AA = 8; gamma_A = 1;  
    s = 1.3e4; n = 3.422e-10; mu = 0.0412; eta1 = 0.1245; g1 = 2.019e7;  r = 0.18; k = 5e8; e1 = 1.101e-7; 
    eta2 = 0.2; g2 = 2.5e7; e2 = 0.4; beta = 5; g3 = 2e7; e3 = 0.3; %gamma = 4;极限环   n = 5e-10
    v = 0.01; b = 0.1; omega = 5; gamma_C = 1; t0 = 0;  tau = 0.05; 
   
    dy = zeros(9,1);  
    dy(1) =  alpha.*y(1)*exp(-v*(t-t0))-((y(7).^m)./(A_c.^m+y(7).^m)).*gamma.*y(1)-b*y(1)*y(6);
    dy(2) =  alpha.*y(2)*exp(-v*(t-t0))-((y(7).^m)./(A_c.^m+y(7).^m)).*gamma.*y(2)-b*y(2)*y(4);
    dy(3) =  alpha.*y(3)*exp(-v*(t-t0))-((y(7).^m)./(A_c.^m+y(7).^m)).*gamma.*y(3)-b*y(3)*y(5);
    dy(4) = ((y(7).^m)./(A_c.^m+y(7).^m)).*gamma.*y(1)*omega - gamma_C*y(4);
    dy(5) = ((y(7).^m)./(A_c.^m+y(7).^m)).*gamma.*y(2)*omega - gamma_C*y(5);
    dy(6) = ((y(7).^m)./(A_c.^m+y(7).^m)).*gamma.*y(3)*omega - gamma_C*y(6);
    dy(7) = (alpha_A+alpha_AA.*((y(7).^m)./(A_c.^m+y(7).^m))).*(y(1)+y(2)+y(3))-gamma_A.*y(7);
    dy(8) = tau * (s-n*y(8)*y(9)-mu*y(8)+(eta1*y(8)*y(9))/(g1+y(9))+(eta2*y(8)*(y(1)+y(2)+y(3)))/(g2+(y(1)+y(2)+y(3))));
    dy(9) = tau * (r*y(9)*(1-y(9)/k)-e1*y(8)*y(9)-(e2*(beta*((y(7)^m)/(A_c^m+y(7)^m))*gamma*(y(1)+y(2)+y(3)))^2*y(9))/((g3)^2+(beta*((y(7)^m)/(A_c^m+y(7)^m))*gamma*(y(1)+y(2)+y(3)))^2)-e3*(y(1)+y(2)+y(3))*y(9));
end



