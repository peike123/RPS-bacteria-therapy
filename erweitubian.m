function r
alpha = 1; m = 2; q_c = 1; gamma = 4; alpha_q = 0.4; alpha_Q = 8; gamma_q = 1;

y0=[0.1 0.5];

[T,Y]=ode45(@AHL,[0 90],y0);
[T1,Y1]=ode45(@tubianAHL,[0 90],y0);

figure(1)
plot(T,Y(:,1),'k',T1,Y1(:,1),'r','LineWidth',2);
hold on;
legend({'Strain without mutation','Strain with mutation'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);
xlabel({'Time'},'FontWeight','bold','FontSize',14);
ylabel({'Bacteria'},'FontWeight','bold','FontSize',14);
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',1.5);  %坐标轴的字体和大小
xlim([0 90]);
% figure(2)
% plot(T,Y(:,2),'k',T1,Y1(:,2),'b','LineWidth',2);
% hold on;
% legend({'\rho_R','\rho_S','\rho_P','\rho_0'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);
% xlabel({'Time'},'FontWeight','bold','FontSize',14);
% ylabel({'AHL'},'FontWeight','bold','FontSize',14);
% set(gca,'FontWeight','bold','FontSize',14,'LineWidth',1.5);  %坐标轴的字体和大小
end


function dy = AHL( t,y )
    alpha = 1; m = 2; q_c = 1; gamma = 4; alpha_q = 0.4; alpha_Q = 8; gamma_q = 1; 

    dy = zeros(2,1);
    dy(1) = alpha.*y(1)-((y(2).^m)./(q_c.^m+y(2).^m)).*gamma.*y(1);
    dy(2) = (alpha_q+alpha_Q.*((y(2).^m)./(q_c.^m+y(2).^m))).*y(1)-gamma_q.*y(2);

end

function dy = tubianAHL( t,y )
%     alpha = 1; m = 2; q_c = 1; gamma = 20; alpha_q = 0.4; alpha_Q = 8; gamma_q = 1;    
    alpha = 1; m = 2; q_c = 1; gamma = 4; alpha_q = 0.4; alpha_Q = 8; gamma_q = 1; 
    v = 0.05; t0 = 0;

    dy = zeros(2,1);
    dy(1) = alpha.*y(1) *exp(-v*(t-t0))-((y(2).^m)./(q_c.^m+y(2).^m)).*gamma.*y(1);
    dy(2) = (alpha_q+alpha_Q.*((y(2).^m)./(q_c.^m+y(2).^m))).*y(1)-gamma_q.*y(2);

end


