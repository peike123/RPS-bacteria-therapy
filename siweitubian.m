function r
%  alpha = 1; m = 2; A_c = 1; gamma = 20; alpha_A = 0.4; alpha_AA = 8; gamma_A = 1;  
%  s = 1.3; n = 5e-10; mu = 0.0412; eta1 = 0.1245; g1 = 2.019e7;  r = 0.18; k = 5e8; e1 = 1.101e-7; %s = 1.3e4; n = 3.422e-10; 
%  eta2 = 0.2; g2 = 2.5e7; e2 = 0.4; beta = 5; g3 = 2e7; e3 = 0.3; v = 0; %gamma = 4;极限环 20稳定平衡点
%   
 alpha = 1; m = 2; A_c = 1; gamma = 4; alpha_A = 0.4; alpha_AA = 8; gamma_A = 1;  
 s = 1.3e4; n = 3.422e-10; mu = 0.0412; eta1 = 0.1245; g1 = 2.019e7;  r = 0.18; k = 5e8; e1 = 1.101e-7; 
 eta2 = 0.2; g2 = 2.5e7; e2 = 0.4; beta = 5; g3 = 2e7; e3 = 0.3; v = 0.05; tau = 0.05;
   
Time = 15000;
y0 = [0.5 0 1e3 1e4];
% y0=[0.5 0 2e4 7e8];
% y0=[0.5 0 1e2 1e8];
[T,Y]=ode45(@si,[0 Time],y0); 
[T1,Y1]=ode45(@situ,[0 Time],y0);  %V=0.05  T=20


figure(1)
plot(T,Y(:,1),'b','LineWidth',2);% .*10.^(-5)
hold on;
plot(T1,Y1(:,1),'k','LineWidth',2);
xlabel({'Time'},'FontWeight','bold','FontSize',14);
ylabel({'Bacteria'},'FontWeight','bold','FontSize',14);
% 
% figure(2)
% plot(T,Y(:,2),'LineWidth',2);
% hold on;
% xlabel({'Time'});
% ylabel({'AHL'});

figure(3)
plot(Y(:,1),Y(:,2),'b','LineWidth',2); 
hold on
plot(Y1(:,1),Y1(:,2),'k','LineWidth',2); 
xlabel({'Bacteria'},'FontWeight','bold','FontSize',14);
ylabel({'AHL'},'FontWeight','bold','FontSize',14);
hold on
%画平衡点
syms x1 x2 x3 x4; % 定义x1 x2 是未知量
eqn1=[alpha*x1-((x2^m)/(A_c^m+x2^m))*gamma*x1==0, (alpha_A+alpha_AA*((x2^m)/(A_c^m+x2^m)))*x1-gamma_A*x2==0]; % 正平衡点满足的方程组,
vars=[x1,x2]; % 定义求解的未知量
[solX1,solX2]=solve(eqn1,vars,'Real',true); % 求解平衡点
id = solX1<0; %去除<0的
solX1(id) = [];
solX2(id) = [];
plot(solX1,solX2,'r*')
hold on
legend({'Strain without mutation','Strain with mutation','Equilibrium point without mutation'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);

figure(4)
plot(T,Y(:,3),'b','LineWidth',2);
hold on;
plot(T1,Y1(:,3),'k','LineWidth',2);
xlabel({'Time'},'FontWeight','bold','FontSize',14);
ylabel({'Effector cells'},'FontWeight','bold','FontSize',14);
legend({'Strain without mutation','Strain with mutation'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);

figure(5)
plot(T,Y(:,4),'b','LineWidth',2);
hold on;
plot(T1,Y1(:,4),'k','LineWidth',2);
xlabel({'Time'},'FontWeight','bold','FontSize',14);
ylabel({'Tumor cells'},'FontWeight','bold','FontSize',14);
legend({'Strain without mutation','Strain with mutation'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);

figure(6)
plot(Y(:,3),Y(:,4),'b','LineWidth',2); 
hold on
plot(Y1(:,3),Y1(:,4),'k','LineWidth',2); 
xlabel({'Effector cells'},'FontWeight','bold','FontSize',14);
ylabel({'Tumor cells'},'FontWeight','bold','FontSize',14);
hold on 
%画正平衡点
AA = ((alpha*A_c^m)/(gamma-alpha))^(1/m);
BB = ((gamma*gamma_A)* AA)/(gamma*alpha_A+alpha*alpha_AA);
syms x3 x4;
eqn1=[s-n*x3*x4-mu*x3+(eta1*x3*x4)/(g1+x4)+(eta2*x3*BB)/(g2+BB)==0,...
    r*x4*(1-x4/k)-e1*x3*x4-(e2*(beta*((AA^m)/(A_c^m+AA^m))*gamma*BB)^2*x4)/((g3)^2+(beta*((AA^m)/(A_c^m+AA^m))*gamma*BB)^2)-e3*BB*x4==0]; % 正平衡点满足的方程组,
vars=[x3,x4]; % 定义求解的未知量
[solX3,solX4]=solve(eqn1,vars); % 求解平衡点
plot(solX3,solX4,'r*')
hold on
legend({'Without mutation','With mutation','Equilibrium point without mutation'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);
% %向量场
% [y3,y4]=meshgrid(linspace(0,8*10^7));
% y1 = ( (alpha*A_c^m)/(gamma-alpha) )^(1/m);
% y2 = (gamma*gamma_A)*( (alpha*A_c^m)/(gamma-alpha) )^(1/m)/(gamma*alpha_A+alpha*alpha_AA);
% streamslice(y3,y4,s-n.*y3.*y4-mu.*y3+(eta1.*y3.*y4)./(g1+y4)+(eta2.*y3.*y1)./(g2+y1) ,r.*y4.*(1-y4./k)-e1.*y3.*y4-(e2.*(beta*((y2.^m)/(A_c.^m+y2.^m)).*gamma.*y1).^2.*y4)/((g3).^2+(beta.*((y2.^m)/(A_c.^m+y2.^m)).*gamma.*y1).^2)-e3.*y1.*y4 );

figure(7)
plot3(Y(:,1),Y(:,3),Y(:,4),'LineWidth',2);
hold on
plot3(Y1(:,1),Y1(:,3),Y1(:,4),'LineWidth',2);
xlabel({'Strains'},'FontWeight','bold','FontSize',14);
ylabel({'Effector cells'},'FontWeight','bold','FontSize',14);
zlabel({'Tumor cells'},'FontWeight','bold','FontSize',14)
legend({'Without mutation','With mutation'},'Location','northeast','FontWeight','bold','FontSize',14,'LineWidth',1.5);
end

function dy = situ( t,y )
    alpha = 1; m = 2; A_c = 1; gamma = 4; alpha_A = 0.4; alpha_AA = 8; gamma_A = 1;  
    s = 1.3e4; n = 3.422e-10; mu = 0.0412; eta1 = 0.1245; g1 = 2.019e7;  r = 0.18; k = 5e8; e1 = 1.101e-7; 
    eta2 = 0.2; g2 = 2.5e7; e2 = 0.4; beta = 5; g3 = 2e7; e3 = 0.3; v = 0.05; tau = 0.05; 
  
    dy = zeros(4,1);
    dy(1) = alpha.*y(1)*exp(-v*t)-((y(2).^m)./(A_c.^m+y(2).^m)).*gamma.*y(1);  
    dy(2) = (alpha_A+alpha_AA.*((y(2).^m)./(A_c.^m+y(2).^m))).*y(1)-gamma_A.*y(2);
    dy(3) = tau * (s-n*y(3)*y(4)-mu*y(3)+(eta1*y(3)*y(4))/(g1+y(4))+(eta2*y(3)*y(1))/(g2+y(1)));
    dy(4) = tau * (r*y(4)*(1-y(4)/k)-e1*y(3)*y(4)-(e2*(beta*((y(2)^m)/(A_c^m+y(2)^m))*gamma*y(1))^2*y(4))/((g3)^2+(beta*((y(2)^m)/(A_c^m+y(2)^m))*gamma*y(1))^2)-e3*y(1)*y(4));
end

function dy = si( t,y )
    alpha = 1; m = 2; A_c = 1; gamma = 4; alpha_A = 0.4; alpha_AA = 8; gamma_A = 1;  
    s = 1.3e4; n = 3.422e-10; mu = 0.0412; eta1 = 0.1245; g1 = 2.019e7;  r = 0.18; k = 5e8; e1 = 1.101e-7; 
    eta2 = 0.2; g2 = 2.5e7; e2 = 0.4; beta = 5; g3 = 2e7; e3 = 0.3; tau = 0.05; 
  
    dy = zeros(4,1);
    dy(1) = alpha.*y(1)-((y(2).^m)./(A_c.^m+y(2).^m)).*gamma.*y(1);
    dy(2) = (alpha_A+alpha_AA.*((y(2).^m)./(A_c.^m+y(2).^m))).*y(1)-gamma_A.*y(2);
    dy(3) = tau * (s-n*y(3)*y(4)-mu*y(3)+(eta1*y(3)*y(4))/(g1+y(4))+(eta2*y(3)*y(1))/(g2+y(1)));
    dy(4) = tau * (r*y(4)*(1-y(4)/k)-e1*y(3)*y(4)-(e2*(beta*((y(2)^m)/(A_c^m+y(2)^m))*gamma*y(1))^2*y(4))/((g3)^2+(beta*((y(2)^m)/(A_c^m+y(2)^m))*gamma*y(1))^2)-e3*y(1)*y(4));
end