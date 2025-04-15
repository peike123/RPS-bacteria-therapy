function bistable_yuanhong()  %画图
alpha = 1; m = 2; A_c = 1; gamma = 20; alpha_A = 0.4; alpha_AA = 8; gamma_A = 1;  
s = 1.3; n = 5e-10; mu = 0.0412; eta1 = 0.1245; g1 = 2.019e7; r = 0.1; k = 5e8; e1 = 1.101e-7; 
eta2 = 0.2; g2 = 2.5e7; e2 = 0.4; beta = 5; g3 = 2e7; e3 = 0.3; tau = 0.05;

Findequilibrium;

ep1 = load('equilibriumpointX.dat');
ep2 = load('equilibriumpointY.dat');

figure(1);
x = [0 2.5e6];
y = [0 7e8];
imagesc(x, y, ep2);
set(gca, 'YDir', 'normal');
alpha(0.75)

grid minor
xlabel({'Effector cells'});
ylabel({'Tumor cells'});

% 添加颜色柱
colormap jet;  % 设置颜色映射方案，可以选择不同的映射方案，如 'parula', 'jet' 等
colorbar;      % 添加颜色柱
caxis([min(ep2(:)) max(ep2(:))]);  % 设置颜色范围，与数据匹配

%画正平衡点
hold on
AA = ((alpha*A_c^m)/(gamma-alpha))^(1/m);
BB = ((gamma*gamma_A)* AA)/(gamma*alpha_A+alpha*alpha_AA);
syms x3 x4;
eqn1=[s-n*x3*x4-mu*x3+(eta1*x3*x4)/(g1+x4)+(eta2*x3*BB)/(g2+BB)==0,...
    r*x4*(1-x4/k)-e1*x3*x4-(e2*(beta*((AA^m)/(A_c^m+AA^m))*gamma*BB)^2*x4)/((g3)^2+(beta*((AA^m)/(A_c^m+AA^m))*gamma*BB)^2)-e3*BB*x4==0]; % 正平衡点满足的方程组,
vars=[x3,x4]; % 定义求解的未知量
[solX3,solX4]=solve(eqn1,vars); % 求解平衡点
plot(solX3,solX4,'k*')
end

function Findequilibrium()  %取点得结果

%构造初值矩阵
xlimt = 0.1:2.5e4:2.5e6;  
ylimt = 0.1:2.5e6:7e8;
ntx = length(xlimt);
nty = length(ylimt);

matrixs1 = zeros(ntx, nty);
matrixs2 = zeros(ntx, nty);

for i = 1:ntx
    for j = 1:nty
        xy = erfenfa(xlimt(i), ylimt(j)); %求解方程
        matrixs1(i, j) = xy(end, 3);  % 将y(3)存放到矩阵中初值的位置
        matrixs2(i, j) = xy(end, 4);  % 将y(4)存放到矩阵中初值的位置
    end
end

dlmwrite('equilibriumpointX.dat', matrixs1');
dlmwrite('equilibriumpointY.dat', matrixs2');

end

function yy = erfenfa(X, Y) %解方程

t0 = 0; tfinal = 1500;
y0 = [0.2 0.5 X Y];  % 初始化 y0, 包括四个初始值
[~, yy] = ode15s(@(t, y) ode(t, y), [t0, tfinal], y0);

end

function dy = ode(~, y)

% 定义参数
alpha = 1; m = 2; A_c = 1; gamma = 20; alpha_A = 0.4; alpha_AA = 8; gamma_A = 1;  
s = 1.3; n = 5e-10; mu = 0.0412; eta1 = 0.1245; g1 = 2.019e7; r = 0.1; k = 5e8; e1 = 1.101e-7; 
eta2 = 0.2; g2 = 2.5e7; e2 = 0.4; beta = 5; g3 = 2e7; e3 = 0.3; tau = 0.05; 

% 初始化导数数组
dy = zeros(4, 1);

% 定义方程组
dy(1) = alpha * y(1) - ((y(2).^m) ./ (A_c.^m + y(2).^m)) .* gamma .* y(1);
dy(2) = (alpha_A + alpha_AA .* ((y(2).^m) ./ (A_c.^m + y(2).^m))) .* y(1) - gamma_A .* y(2);
dy(3) = tau * (s-n*y(3)*y(4)-mu*y(3)+(eta1*y(3)*y(4))/(g1+y(4))+(eta2*y(3)*y(1))/(g2+y(1)));
dy(4) = tau * (r*y(4)*(1-y(4)/k)-e1*y(3)*y(4)-(e2*(beta*((y(2)^m)/(A_c^m+y(2)^m))*gamma*y(1))^2*y(4))/((g3)^2+(beta*((y(2)^m)/(A_c^m+y(2)^m))*gamma*y(1))^2)-e3*y(1)*y(4));

end
