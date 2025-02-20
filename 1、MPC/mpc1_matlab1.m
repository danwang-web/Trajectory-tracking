
%只考虑3个状态量的MPC
%% MPC系列matlab.m文件写轨迹跟踪代码
%清屏函数
clc;
clear;
close all;
%% 控制器参数设计
Nx=3;%状态量个数
Nu =2;%控制量个数
Np =20;%预测步长
Nc=10;%控制步长
Row=10;%松弛因子
X0=[0 0 0];%初始位置x，y，yaw都是0
%[Nr,Nc]=size(Xout);%
C=[1 0 0;0 1 0;0 0 1];%输出三个状态量
L=2.6;%轴距2.6
vd1=10;%车辆参数速度1m/s
vd2=0;%车辆参考前轮转角0rad
%% 参考轨迹生成
N=201;%仿真步长
T=0.1;%采样时间
load road_ref.mat;%导入道路数据
Xref=zeros(N,3);
Tout=zeros(N,1);
Xref(:,1)=out.xr;%参考X
Xref(:,2)=out.yr;%参考Y
Xref(:,3)=out.yawr;%参考Yaw
Tout(:,1)=T*out.tout;%仿真时间

% Xref=zeros(N,3);%定义参考轨迹状态量矩阵
% Tout=zeros(N,1);%定义时间矩阵
% for k=1:1:N
%     Xref(k,1)=k*vd1*T;
%     Xref(k,2)=2;
%     Xref(k,3)=0;
%     Tout(k,1)=(k-1)*T;
%  end
%% MPC主体
x_real=zeros(N,Nx);%车辆的实际状态量N行，Nx列
x_piao=zeros(N,Nx);%车辆的状态量误差N行，Nx列
u_real=zeros(N,Nu);%车辆的实际控制量N行，Nc列
u_piao=zeros(N,Nu);%车辆的控制量误差N行，Nc列
x_real(1,:)=X0;%将车辆的初始位置放入车辆实际状态量矩阵的第一行
x_piao(1,:)=X0-x_real(1,:);%对状态误差的第一行进行赋值
X_PIAO=zeros(N,Nx*Np);%储存每一个时刻预测时域内状态量误差
XXX=zeros(N,Nx*Np);%用于保持每个时刻预测的所有状态值
% kesi=[x_piao,u_piao];%定义误差矩阵状态量误差和控制量误差都有
Q=eye(Np*Nx);%状态量权重矩阵Q
R=eye(Nc*Nu);%控制量权重矩阵R
for i=1:1:N
    yawref=Xref(i,3);%将参考轨迹的每一时刻的第三个信号赋值为参考横摆角
    A1=[1    0   -vd1*sin(yawref)*T;
        0    1   vd1*cos(yawref)*T;
        0    0   1];%运动学误差模型状态矩阵
    B1=[cos(yawref)*T   0;
        sin(yawref)*T   0;
        tan(vd2)*T/L vd1*T/L/(cos(vd2))^2];%运动学误差模型控制矩阵
    %%公式4-10
    %更新矩阵
    PHI_cell=cell(Np,1);
    THETA_cell=cell(Np,Nc);
    for k=1:1:Np
        PHI_cell{k,1}=C*A1^k;
        for j=1:1:Nc
            if  j<=k
                THETA_cell{k,j}=C*A1^(k-j)*B1;
            else
                THETA_cell{k,j}=zeros(Nx,Nu);
            end
        end
    end
    PHI=cell2mat(PHI_cell);
    THETA=cell2mat(THETA_cell);
    %%公式4-19
    H=THETA'*Q*THETA+R;
    f=2*THETA'*Q*PHI*x_piao(i,:)';
    A_cons=[];
    B_cons=[];
    lb=[-0.2;-0.44];
    ub=[0.2;0.44];
    tic
    options=optimset('Algorithm','interior-point-convex');
    [X,fval(i,1),exitflag(i,1),output(i,1)]=quadprog(H,f,A_cons,B_cons,[],[],lb,ub);
    toc
    X_PIAO(i,:)=(PHI*x_piao(i,:)'+THETA*X(1:Np,1))';%保存每一计算步长预测时域的未来状态变化量
    if i<N
         for j=1:1:Np
             XXX(i,1+3*(j-1))=X_PIAO(i,1+3*(j-1))+Xref(i,1);
             XXX(i,2+3*(j-1))=X_PIAO(i,2+3*(j-1))+Xref(i,2);
             XXX(i,3+3*(j-1))=X_PIAO(i,3+3*(j-1))+Xref(i,3);
         end
    else
         for j=1:1:Np
             XXX(i,1+3*(j-1))=X_PIAO(i,1+3*(j-1))+Xref(N,1);
             XXX(i,2+3*(j-1))=X_PIAO(i,2+3*(j-1))+Xref(N,2);
             XXX(i,3+3*(j-1))=X_PIAO(i,3+3*(j-1))+Xref(N,3);
         end
    end
    u_piao(i,1)=X(1,1);%每一步速度控制量增量
    u_piao(i,2)=X(2,1);%每一步转角控制量增量

    X00=x_real(i,:);
    vd11=vd1+u_piao(i,1);%实际每一步的速度控制量
    vd22=vd2+u_piao(i,2);%实际每一步的转角控制量
    XOUT=dsolve('Dx-vd11*cos(z)=0','Dy-vd11*sin(z)=0','Dz-tan(vd22)/L*vd11=0','x(0)=X00(1)','y(0)=X00(2)','z(0)=X00(3)');
     t=T;%求解器采样时间
     x_real(i+1,1)=eval(XOUT.x);%求解微分方程纵向位移X
     x_real(i+1,2)=eval(XOUT.y);%求解微分方程横向位移Y
     x_real(i+1,3)=eval(XOUT.z);%求解微分方程航向角Yaw
     %计算每一步的状态量偏差
     if(i<N)
         x_piao(i+1,:)=x_real(i+1,:)-Xref(i+1,:);
     end
    u_real(i,1)=vd1+u_piao(i,1);
    u_real(i,2)=vd2+u_piao(i,2);
    
    figure(1);
    %轨迹图
    hold on;
    
    title('跟踪结果对比');
    xlabel('横向位置X');
    axis([-1 vd1*N*T -1 4]);
    ylabel('纵向位置Y');
    hold on;
    for k=1:1:Np
         X1(i,k+1)=XXX(i,1+3*(k-1));%储存实际每一步的速度控制量
         Y1(i,k+1)=XXX(i,2+3*(k-1));%储存实际每一步的转角控制量
         plot(X1(i,:),Y1(i,:),'y');
         hold on;
    end
    X(i,1)=x_real(i,1);
    Y(i,1)=x_real(i,2);
    plot(X(i,:),Y(i,:),'r','LineWidth',1,'Marker','*');
    hold on;
    plot(Xref(i,1),Xref(i,2),'b-','LineWidth',1,'Marker','o');
    title('跟踪结果对比');
    xlabel('横向位置X');
    ylabel('纵向位置Y');


    figure(2)
    %偏差图
    subplot(3,1,1); 
    plot(Tout(i),x_piao(i,1),'k--','Marker','o');
    hold on;
    %grid on;
    title('纵向误差');
    xlabel('采样时间T');
    ylabel('纵向误差ΔX(m)')
    subplot(3,1,2);
    plot(Tout(i),x_piao(i,2),'k--','Marker','+');
    hold on;
    %grid on;
    title('横向误差');
    xlabel('采样时间T');
    ylabel('横向误差ΔY(m)')
    subplot(3,1,3);
    plot(Tout(i),x_piao(i,3),'k--','Marker','x');
    hold on;
    %grid on;
    title('航向角误差');
    xlabel('采样时间T');
    ylabel('航向误差ΔYaw(rad)')


    figure(3)
    %状态量
    subplot(3,1,1);
    plot(Tout(i),u_real(i,2),'k--','Marker','o');
    hold on;
    %grid on;
    title('前轮转角');
    xlabel('采样时间T');
    ylabel('前轮转角(rad)')
    subplot(3,1,2);
    plot(Tout(i),x_real(i,2),'k--','Marker','o');
    hold on;
    %grid on;
    title('横向位移');
    xlabel('采样时间T');
    ylabel('横向位移Y(m)')
    subplot(3,1,3);
    plot(Tout(i),x_real(i,3),'k--','Marker','o');
    hold on;
    %grid on;
    title('航向角');
    xlabel('采样时间T');
    ylabel('航向角Yaw(rad)')


    figure(4)
    %状态量的微分
    subplot(3,1,1);
    x_dt(i+1,1)=(x_real(i+1,1)-x_real(i,1))/T;
    y_dt(i+1,1)=(x_real(i+1,2)-x_real(i,2))/T;
    yaw_dt(i+1,1)=(x_real(i+1,3)-x_real(i,3))/T;
    plot(Tout(i),x_dt(i+1,1),'k--','Marker','o');
    hold on;
    %grid on;
    title('纵向速度');
    xlabel('采样时间T');
    ylabel('纵向速度vx(m/s)')
    subplot(3,1,2);
    plot(Tout(i),y_dt(i+1,1),'k--','Marker','o');
    hold on;
    %grid on;
    title('横向速度');
    xlabel('采样时间T');
    ylabel('横向速度vy(m/s)')
    subplot(3,1,3);
    plot(Tout(i),yaw_dt(i+1,1),'k--','Marker','o');
    hold on;
    %grid on;
    title('横摆角速度');
    xlabel('采样时间T');
    ylabel('横摆角速度wr(rad/s)')
end
