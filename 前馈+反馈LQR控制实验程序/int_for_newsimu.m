%% 三自由度直升机控制实验参数设置

clc;
clear;
load('desire_case1.mat');
Sample_Time_Net = 0.005;

mh = 1.8; %直升机本体质量 单位：kg
loc = 0.78; %直升机本体到旋转轴的距离 单位：m

mc = 3.433; %平衡块的质量 单位：kg
lob = 0.33; %平衡快到支点的距离 单位：m

mb = 0.67948;
loa = 0.2296;

m1 = 0.552; %直升机单个推进器的质量 单位：kg
m2 = 0.552; %直升机单个推进器的质量 单位：kg  后面没有用

ldf = 0.17; %直升机推进器到翻转轴的距离 单位：m
lcd = 0.10; %直升机两个推进器连接杠中心到翻转轴的距离 单位：m

I_t = 0.0319;  %直升机俯仰角转动惯量
I_e = 1.469;   %直升机升降角转动惯量
I_psi = 1.469; %直升机偏航角转动惯量
g = 9.81;
m = (-mc*lob + mh *loc + mb*loa)/loc;%有效质量
%m = (-mc*lob + mh *loc)/loc;%有效质量


%% 升降通道辨识数据

a1 = 0.2634;  
a2 = 1.6279;
b1 = 0.531;
ce = 0.3010;

%% 俯仰通道辨识数据

a3 = 25.6488;
b2 = 5.3292;
ct = 0.3410; 

%% 偏航通道辨识数据
b3 = 0.531; 
b4 = 0.1157;
cp = 0.4804;

%% 速度测量二阶模块参数设置
%              w^2 s
%  D(s)=--------------------
%        s^2 + 2wn s + w^2 
% w= 125.66;
w= 100;
n=0.9;
%

%% 定义平衡点
% theta = 0;
theta = 0*pi/180;
theta_d = theta;
epsilon = 0;
epsilon_d = epsilon;
% psi = 45*pi/180;
psi = 0;
psi_d = psi;



%% 计算前馈合力差力

fs = (a1*cos(theta)*sin(epsilon)+a2*cos(epsilon))/(b1*cos(theta));
fd = a3*cos(epsilon)*sin(theta)/b2;


%% 线性化模型
d1 = -b3.*sin(theta).*sin(epsilon).*fs - b4.*sin(theta).*cos(epsilon).*fd;
d2 = b3.*cos(theta).*cos(epsilon).*fs - b4.*cos(theta).*sin(epsilon).*fd;
d3 = -a1.*cos(theta).*cos(epsilon) + a2.*sin(epsilon);
d4 = a1.*sin(theta).*sin(epsilon) - b1.*sin(theta).*fs;
d5 = a3.*sin(theta).*sin(epsilon);
d6 = -a3.*cos(theta).*cos(epsilon);
d7 = b3.*sin(theta).*cos(epsilon);
d8 = -b4.*sin(theta).*sin(epsilon);
d9 = b1.*cos(theta);

% 三通道
% A = [0 1 0 0 0 0;
%      0 -cp d1 0 d2 0;
%      0 0 0 1 0 0;
%      0 0 d3 -ce d4 0;
%      0 0 0 0 0 1;
%      0 0 d5 0 d6 -ct;];
%  B = [0 0;
%      d7 d8;
%      0 0;
%      d9 0;
%      0 0;
%      0 b2;];
 
 % 升降-俯仰双通道模型
 A = [0 1 0 0;
      d3 -ce d4 0;
      0 0 0 1;
      d5 0 d6 -ct;]; 
  
 B = [0 0;
      d9 0;
      0 0;
      0 b2;]; 
 
 
%%  LQR

% Q = [100 0 0 0 0 0;
%      0 1 0 0 0 0;
%      0 0 100 0 0 0;
%      0 0 0 10 0 0;
%      0 0 0 0 400 0;
%      0 0 0 0 0 1];  %%三自由度 偏航无正弦

%三通道
%  Q = [30 0 0 0 0 0;
%      0 01 0 0 0 0;
%      0 0 180 0 0 0;
%      0 0 0 10 0 0;
%      0 0 0 0 400 0;
%      0 0 0 0 0 1];  %偏航+升降 有正弦

% 升降-俯仰双通道
Q = [180 0 0 0;
     0 10 0 0;
     0 0 400 0;
     0 0 0 1]; 
 
R =[1 0;
    0 1];
[K,P,e] = lqr(A,B,Q,R);

K_psi2theta = [1.8,0.2];
