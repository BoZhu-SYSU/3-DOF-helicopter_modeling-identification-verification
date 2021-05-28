
clear;
load('ModelInputSin.mat');
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

I_t = 0.0319;
I_e = 1.469;
I_psi = 1.469;
g = 9.81;
m = (-mc*lob + mh *loc + mb*loa)/loc;%有效质量
%m = (-mc*lob + mh *loc)/loc;%有效质量

q0 = [0 0 0]';
dot_q0 = [0 0 0]';

Ts = 0.005;  %每步采样步长
Time = 50;

%% 升降有输入的情况 
a1 = 0.2634;  %%
a2 = 1.6279;
b1 = 0.531;
ce = 0.10; 

%% 俯仰

a3 = 25.6488;
b2 =  5.3292;
ct = 0.3410;  

%% 偏航A12
b3 =  0.531;  
b4 =  0.1157;   
cp = 0.4804; 

cp1 = 0.274;
ce1 = 0.053;
ct1 = 0.048;
