clear;
clc;
Sample_Time_Net = 0.005;
load('Ai17.mat');

% 数据列信息 time ele pitch travel
Ts = 0.005;
l = 2;
     
e = -exqforInput1(:,l)'; 
len = length(e);

theta = -exqforInput1(1:len,3)';  %theta序列
psi = -exqforInput1(1:len,4)';    %psi序列

ele= -exqforInput1(1:len,2)';     %ele序列
len = length(ele);
figure(1)
plot(0:Ts:(len-1)*Ts,ele);
 xlim([0,45]);
xlabel('\fontname{宋体}\fontsize{15}时间\fontname{Times New Roman}\fontsize{15}(s)');
ylabel('\fontname{宋体}\fontsize{15}升降角\fontname{Times New Roman}\fontsize{15}(rad/s)');
le = legend('${\epsilon}$');
set(le,'Interpreter','latex','FontSize',14);

dot_ex = -ex_dotqforInput(:,2);
dot_e = smooth(dot_ex, 60)'; 
figure(2)
plot(0:Ts:(len-1)*Ts,dot_ex);
hold on;
plot(0:Ts:(len-1)*Ts,dot_e,'color',[0 1 0]);

 xlim([0,45]);
xlabel('\fontname{宋体}\fontsize{15}时间\fontname{Times New Roman}\fontsize{15}(s)');
ylabel('\fontname{宋体}\fontsize{15}升降角速度\fontname{Times New Roman}\fontsize{15}(rad/s)');
le = legend('$\dot{\epsilon}$','$\dot{\bar{\epsilon}}$');
set(le,'Interpreter','latex','FontSize',14);
hold off;

% ddot_ex= -ex_ddotqforInput(:,2);            %FTC估计角加速度

for i = 1:(len-1)                                  %差分算角加速度
    ddot_ex(i)= (dot_e(i+1) - dot_e(i))/Ts;
end
ddot_ex(len)= ddot_ex(len-1);


ddot_e = smooth(ddot_ex, 60)';
figure(3)
plot(0:Ts:(len-1)*Ts,ddot_ex);
hold on;
plot(0:Ts:(len-1)*Ts,ddot_e,'color',[0 1 0]);
 xlim([0,45]);
xlabel('\fontname{宋体}\fontsize{15}时间\fontname{Times New Roman}\fontsize{15}(s)');
ylabel('\fontname{宋体}\fontsize{15}升降角加速度\fontname{Times New Roman}\fontsize{15}(rad/s^2)');
le = legend('$\ddot{\epsilon}$','$\ddot{\bar{\epsilon}}$');
set(le,'Interpreter','latex','FontSize',14);
hold off;
fs = 2.8*ones(1,len);


