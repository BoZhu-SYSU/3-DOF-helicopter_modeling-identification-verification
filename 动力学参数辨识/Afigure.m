close all
load('experiment_data.mat');
load('ModelInputSin.mat');
%% 建模实验数据和仿真数据对比
%输入信号
figure
plot(controlInput2(:,1),controlInput2(:,2),controlInput2(:,1),controlInput2(:,3),'LineWidth',1.5)
xlabel('\fontname{宋体}\fontsize{18}时间\fontname{Times New Roman}\fontsize{18}(s)');
ylabel('\fontname{宋体}\fontsize{18}升力\fontname{Times New Roman}\fontsize{18}(N)');
h = legend({'f_1','f_2'}, 'Fontname','Times New Roman','FontSize',18);
title('\fontname{Times New Roman}\fontsize{18}a.\fontname{宋体}\fontsize{18}系统输入');
set(h,'Orientation','horizon');
set(h,'Box','off');
set(gca,'FontSize',15,'Fontname', 'Times New Roman');
set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.4]); 

figure
plot(exqforInput1(:,1),exqforInput1(:,4),'b',exqforInput1(:,1),SingleModelq(:,2),'g',exqforInput1(:,1),MultiModelq_undenti(:,2),'r', 'LineWidth',1);
% grid on;
%xlim([0,45]); 
xlabel('\fontname{宋体}\fontsize{18}时间\fontname{Times New Roman}\fontsize{18}(s)');
ylabel('\fontname{宋体}\fontsize{18}偏航角\fontname{Times New Roman}\fontsize{18}(rad)');
% h = legend({'实验','单刚体模型','多刚体模型'},'Fontname','宋体','FontSize',18);
h = legend({'实验','SRB模型','MRB模型'},'Fontname','宋体','FontSize',18);

set(h,'Orientation','horizon');
set(h,'Box','off');
title('\fontname{Times New Roman}\fontsize{18}b.\fontname{宋体}\fontsize{18}偏航角');
set(gca,'FontSize',15,'Fontname', 'Times New Roman');
set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.4]);

figure
plot(exqforInput1(:,1),exqforInput1(:,2),'b',exqforInput1(:,1),SingleModelq(:,3),'g',exqforInput1(:,1),MultiModelq_undenti(:,3),'r', 'LineWidth',1);
% grid on;
%xlim([0,45]); 
xlabel('\fontname{宋体}\fontsize{18}时间\fontname{Times New Roman}\fontsize{18}(s)');
ylabel('\fontname{宋体}\fontsize{18}升降角\fontname{Times New Roman}\fontsize{18}(rad)');
% h = legend({'实验','单刚体模型','多刚体模型'},'Fontname','宋体','FontSize',18);
h = legend({'实验','SRB模型','MRB模型'},'Fontname','宋体','FontSize',18);
set(h,'Orientation','horizon');
set(h,'Box','off');
title('\fontname{Times New Roman}\fontsize{18}c.\fontname{宋体}\fontsize{18}升降角');
set(gca,'FontSize',15,'Fontname', 'Times New Roman');
set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.4]);

figure
plot(exqforInput1(:,1),-exqforInput1(:,3),'b',exqforInput1(:,1),SingleModelq(:,4),'g',exqforInput1(:,1),MultiModelq_undenti(:,4),'r', 'LineWidth',1);
% grid on;
hold on;
xlabel('\fontname{宋体}\fontsize{18}时间\fontname{Times New Roman}\fontsize{18}(s)');
ylabel('\fontname{宋体}\fontsize{18}俯仰角\fontname{Times New Roman}\fontsize{18}(rad)');
% h = legend({'实验','单刚体模型','多刚体模型'},'Fontname','宋体','FontSize',18);
h = legend({'实验','SRB模型','MRB模型'},'Fontname','宋体','FontSize',18);
set(h,'Orientation','horizon');
set(h,'Box','off');
title('\fontname{Times New Roman}\fontsize{18}d.\fontname{宋体}\fontsize{18}俯仰角');
set(gca,'FontSize',15,'Fontname', 'Times New Roman');
set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.4]);



%% 模型辨识画图
figure
plot(controlInput2(:,1),controlInput2(:,2),controlInput2(:,1),controlInput2(:,3))   % 力
%xlim([0,25])
xlabel('\fontname{宋体}\fontsize{15}时间\fontname{Times New Roman}\fontsize{15}(s)');
ylabel('\fontname{宋体}\fontsize{15}升力\fontname{Times New Roman}\fontsize{15}(N)');
h = legend({'$f_1$','$f_2$'}, 'Interpreter','latex','FontSize',18);
set(h,'Orientation','horizon');
set(h,'Box','off');
set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.4]);

% 对比

figure
plot(exqforInput1(:,1),exqforInput1(:,4),'b',exqforInput1(:,1),MultiModelq_undenti(:,2),'g',exqforInput1(:,1),MultiModelq(:,2),'r', 'LineWidth',1.5);
% grid on;
%xlim([0,45]); 
xlabel('\fontname{宋体}\fontsize{18}时间\fontname{Times New Roman}\fontsize{18}(s)');
ylabel('\fontname{宋体}\fontsize{18}偏航角\fontname{Times New Roman}\fontsize{18}(rad)');
%符号可以识别，公式不行
h = legend({'实验','辨识前','辨识后'},'Fontname','宋体','FontSize',18);
set(h,'Orientation','horizon');
set(h,'Box','off');
title('\fontname{Times New Roman}\fontsize{18}a.\fontname{宋体}\fontsize{18}偏航角');
set(gca,'FontSize',15,'Fontname', 'Times New Roman');
set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.4]);
% % 
figure
plot(exqforInput1(:,1), exqforInput1(:,2),'b',exqforInput1(:,1),MultiModelq_undenti(:,3),'g',exqforInput1(:,1),MultiModelq(:,3),'r', 'LineWidth',1.5);
% grid on;
%xlim([0,45]); 
xlabel('\fontname{宋体}\fontsize{18}时间\fontname{Times New Roman}\fontsize{18}(s)');
ylabel('\fontname{宋体}\fontsize{18}升降角\fontname{Times New Roman}\fontsize{18}(rad)');
%符号可以识别，公式不行
h = legend({'实验','辨识前','辨识后'},'Fontname','宋体','FontSize',18);
set(h,'Orientation','horizon');
set(h,'Box','off');
title('\fontname{Times New Roman}\fontsize{18}b.\fontname{宋体}\fontsize{18}升降角');
set(gca,'FontSize',15,'Fontname', 'Times New Roman');
set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.4]);
% % 
figure
plot(exqforInput1(:,1),-exqforInput1(:,3),'b',exqforInput1(:,1),MultiModelq_undenti(:,4),'g',exqforInput1(:,1),MultiModelq(:,4),'r', 'LineWidth',1.5);
% grid on;
hold on;
xlabel('\fontname{宋体}\fontsize{18}时间\fontname{Times New Roman}\fontsize{18}(s)');
ylabel('\fontname{宋体}\fontsize{18}俯仰角\fontname{Times New Roman}\fontsize{18}(rad)');
h = legend({'实验','辨识前','辨识后'},'Fontname','宋体','FontSize',18);
%set(h,'Orientation','horizon');
set(h,'Box','off');
title('\fontname{Times New Roman}\fontsize{18}c.\fontname{宋体}\fontsize{18}俯仰角');
set(gca,'FontSize',15,'Fontname', 'Times New Roman');
set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.4]);
axes('position',[0.18,0.3,0.3,0.25]); %第二组
% axes('position',[0.2,0.56,0.3,0.25]); %第一组
plot(exqforInput1(:,1),-exqforInput1(:,3),'b',exqforInput1(:,1),MultiModelq(:,4),'r', 'LineWidth',1);
% xlim([40,50]);  %第三组0.9863
xlim([25,30]); %第二组
% xlim([30,35]); %第一组
