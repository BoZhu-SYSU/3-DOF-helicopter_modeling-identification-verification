close all
load('data0817三通道升降正弦偏航匀速实验.mat')
% size1 = size(elevation_real_data,1);
%% 对比

% 没给期望
% epsilon_d_data = zeros(size1,1);
% for i=1:size1
% epsilon_d_data(i,1) = epsilon_d;
% end
% 
% theta_d_data = zeros(size1,1);
% for i=1:size1
% theta_d_data(i,1) = theta_d;
% end
% 
% psi_d_data = zeros(size1,1);
% for i=1:size1
% psi_d_data(i,1) = psi_d;
% end


%数据中有期望
epsilon_d_data(:,1) = ele_d_data.signals.values(1,1,:);
theta_d_data(:,1) = pitch_d_data.signals.values(1,1,:);
psi_d_data(:,1) = travel_d_data.signals.values(1,1,:);


figure
%仿真
plot(elevation_real_data(:,1),epsilon_d_data(:,1),'b','LineWidth',1);
hold on;
plot(elevation_simu_data(:,1),elevation_simu_data(:,2),'--g','LineWidth',2);
plot(elevation_real_data(:,1),elevation_real_data(:,2),'r','LineWidth',1.5);
%实验
% plot(elevation_real_data(:,1),epsilon_d_data(:,1),elevation_real_data(:,1),elevation_real_data(:,2));
grid on;
%xlim([0,45]); 
xlabel('\fontname{宋体}\fontsize{15}时间\fontname{Times New Roman}\fontsize{15}(s)');
ylabel('\fontname{宋体}\fontsize{15}升降角\fontname{Times New Roman}\fontsize{15}(rad)');
%符号可以识别，公式不行
legend({'期望','仿真','实验'},'Fontname','宋体','FontSize',12);
% legend({'期望','实验'},'Fontname','宋体','FontSize',12);

% % 
figure
%仿真
% plot(pitch_real_data(:,1),theta_d_data(:,1),'b','LineWidth',1);
plot(theta_d(:,1),-theta_d(:,2),'b','LineWidth',1);
hold on;
% plot(pitch_simu_data(:,1),pitch_simu_data(:,2),'--g','LineWidth',2);
plot(pitch_real_data(:,1),-pitch_real_data(:,2),'r','LineWidth',1.5);

%实验
% plot(pitch_real_data(:,1), theta_d_data(:,1),pitch_real_data(:,1),pitch_real_data(:,2));
grid on;
%xlim([0,45]); 
xlabel('\fontname{宋体}\fontsize{15}时间\fontname{Times New Roman}\fontsize{15}(s)');
ylabel('\fontname{宋体}\fontsize{15}俯仰角\fontname{Times New Roman}\fontsize{15}(rad)');
%符号可以识别，公式不行
% legend({'期望','仿真','实验'},'Fontname','宋体','FontSize',12);
legend({'期望','实际'},'Fontname','宋体','FontSize',12);

figure
%仿真
plot(tralvel_real_data(:,1), psi_d_data(:,1),'b','LineWidth',1);
hold on;
plot(travel_simu_data(:,1),travel_simu_data(:,2),'--g','LineWidth',2);
plot(tralvel_real_data(:,1),tralvel_real_data(:,2),'r','LineWidth',1.5);
%实验
% plot(tralvel_real_data(:,1), psi_d_data(:,1),tralvel_real_data(:,1),-tralvel_real_data(:,2));
grid on;
xlabel('\fontname{宋体}\fontsize{15}时间\fontname{Times New Roman}\fontsize{15}(s)');
ylabel('\fontname{宋体}\fontsize{15}偏航角\fontname{Times New Roman}\fontsize{15}(rad)');
legend({'期望','仿真','实验'},'Fontname','宋体','FontSize',12);
% legend({'期望','实际'},'Fontname','宋体','FontSize',12);


% figure
% plot(fs_d(:,1),fs_d(:,2),fd_d(:,1),fd_d(:,2));
% grid on;
% xlabel('\fontname{宋体}\fontsize{15}时间\fontname{Times New Roman}\fontsize{15}(s)');
% ylabel('\fontname{宋体}\fontsize{15}前馈输入\fontname{Times New Roman}\fontsize{15}(N)');
% legend({'合力','差力'},'Fontname','宋体','FontSize',12);
% ylim([-2,6]);

%前馈输入
figure
plot(fs_data(:,1),fs_data(:,2),fd_data(:,1),fd_data(:,2));
grid on;
xlabel('\fontname{宋体}\fontsize{15}时间\fontname{Times New Roman}\fontsize{15}(s)');
ylabel('\fontname{宋体}\fontsize{15}输入\fontname{Times New Roman}\fontsize{15}(N)');
legend({'合力','差力'},'Fontname','宋体','FontSize',12);
ylim([-2,6]);







