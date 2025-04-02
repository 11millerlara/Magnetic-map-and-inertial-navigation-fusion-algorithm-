% 磁地图理论参考仿真器
clear all; clc; close all;

load('./BPIT/fx.mat')
load('./BPIT/xfx.mat')
addpath('./multi_algorithm')
addpath('./cal_pos')
disp("*1*输入初始数据******************************************************")
set(0,'defaultfigurecolor','w')%显示背景设置为白色
%% 读取参考轨迹
load('./magnetic_map_data/predict_data.mat');
magnetic_mod = predict_data(:,:,3)';
min_lon = 0;min_lat = 0;
%真实轨迹矩阵_闭合曲线
[true_pos1,init_pos1_lat,init_pos1_lon,const_i,const_j] = cal_indoor_pos2_10();
%[true_pos1,init_pos1_lat,init_pos1_lon] = cal_33_pos2_3();
number = size(true_pos1,1);
magnetic_pos1=ones(number,1);
for i=1:number
    lat_point = true_pos1(i,1)*10+1;    %读取轨迹坐标的在磁地图中的行数和列数
    lon_point = true_pos1(i,2)*10+1;
    magnetic_pos1(i,1)=magnetic_mod(lat_point,lon_point);
end

%Data1 = Data(1:100,:);
%true_pos1 = Data1(:,1:2);

init_pos1_lat = true_pos1(1,1);
init_pos1_lon = true_pos1(1,2);
%magnetic_pos1 = Data1(:,4);

N = size(magnetic_mod,1);
M = size(magnetic_mod,2);
dataset = zeros(N,M,1);
dataset(:,:,1) = magnetic_mod;


%% 生成轨迹数据近似的步长数据矩阵和方向数据矩阵——曲折闭合曲线
M =size(true_pos1,1);
sl1=zeros(M-1,1);
deg1=zeros(M-1,1);
for i=1:M-1
    sl1(i)=sqrt((true_pos1(i+1,1)-true_pos1(i,1)).^2+(true_pos1(i+1,2)-true_pos1(i,2)).^2);
    y=true_pos1(i+1,1)-true_pos1(i,1);
    x=true_pos1(i+1,2)-true_pos1(i,2);
    theta = atan2d(y, x); % 两点之间的方向角，单位为度
    deg1(i) = theta; % 将方向角（度）存储到 deg1 中
end
% 生成与原始矩阵相同大小的随机矩阵作为干扰项
random_noise_sl1 = rand(size(sl1));
random_noise_deg1 = rand(size(deg1));
% 设置干扰系数
interference_coefficient_sl1 =0.01;
interference_coefficient_deg1 =0.1;
% 计算干扰后的矩阵
sl1 = sl1 + interference_coefficient_sl1 * random_noise_sl1;
deg1 = deg1 + interference_coefficient_deg1 * random_noise_deg1;


% 初始化路径矩阵，预先分配空间
n_steps = length(sl1);
ins_pos = zeros(n_steps+1, 2); % 第一行为初始位置，后面为每一步的位置
ins_pos(1, :) = [init_pos1_lat, init_pos1_lon]; % 存储初始位置

% 将角度转换为弧度
radians = deg1*pi/180;
%deg1 = radians;
% 计算每一步的轨迹位置
for i = 1:n_steps
    % 当前步长和方向
    step_length = sl1(i);
    angle = radians(i);
    
    % 计算新位置
    delta_y = step_length * sin(angle); % y方向增量
    delta_x = step_length * cos(angle); % x方向增量
    
    % 更新位置
    new_y = ins_pos(i, 1) + delta_y;
    new_x = ins_pos(i, 2) + delta_x;
    
    % 存储新位置
    ins_pos(i+1, :) = [new_y, new_x];
end



deg1 = deg1*pi/180;
%% 利用粒子滤波算法计算地磁匹配结果
number_mag = 1; 
disp("*6*分别利用4种算法计算地磁匹配结果数据***************************************")
disp("*6.1*计算第一个轨迹的地磁匹配结果数据")
[pos_PF_1,pos_AOFA_1,pos_EKF_1,pos_AIPF_1,rsm_pf1_1,rsm_pf2_1,rsm_pf3_1,rsm_pf4_1,mean_tm_1,mean_rsm_pf_1,mean_mean_pf_1]=...
    prediction_pos1(init_pos1_lon,init_pos1_lat,sl1,deg1,magnetic_pos1(:,[1:number_mag]),dataset(:,:,[1:number_mag]),fx,xfx,true_pos1,min_lat,min_lon);

%        轨迹1             
stdpf1_1 = std(rsm_pf1_1);  
stdpf2_1 = std(rsm_pf2_1);  
stdpf3_1 = std(rsm_pf3_1);  
stdpf4_1 = std(rsm_pf4_1); 

[min_value_1_1, index_1_1] = min(rsm_pf1_1); % 找到最小值及其索引
pos_PF_1 = pos_PF_1(:,:,index_1_1);


[min_value_2_1, index_2_1] = min(rsm_pf2_1); % 找到最小值及其索引
pos_AOFA_1 = pos_AOFA_1(:,:,index_2_1);


[min_value_3_1, index_3_1] = min(rsm_pf3_1); % 找到最小值及其索引
pos_EKF_1 = pos_EKF_1(:,:,index_3_1);


[min_value_4_1, index_4_1] = min(rsm_pf4_1); % 找到最小值及其索引
pos_AIPF_1 = pos_AIPF_1(:,:,index_4_1);



%平均匹配误差指的是，n次实验的均方根误差的均值；稳定性指n次实验的标准差
disp("*7*分析：****************************************************************")
disp("                        PF       AOFA    EKF      AIPF            ")
fprintf("轨迹1：平均计算时间 ： %f       %f      %f       %f\n",mean_tm_1(1),mean_tm_1(2),mean_tm_1(3),mean_tm_1(4));

disp("                                                              ")
fprintf("轨迹1：平均匹配误差 ： %f       %f      %f       %f\n",mean_rsm_pf_1(1),mean_rsm_pf_1(2),mean_rsm_pf_1(3),mean_rsm_pf_1(4));

disp("                                                              ")
%fprintf("轨迹1：稳定性 ：      %f       %f      %f       %f\n",stdpf1_1,stdpf2_1,stdpf3_1,stdpf4_1);
%fprintf("轨迹2：稳定性 ：      %f       %f      %f       %f\n",stdpf1_2,stdpf2_2,stdpf3_2,stdpf4_2);
%fprintf("轨迹3：稳定性：       %f       %f      %f       %f\n",stdpf1_3,stdpf2_3,stdpf3_3,stdpf4_3);

fprintf("轨迹1：平均误差均值： %f       %f      %f       %f\n",mean_mean_pf_1(1),mean_mean_pf_1(2),mean_mean_pf_1(3),mean_mean_pf_1(4));


%fprintf("轨迹2：阱的大小： %f \n",rmse_matrix(const_j*10+1,const_i*10+1));
%绘图展示四种算法num次计算得到的轨迹与真实轨迹的误差值图

plotResult(rsm_pf1_1,rsm_pf2_1,rsm_pf3_1,rsm_pf4_1);
hold off;

N = size(true_pos1,1);
Data1 = zeros(N,4);
Data1(:,1:2) = true_pos1;
Data1(:,4) = magnetic_pos1;
%% 创建磁地图的三维表面
X0 = predict_data(:,:,1)';
Y0 = predict_data(:,:,2)';
Z0 = predict_data(:,:,3)';
figure(7);
mesh(Y0 ,X0 ,Z0)
hold on
% 绘制真实轨迹路径
h1=plot3(Data1(:,2) ,Data1(:,1) ,Data1(:,4),'k*-','MarkerFaceColor','k');
% 绘制AIPF轨迹路径
h2=plot3(pos_AIPF_1(:,2) , pos_AIPF_1(:,1) , Data1(:,4),'r-o', 'LineWidth', 0.2);
% 绘制PF轨迹路径
h3=plot3( pos_PF_1(:,2) , pos_PF_1(:,1) ,Data1(:,4), 'g-*', 'LineWidth', 0.2); 
% 绘制AOFA轨迹路径
h4=plot3(pos_AOFA_1(:,2) , pos_AOFA_1(:,1) , Data1(:,4),'b-s', 'LineWidth', 0.2); 
% 绘制EKF轨迹路径
h5=plot3(pos_EKF_1(:,2) , pos_EKF_1(:,1) , Data1(:,4),'m-d', 'LineWidth', 0.2); 

% 设置图例
legend([h1,h2,h3,h4,h5],'True', 'AIPF', 'PF', 'AOFA', 'EKPF');

view(3)
xlabel('Normalized longitude','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
ylabel('Normalized latitude','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
zlabel('Magnetic anomaly(nT)','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
colorbar

figure(8);

hold on

% 绘制真实轨迹路径
h1=plot(Data1(:,2) ,Data1(:,1) ,'k*-','MarkerFaceColor','k');
% 绘制AIPF轨迹路径
h2=plot(pos_AIPF_1(:,2) , pos_AIPF_1(:,1) ,'r-o', 'LineWidth', 0.2);
% 绘制PF轨迹路径
h3=plot( pos_PF_1(:,2) , pos_PF_1(:,1) , 'g-*', 'LineWidth', 0.2); 
% 绘制AOFA轨迹路径
h4=plot(pos_AOFA_1(:,2) , pos_AOFA_1(:,1) ,'b-s', 'LineWidth', 0.2); 
% 绘制EKF轨迹路径
h5=plot(pos_EKF_1(:,2) , pos_EKF_1(:,1) ,'m-d', 'LineWidth', 0.2); 

% 映射磁地图，并设置透明度
surf(Y0, X0 , Z0, magnetic_mod, 'EdgeColor', 'none', 'FaceAlpha', 0.6); % 将透明度设为0.6
%mesh(Y0,X0,Z0)
% 设置图例
legend([h1,h2,h3,h4,h5],'True', 'AIPF', 'PF', 'AOFA', 'EKPF');


xlabel('Normalized longitude','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
ylabel('Normalized latitude','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
zlabel('Magnetic anomaly(nT)','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');

% 计算每个算法估计位置与真实位置之间的误差
error_PF = sqrt(sum((pos_PF_1  - true_pos1 ).^2, 2));
error_AOFA = sqrt(sum((pos_AOFA_1  - true_pos1 ).^2, 2));
error_EKF = sqrt(sum((pos_EKF_1  - true_pos1 ).^2, 2));
error_AIPF = sqrt(sum((pos_AIPF_1  - true_pos1 ).^2, 2));

% 绘制CDF图
figure(9);
hold on;

% 绘制各算法误差的CDF
cdfplot(error_PF);
cdfplot(error_AOFA);
cdfplot(error_EKF);
cdfplot(error_AIPF);

% 设置图例
legend('PF', 'AOFA', 'EKPF','AIPF', 'Location', 'best');

% 设置轴标签和标题
xlabel('Error','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
ylabel('Cumulative Distribution Function (CDF)','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');
title(' ','FontSize',12,'Fontname', 'Times New Roman','FontWeight','bold');

% 设置网格
grid on;

hold off;







