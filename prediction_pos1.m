function [pos1,pos2,pos3,pos4,rsm_pf1,rsm_pf2,rsm_pf3,rsm_pf4,mean_tm,mean_rsm_pf,mean_mean_pf]= prediction_pos1(init_lon,init_lat,sl1,deg1,mag,dataset,fx,xfx,tr, min_lat ,min_lon) 
    num=10;                      %仿真次数
    %rot_r=2.56;
    particle_num=100;            %粒子滤波使用粒子数

    rsm_pf1=zeros(num,1);        %各种方法不同仿真次数的均方根误差
    rsm_pf2=zeros(num,1);
    rsm_pf3=zeros(num,1);
    rsm_pf4=zeros(num,1);

    mean_pf1=zeros(num,1);        %各种方法不同仿真次数的均方根误差
    mean_pf2=zeros(num,1);
    mean_pf3=zeros(num,1);
    mean_pf4=zeros(num,1);

    tm1=zeros(num,1);            %各种方法所用时间
    tm2=zeros(num,1);
    tm3=zeros(num,1);
    tm4=zeros(num,1);
    
    i = 1;
    for i=1:num
        %分别使用下面四种方法得到地磁匹配后的修正位置信息，以及算法运行的时间
       fprintf("****四种算法第%d次测试结果数据展示：\n",i)
        tic 
        t1=clock;
        %传统粒子滤波方法
        pos1(:,:,i)=ParticalFilter1(particle_num,init_lat,init_lon,sl1,deg1,0.1*sl1,0.26817*ones(length(sl1),1),mag,dataset, min_lat ,min_lon);
        %pos1(:,:,i)=ParticalFilter1(particle_num,init_lat,init_lon,sl1,deg1,0*sl1,0*ones(length(sl1),1),mag,dataset);
        t2=clock;
        tm1(i)=etime(t2,t1);
        toc
        
        tic
        % disp(i)
        t1=clock;
        %淼鑫师兄AOFA改进粒子滤波
        pos2(:,:,i)=AOFAParticalFilter1(100,init_lat,init_lon,sl1,deg1,0.1*sl1,0.26817*ones(length(sl1),1),mag,dataset, min_lat ,min_lon);
        t2=clock;
        tm2(i)=etime(t2,t1);
        toc
        
        tic
        t1=clock;
        pos3(:,:,i)=EKPF2(particle_num,init_lat,init_lon,sl1,deg1,0.1*sl1,0.26817*ones(length(sl1),1),mag,dataset, min_lat ,min_lon);
        t2=clock;
        tm3(i)=etime(t2,t1);
        toc
        

        tic
        t1=clock;
        %集成粒子滤波方法
        %pos4(:,:,i) = AIPF4(init_lat,init_lon,sl1,deg1,0.1*sl1,0.26817*ones(length(sl1),1),mag,dataset,fx,xfx);%IPF分布表计算，IPF2积分计算
        pos4(:,:,i) = AIPF4(init_lat,init_lon,sl1,deg1,0.1*sl1,0.26817*ones(length(sl1),1),mag,dataset,fx,xfx, min_lat ,min_lon);
        t2=clock;
        tm4(i)=etime(t2,t1);
        toc
        

   
        %%%%%%%%%各种方法的均方根误差计算%%%%%%%%%%%%%%%%%%%%%%%        
        [rsm1,MEAN1]=CountRSM(pos1(1:end,:,i),tr(1:end,:));
        [rsm2,MEAN2]=CountRSM(pos2(1:end,:,i),tr(1:end,:));
        [rsm3,MEAN3]=CountRSM(pos3(1:end,:,i),tr(1:end,:));
        [rsm4,MEAN4]=CountRSM(pos4(1:end,:,i),tr(1:end,:));

        rsm_pf1(i)=rsm1;
        rsm_pf2(i)=rsm2;
        rsm_pf3(i)=rsm3;
        rsm_pf4(i)=rsm4;

        mean_pf1(i)=MEAN1;
        mean_pf2(i)=MEAN2;
        mean_pf3(i)=MEAN3;
        mean_pf4(i)=MEAN4;
   end
    mean_tm1 = mean(tm1); %平均运行时间计算
    mean_tm2 = mean(tm2);
    mean_tm3 = mean(tm3);
    mean_tm4 = mean(tm4);
    mean_tm = [mean_tm1,mean_tm2,mean_tm3,mean_tm4];

    mean_rsm_pf1 = mean(rsm_pf1);%n次实验平均均方根误差
    mean_rsm_pf2 = mean(rsm_pf2);
    mean_rsm_pf3 = mean(rsm_pf3);
    mean_rsm_pf4 = mean(rsm_pf4);
    mean_rsm_pf = [mean_rsm_pf1,mean_rsm_pf2,mean_rsm_pf3,mean_rsm_pf4];
    
    mean_mean_pf1 = mean(mean_pf1);%n次实验平均均方根误差
    mean_mean_pf2 = mean(mean_pf2);
    mean_mean_pf3 = mean(mean_pf3);
    mean_mean_pf4 = mean(mean_pf4);
    mean_mean_pf = [mean_mean_pf1,mean_mean_pf2,mean_mean_pf3,mean_mean_pf4];
end