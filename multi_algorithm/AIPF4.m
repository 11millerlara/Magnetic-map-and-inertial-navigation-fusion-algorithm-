function pos = IPF4(init_x,init_y,sl1,deg1,sigma1,sigma2,mag1,dataset,fx,xfx, min_lat ,min_lon)
%%%%%%%%初始化操作%%%%%%%%%%%%%%%%%%%%%%
v_l=0;                  %初始速度
v_d=0;
sizedataset = size(mag1,2);
step=length(sl1);       %路径的总步数，即采集的数据点数目
pos=zeros(step+1,2);    %用于存储估计的位置信息
pos(1,1)=init_x;        %第一行表示起始位置（init_x,init_y）
pos(1,2)=init_y;
mmp=0.1;                %地磁场地图精度，设置为0.1
n=1;                  
pp=mmp/n;               %特征长度
%P=diag([3^2]);
P=3^2*eye(sizedataset); %P是一个协方差矩阵，初始化对角线上的元素

%%%%%%%%%迭代处理每一步的位置和方向数据，同时更新速度估计，最终得到存储每一步位置信息的pos矩阵%%%%%%%%%%%%
    for s=1:step
        sl=sl1(s);          %从s11数组中获取当前步长数据
        deg=deg1(s);        %从deg1数组中获取当前方向数据
        mag=mag1(s,:);      %从mag1矩阵中获取当前的地磁信息
           
        v_l=sqrt(v_l+(max(0.05,0.1*sl))^2);%更新线性速度估计v_1
        v_d=0.26817;        %设定v_d方向速度估计为常数值
         
        x=pos(s,1);         %从pos数组中获取当前步骤的位置坐标
        y=pos(s,2);
        
        sigma1=v_l;         %将线性速度估计和方向速度估计存储到sigma1和sigma2变量中
        sigma2=v_d;
        %row为航向边划分区域的个数；col为步长边划分区域的个数
        row=ceil((sl+3*sigma1)*6*sigma2/pp);
        col=ceil(6*sigma1/pp);
            
        row_p=6*sigma2/row; % row_p = d/(l0 + 3sigma1);
        col_p=6*sigma1/col; % col_p = d;
        
        pn=row*col;         %计算网格中的总区域数
        
        %分别创建用于存储概率p，步长估计sl_e，方向估计deg_e的数组，它们的长度为pn
        p=zeros(pn,1);
        sl_e=zeros(pn,1);
        deg_e=zeros(pn,1);
        for m=1:10 
                for i=1:row
                    for j=1:col
                        %分别计算当前网格区域的行边界row_l,row_s和列边界col_1,col_s
                        row_l=-3+(i-1)*row_p/sigma2;
                        row_s=-3+i*row_p/sigma2;
                        col_l=-3+(j-1)*col_p/sigma1;
                        col_s=-3+j*col_p/sigma1;
                        
                        num=(i-1)*col+j;%计算当前网格单元的索引，用于数组中存储相关参数
                        
                        %计算概率p，计算步长估计sl_e，计算方向估计deg_e
                        p(num)=Integrate(row_l,row_s,fx)*Integrate(col_l,col_s,fx);
                        sl_e(num)=(sigma1*fxIntegrate(col_l,col_s,xfx)+sl*Integrate(col_l,col_s,fx))*Integrate(row_l,row_s,fx)/p(num);
                        deg_e(num)=(sigma2*fxIntegrate(row_l,row_s,xfx)+deg*Integrate(row_l,row_s,fx))*Integrate(col_l,col_s,fx)/p(num);
                        
                        %这两行代码基于当前位置（x,y），步长估计和方向估计计算得到位置估计的更新
                        x_t=x+sl_e(num)*sin(deg_e(num));
                        y_t=y+sl_e(num)*cos(deg_e(num));
                        
                        index=Find3(dataset,x_t,y_t, min_lat ,min_lon);%根据更新后的位置坐标在dataset中查找相关数据点的索引
                        
                        if(index(1)==-1||index(2)==-1)%没有找到
                            p(num)=0;                            %将p设置为0，表示在该位置的概率为0
                        else
                            
                            %计算临近数据点与观测数据 (mag) 之间的差异，并基于协方差矩阵 (P) 计算了概率。######
                            %tmp=(exp(-0.5*(mag-reshape(dataset(index(1),index(2),4:6),1,3))*inv(P)*(mag-reshape(dataset(index(1),index(2),4:6),1,3))'));
                            tmp=(exp(-0.5*(mag-reshape(dataset(index(1),index(2),1:sizedataset),1,sizedataset))*inv(P)*(mag-reshape(dataset(index(1),index(2),1:sizedataset),1,sizedataset))'));
                            p(num)=p(num)*tmp; %基于数据点与观测数据mag之间的关系tmp更新概率p
                           
                            if(isnan(p(num)))%检查概率p是否是NaN（无效）,如果是NaN，则将其设置为0
                                p(num)=0;      
                            end                
                        end             
                    end
                end
        
                sum_p=sum(p,1);     %计算概率向量 p 中所有元素的总和，得到 sum_p。
                p=p/sum_p;          %将每个概率值除以总和来将概率向量 p 归一化，确保它们的总和为1。
                
                %创建两个零向量，sl_n 和 deg_n，用于存储新的步长和方向估计。
                sl_n=zeros(pn,1);   
                deg_n=zeros(pn,1);
               
                %创建向量 c，并将其第一个元素设置为 p(1)。 c 向量用于计算累积概率。
                c=zeros(pn,1);
                c(1)=p(1);
                
                %计算 c 向量的累积概率。
                for j=2:pn
                    c(j)=c(j-1)+p(j);
                end
            
                %用于生成新的步长和方向估计。
                for j=1:pn
                    a=unifrnd(0,1);             %生成一个0到1之间的随机数 a。
            
                    for k=1:pn                  %这个循环用于查找随机数 a 在累积概率向量 c 中对应的位置 k，以确定新的步长和方向估计。
                        %检查随机数 a 是否小于累积概率 c(k)，如果满足条件，则将该位置 k 的步长和方向估计分配给 sl_n(j) 和 deg_n(j)。
                        if(a<c(k))                 
                            sl_n(j)=sl_e(k);
                            deg_n(j)=deg_e(k);
                            break;
                        end            
                    end
            
                end
                
                %重新初始化概率向量 p，将其所有元素均设置为 1/pn，以准备下一轮的更新。
                p=1/pn*ones(pn,1);
            
                %计算新的步长 sl 和方向 deg 估计，基于加权平均，其中权重由概率向量 p 决定。
                sl=sl_n.'*p;
                deg=deg_n.'*p;

                pos_test = pos(s,:)+[sl*sin(deg),sl*cos(deg)];
                if pos(s, :) == pos_test
                    row = row+2;
                    col = col+2;
                    row_p=6*sigma2/row;%分别计算行数和列数的乘积
                    col_p=6*sigma1/col;
                    pn=row*col;        %计算网格中的总点数
                else
                    break
                end
        end
        %计算新的线性速度估计方差 v_l 和方向速度估计方差 v_d。
        v_l=((sl_n-sl).^2).'*p;
        v_d=((deg_n-deg).^2).'*p;
        pos(s+1,:)=pos(s,:)+[sl*sin(deg),sl*cos(deg)];%在当前位置pos(s,:)上基于新的步长和方向估计，以获得下一个位置           
    end
end