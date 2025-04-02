function result=Integrate(little_limit,superior_limit,data)
%little_limit：化为标准正态分布后的下限
%superior_limit：化为标准正态分布后的上限
%data:数据
%result：积分结果
ll=round(little_limit*1000);
sl=round(superior_limit*1000);
if ll>=0
    result=data(sl+1)-data(ll+1);
elseif sl>=0&&ll<0
    result=data(sl+1)+data(abs(ll)+1)-1;
elseif sl<0
    result=data(abs(ll)+1)-data(abs(sl)+1);
end
% result=result*100;