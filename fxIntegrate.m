function result=fxIntegrate(little_limit,superior_limit,data)
%little_limit：化为标准正态分布后的下限
%superior_limit：化为标准正态分布后的上限
%data:数据
%result：积分结果
ll=abs(round(little_limit*1000))+1;
sl=abs(round(superior_limit*1000))+1;
result=data(sl)-data(ll);
