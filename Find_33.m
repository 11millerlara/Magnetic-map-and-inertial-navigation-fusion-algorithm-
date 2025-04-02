function result=Find_33(DataSet,x,y)
result=[-1,-1];
M=size(DataSet,1);
N=size(DataSet,2);
%x=round((x-0.3416)/0.001466)+1;
%y=round((y-0.9)/0.00095)+1;
x=round((x-34.2)/0.1)+1;
y=round((y-90)/0.1)+1;
if(x>M||x<=0||y>N||y<=0)
    return
else
    result=[x,y];
end
