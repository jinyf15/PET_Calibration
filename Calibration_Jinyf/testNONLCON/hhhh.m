clc;
clear all
tmax=440;
cmin=160;
A=[2 0 0 1;0 -2 -2 0;-1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 -1];
b=[tmax;-cmin;-10;-10;-10;-50];
x0=[20;50;50;400];
myfun=@(x)(3000*x(1)+2*x(4)*x(2)+2*x(4)*x(3));
nonlcon = @mycon;
[x,fval,exitflag,output,lambda,grad,hession]=fmincon(myfun,x0,A,b,[],[],[],[],nonlcon);