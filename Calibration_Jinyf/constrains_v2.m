function [A,b,Aeq,beq,lb,ub]= constrains_v2(x,fitIndex)
A=[];
b=[];
% A=zeros(size(x,2),12);
% b=ones(1,12);
% for i = 1:3
%     for j = 0:2:2
%         A(i+j*30,i+j*3) = 1;
%         A(i+j*30+30,i+j*3) = -1;
%         A(i+j*30,i+j*3+3) = -1;
%         A(i+j*30+30,i+j*3+3) = 1;
%     end
% end
Aeq=[];
beq=[];

ub = zeros(1,size(x,2));
lb = zeros(1,size(x,2));

lb(121:123) = x(121:123) - 10;
ub(121:123) = x(121:123) + 10;

for i = 0:3
    lb(1+i*30:3+i*30) = x(1+i*30:3+i*30) - 2;
    ub(1+i*30:3+i*30) = x(1+i*30:3+i*30) + 2;
    lb(4+i*30:6+i*30) = x(4+i*30:6+i*30) - 1/180*pi;
    ub(4+i*30:6+i*30) = x(4+i*30:6+i*30) + 1/180*pi;
    lb(7+i*30:18+i*30) = x(7+i*30:18+i*30) - 0.5;
    ub(7+i*30:18+i*30) = x(7+i*30:18+i*30) + 0.5;
    lb(19+i*30:30+i*30) = x(19+i*30:30+i*30) - 0.5/180*pi;
    ub(19+i*30:30+i*30) = x(19+i*30:30+i*30) + 0.5/180*pi;
    
end

lb(fitIndex==0) = [];
ub(fitIndex==0) = [];