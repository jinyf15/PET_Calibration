clc
close all
clear all
iSave = 0;
%% unit of length: mm
syspara.NANG=1;
syspara.DM_DET(1:4)=20;syspara.DN_DET(1:4)=64;
syspara.DDZ_DET(1:4)=0.5;syspara.DDX_DET(1:4)=2*20.9./syspara.DN_DET;syspara.DDY_DET(1:4)=2*20.9./syspara.DN_DET;
syspara.ptRad=1;
syspara.DANG=2*pi/syspara.NANG;
syspara.NDET_DET=syspara.DN_DET.*syspara.DN_DET;
syspara.SDX=0.05;syspara.SDY=0.05;syspara.SDZ=0.05;
syspara.NSX=40; syspara.NSY=40; syspara.NSZ=40;NIMG=syspara.NSX*syspara.NSY*syspara.NSZ;
syspara.mu_det(1:4)=0.058; syspara.fn=6;
%convert_to_D1('..\data\','CombinedAllEvents_100ns.txt','D1.bin',syspara);
%convert_to_D2('..\data\','D1.bin','D2.bin');
% read experiment data
subset = 8;
for na=0:syspara.NANG-1
    foldername = '../data/';
    filename = 'D2.bin';  
    total_data{na+1}=readdata(foldername,filename,syspara);
    %syspara.mdata{na+1}=[0,0,1,1024*20+32;0,0,1024*20,1024*20*3+1]';
end
[syspara.dxx0, syspara.dyy0] = initDetPixelPosi_v2(syspara.DN_DET,syspara.DDX_DET,syspara.DDY_DET);
%system parameters
[fitX,fitIndex,fullX] = initX_v2;
x = fullX;
det_c = zeros(3,16);
PET_c = zeros(3,4);
eulerAng_PET = zeros(3,4);
eulerAng_det = zeros(3,16);
Mtr_PET = zeros(3,3,4);
Mtr_det = zeros(3,3,16);

for i = 0:3
    PET_c(:,i+1) = x(1+i*30:3+i*30);
    eulerAng_PET(:,i+1) = x(4+i*30:6+i*30);
    det_c(:,i*4+1:i*4+4) = reshape(x(7+i*30:18+i*30),3,4);
    eulerAng_det(:,i*4+1:i*4+4) = reshape(x(19+i*30:30+i*30),3,4);
	Mtr_PET(:,:,i+1) = euler2TranMatrix(eulerAng_PET(1,i+1),eulerAng_PET(2,i+1),eulerAng_PET(3,i+1));
    for j = 1:4
        Mtr_det(:,:,j+i*4) = euler2TranMatrix(eulerAng_det(1,i*4+j),eulerAng_det(2,i*4+j),eulerAng_det(3,i*4+j));
    end
end

s =syspara.DDX_DET.*syspara.DDY_DET;
dp = zeros(3,max(syspara.NDET_DET.*syspara.DM_DET),4);
dpc= zeros(3,max(syspara.NDET_DET)*max(syspara.DM_DET)*4,3);
TotalP = syspara.NDET_DET .* syspara.DM_DET;
% PET pixels location (positions are determined by adding shift of the PET center and the rotation of the PET plane)

for i = 1:4   
    for j = 1:4
        dp(1,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i) = repmat(syspara.dxx0(:,i),syspara.DM_DET(i),1)';
        dp(2,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i) = repmat(syspara.dyy0(:,i),syspara.DM_DET(i),1)';
        dp(3,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i) = reshape(repmat(syspara.DDZ_DET(i)*(0.5:syspara.DM_DET(i)-0.5),syspara.NDET_DET(i)/4,1),1,TotalP(i)/4);
        dp(:,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i) = Mtr_det(:,:,j+(i-1)*4)*dp(:,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i)+det_c(:,j+(i-1)*4);     
    end
    dp(:,1:TotalP(i),i) = Mtr_PET(:,:,i)*dp(:,1:TotalP(i),i)+PET_c(:,i);
end
na = 0;
data=total_data{na+1};

SDX = 100;SDY = 100;SDZ = 100;
Backprojection = zeros(SDX,SDY,SDZ);
l = size(data,2);
for i = 1:l
    ind = find(data(:,i)~=0);
    dn1=data(ind(1),i);
    dn2=data(ind(2),i);
    x1 = dp(:,dn1,ind(1));
    x2 = dp(:,dn2,ind(2));
    if ind(1)==1 && ind(2)==2
        scan = 'Y';
    else 
        scan = 'X';
    end
    for j = 1:100
        index = (j-50.5)/2;
        [x,y,z] = CalIntersect(x1,x2,scan,index);
        Backprojection(x,y,z) = Backprojection(x,y,z)+1;
    end
end
location=find(Backprojection==max(Backprojection(:)));
x = (mod(location,100)-50.5)/2
y = (floor(mod(location,10000)/100)-50.5)/2
z = (floor(location/10000)-50.5)/2