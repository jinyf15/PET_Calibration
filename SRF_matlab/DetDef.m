filename = '../data/DetDef';
load('../data/result13_MT/results8iter20200604T073610.mat');
fid = fopen(filename,'wb');
fclose(fid);
fid = fopen(filename,'ab');
syspara.NANG=1;
syspara.DM_DET(1:4)=20;syspara.DN_DET(1:4)=64;
syspara.DDZ_DET(1:4)=0.5;syspara.DDX_DET(1:4)=2*20.9./syspara.DN_DET;syspara.DDY_DET(1:4)=2*20.9./syspara.DN_DET;
syspara.ptRad=1;
syspara.DANG=2*pi/syspara.NANG;
syspara.NDET_DET=syspara.DN_DET.*syspara.DN_DET;
syspara.SDX=0.05;syspara.SDY=0.05;syspara.SDZ=0.05;
syspara.NSX=40; syspara.NSY=40; syspara.NSZ=40;NIMG=syspara.NSX*syspara.NSY*syspara.NSZ;
syspara.mu_det(1:4)=0.058; syspara.fn=6;
[syspara.dxx0, syspara.dyy0] = initDetPixelPosi_v2(syspara.DN_DET,syspara.DDX_DET,syspara.DDY_DET);

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

source_c0=[x(121);x(122);x(123)];
intersection = [x(124:125),0];
theta = x(126);
phi = x(127);
RotAxisVector = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
dp = zeros(3,max(syspara.NDET_DET.*syspara.DM_DET),4);

TotalP = syspara.NDET_DET .* syspara.DM_DET;
SideP = syspara.DN_DET .*syspara.DM_DET;
% PET pixels location (positions are determined by adding shift of the PET center and the rotation of the PET plane)
% corner position
dpc = [-10.45,-10.45,-10.45,-10.45,10.45,10.45,10.45,10.45;
       -10.45,10.45,-10.45,10.45,10.45,-10.45,-10.45,10.45;
       0,0,10,10,0,0,10,10];
tmp = dpc;
dpc = repmat(dpc,1,1,16);
NormalVector = zeros(3,6);
NormalVector(1,1:2) = 1;
NormalVector(2,3:4) = 1;
NormalVector(3,5:6) = 1;
NormalVector = repmat(NormalVector,1,1,16);
pp = zeros(3,max(SideP)*2+max(syspara.DN_DET)^2/4,16);
for i = 1:4
    pp(1,1:SideP(i)/2,4*(i-1)+1:4*(i-1)+4) = -10.45*ones(1,SideP(i)/2,4);
    pp(2,1:SideP(i)/2,4*(i-1)+1:4*(i-1)+4) = repmat(syspara.DDY_DET(i)*(-(syspara.DN_DET(i)/2-1)/2:(syspara.DN_DET(i)/2-1)/2),1,syspara.DM_DET(i),4);
    pp(3,1:SideP(i)/2,4*(i-1)+1:4*(i-1)+4) = reshape(repmat(syspara.DDZ_DET(i)*(0.5:syspara.DM_DET(i)-0.5),syspara.DN_DET(i)/2,1,4),1,syspara.DN_DET(i)/2*syspara.DM_DET(i),4);
    pp(1,SideP(i)/2+1:SideP(i),4*(i-1)+1:4*(i-1)+4) = 10.45*ones(1,SideP(i)/2,4);
    pp(2,SideP(i)/2+1:SideP(i),4*(i-1)+1:4*(i-1)+4) = repmat(syspara.DDY_DET(i)*(-(syspara.DN_DET(i)/2-1)/2:(syspara.DN_DET(i)/2-1)/2),1,syspara.DM_DET(i),4);
    pp(3,SideP(i)/2+1:SideP(i),4*(i-1)+1:4*(i-1)+4) = reshape(repmat(syspara.DDZ_DET(i)*(0.5:syspara.DM_DET(i)-0.5),syspara.DN_DET(i)/2,1,4),1,syspara.DN_DET(i)/2*syspara.DM_DET(i),4);
    pp(1,SideP(i)/2*2+1:SideP(i)/2*3,4*(i-1)+1:4*(i-1)+4) = repmat(syspara.DDX_DET(i)*(-(syspara.DN_DET(i)/2-1)/2:(syspara.DN_DET(i)/2-1)/2),1,syspara.DM_DET(i),4);
    pp(2,SideP(i)/2*2+1:SideP(i)/2*3,4*(i-1)+1:4*(i-1)+4) = -10.45*ones(1,SideP(i)/2,4);
    pp(3,SideP(i)/2*2+1:SideP(i)/2*3,4*(i-1)+1:4*(i-1)+4) = reshape(repmat(syspara.DDZ_DET(i)*(0.5:syspara.DM_DET(i)-0.5),syspara.DN_DET(i)/2,1,4),1,syspara.DN_DET(i)/2*syspara.DM_DET(i),4);
    pp(1,SideP(i)/2*3+1:SideP(i)*2,4*(i-1)+1:4*(i-1)+4) = repmat(syspara.DDX_DET(i)*(-(syspara.DN_DET(i)/2-1)/2:(syspara.DN_DET(i)/2-1)/2),1,syspara.DM_DET(i),4);
    pp(2,SideP(i)/2*3+1:SideP(i)*2,4*(i-1)+1:4*(i-1)+4) = 10.45*ones(1,SideP(i)/2,4);
    pp(3,SideP(i)/2*3+1:SideP(i)*2,4*(i-1)+1:4*(i-1)+4) = reshape(repmat(syspara.DDZ_DET(i)*(0.5:syspara.DM_DET(i)-0.5),syspara.DN_DET(i)/2,1,4),1,syspara.DN_DET(i)/2*syspara.DM_DET(i),4);
    pp(1,SideP(i)*2+1:SideP(i)*2+syspara.NDET_DET(i)/4,4*(i-1)+1:4*(i-1)+4) = repmat(syspara.DDX_DET(i)*(-(syspara.DN_DET(i)/2-1)/2:(syspara.DN_DET(i)/2-1)/2),1,syspara.DN_DET(i)/2,4);
    pp(2,SideP(i)*2+1:SideP(i)*2+syspara.NDET_DET(i)/4,4*(i-1)+1:4*(i-1)+4) = reshape(repmat(syspara.DDY_DET(i)*(-(syspara.DN_DET(i)/2-1)/2:(syspara.DN_DET(i)/2-1)/2),syspara.DN_DET(i)/2,1,4),1,syspara.NDET_DET(i)/4,4);
    pp(3,SideP(i)*2+1:SideP(i)*2+syspara.NDET_DET(i)/4,4*(i-1)+1:4*(i-1)+4) = zeros(1,syspara.NDET_DET(i)/4,4);
end

for i = 1:4   
    for j = 1:4
        dp(1,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i) = repmat(syspara.dxx0(:,i),syspara.DM_DET(i),1)';
        dp(2,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i) = repmat(syspara.dyy0(:,i),syspara.DM_DET(i),1)';
        dp(3,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i) = reshape(repmat(syspara.DDZ_DET(i)*(0.5:syspara.DM_DET(i)-0.5),syspara.NDET_DET(i)/4,1),1,TotalP(i)/4);
        dp(:,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i) = Mtr_det(:,:,j+(i-1)*4)*dp(:,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i)+det_c(:,j+(i-1)*4);
        dpc(:,:,j+(i-1)*4) = Mtr_PET(:,:,i)*(Mtr_det(:,:,j+(i-1)*4)*dpc(:,:,j+(i-1)*4)+det_c(:,j+(i-1)*4))+PET_c(:,i);
        pp(:,:,j+(i-1)*4) = Mtr_PET(:,:,i)*(Mtr_det(:,:,j+(i-1)*4)*pp(:,:,j+(i-1)*4)+det_c(:,j+(i-1)*4))+PET_c(:,i);
        Mtr_det(:,:,j+(i-1)*4) = Mtr_PET(:,:,i) * Mtr_det(:,:,j+(i-1)*4);
        NormalVector(:,:,j+(i-1)*4) = Mtr_det(:,:,j+(i-1)*4) * NormalVector(:,:,j+(i-1)*4);
    end
    dp(:,1:TotalP(i),i) = Mtr_PET(:,:,i)*dp(:,1:TotalP(i),i)+PET_c(:,i);
end
fwrite(fid,Mtr_det,'double');
fwrite(fid,NormalVector,'double');
for i = 1:16
        fwrite(fid,dpc(:,tmp(1,:)<0,i),'double');
        fwrite(fid,dpc(:,tmp(1,:)>0,i),'double');
        fwrite(fid,dpc(:,tmp(2,:)<0,i),'double');
        fwrite(fid,dpc(:,tmp(2,:)>0,i),'double');
        fwrite(fid,dpc(:,tmp(3,:)==0,i),'double');
        fwrite(fid,dpc(:,tmp(3,:)==10,i),'double');
end
fwrite(fid,pp,'double');
fwrite(fid,dp,'double');
fclose(fid);        
        

        
