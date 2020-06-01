function total_data = ExcludeChanceCoin(syspara, total_data, fullX)

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
%source_c0=[0;R;0];
source_c0=[x(121);x(122);x(123)];
intersection = [x(124:125),0];
theta = x(126);
phi = x(127);
RotAxisVector = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

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
        if i < 4
            dpc(1,1+TotalP(i)*(j-1):4:TotalP(i)*j,i)=dp(1,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i)-syspara.DDX_DET(i)/2;
            dpc(1,2+TotalP(i)*(j-1):4:TotalP(i)*j,i)=dp(1,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i)+syspara.DDX_DET(i)/2;
            dpc(1,3+TotalP(i)*(j-1):4:TotalP(i)*j,i)=dp(1,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i)+syspara.DDX_DET(i)/2;
            dpc(1,4+TotalP(i)*(j-1):4:TotalP(i)*j,i)=dp(1,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i)-syspara.DDX_DET(i)/2;
            dpc(2,1+TotalP(i)*(j-1):4:TotalP(i)*j,i)=dp(2,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i)-syspara.DDY_DET(i)/2;
            dpc(2,2+TotalP(i)*(j-1):4:TotalP(i)*j,i)=dp(2,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i)-syspara.DDY_DET(i)/2;
            dpc(2,3+TotalP(i)*(j-1):4:TotalP(i)*j,i)=dp(2,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i)+syspara.DDY_DET(i)/2;
            dpc(2,4+TotalP(i)*(j-1):4:TotalP(i)*j,i)=dp(2,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i)+syspara.DDY_DET(i)/2;
            dpc(3,1+TotalP(i)*(j-1):TotalP(i)*j,i)=reshape(repmat(dp(3,1+TotalP(i)/4*(j-1):TotalP(i)/4*j,i),4,1),TotalP(i),1);
            % rotate and shift
            dpc(:,1+TotalP(i)*(j-1):TotalP(i)*j,i) = Mtr_PET(:,:,i)*dpc(:,1+TotalP(i)*(j-1):TotalP(i)*j,i)+PET_c(:,i);
        end
    end
    dp(:,1:TotalP(i),i) = Mtr_PET(:,:,i)*dp(:,1:TotalP(i),i)+PET_c(:,i);
end

pt=zeros(2,syspara.fn*syspara.fn,4);
pt(1,:,:) = repmat(repmat(-(syspara.fn-1)/2:(syspara.fn-1)/2,1,syspara.fn)',1,4);
pt(2,:,:) = repmat(reshape(repmat(-(syspara.fn-1)/2:(syspara.fn-1)/2,syspara.fn,1),1,syspara.fn*syspara.fn),4,1)';
pt(1,:,:) = pt(1,:,:).*reshape(repmat(syspara.DDX_DET,syspara.fn*syspara.fn,1),1,syspara.fn*syspara.fn,4)/syspara.fn;
pt(2,:,:) = pt(2,:,:).*reshape(repmat(syspara.DDY_DET,syspara.fn*syspara.fn,1),1,syspara.fn*syspara.fn,4)/syspara.fn;

reqpara = zeros(5,4);
reqpara(1,:) = syspara.DM_DET;
reqpara(2,:) = syspara.DDZ_DET;
reqpara(3,:) = syspara.NDET_DET;
reqpara(4,:) = syspara.mu_det;
reqpara(5,:) = s;

diameter = max(syspara.ptRad/4/2.355,0.25/2.355);
for na=0:syspara.NANG-1
    source_c = PointRotate(source_c0, intersection, RotAxisVector,na*syspara.DANG);
    sp=zeros(syspara.NSX*syspara.NSY*syspara.NSZ,3);
    sp(:,1)=syspara.SDX*reshape(repmat(repmat(-(syspara.NSX-1)/2:(syspara.NSX-1)/2,1,syspara.NSY),1,syspara.NSZ),syspara.NSX*syspara.NSY*syspara.NSZ,1);
    sp(:,2)=syspara.SDY*reshape(repmat(repmat(-(syspara.NSY-1)/2:(syspara.NSY-1)/2,1,syspara.NSX),syspara.NSZ,1),syspara.NSX*syspara.NSY*syspara.NSZ,1);
    sp(:,3)=syspara.SDZ*reshape(repmat(repmat(-(syspara.NSZ-1)/2:(syspara.NSZ-1)/2,syspara.NSX,1),syspara.NSY,1,1),syspara.NSX*syspara.NSY*syspara.NSZ,1);
    sp(sp(:,1).^2+sp(:,2).^2+sp(:,3).^2>syspara.ptRad^2,:)=[];
    spc= 1e5/(diameter)*exp(-(sp(:,1).^2+sp(:,2).^2+sp(:,3).^2)/(2*(diameter)^2));
    sp=sp+source_c';
    
    % start to calculate probability
    data=total_data{na+1};
    l = size(data,2);
    index = 1;
    while index <= l
        data=total_data{na+1}(:,index);
        logpr = forward_proj_multi_angle_prob(data, size(data,2), size(sp,1), dp, dpc, sp(:,1),sp(:,2),sp(:,3),...
            Mtr_PET, pt, reqpara, syspara.fn, spc);
        if logpr == -100 % chance coincidence
            total_data{na+1}(:,index) = [];
            index = index-1;
            l=l-1;
        end
        index = index+1;
    end
        
end
end







