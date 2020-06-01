function logpr=forward_proj_PET_multiang_new(x, syspara,fitIndex,fullX)
fullX(fitIndex==1) = x;
x = fullX;
R = x(1); % rotation radius
det_c = zeros(3,41);
eulerAng = zeros(3,41);
Mtr = zeros(3,3,4);
gap = zeros(4,41);
for i = 0:3
    det_c(:,i+1) = x(2+i*10:4+i*10);
    eulerAng(:,i+1) = x(5+i*10:7+i*10);
	Mtr(:,:,i+1) = euler2TranMatrix(eulerAng(1,i+1),eulerAng(2,i+1),eulerAng(3,i+1));
    gap(:,i+1) = x(8+i*10:11+i*10);
end
%source_c0=[0;R;0];
source_c0=[x(1);x(42);x(43)];
s =syspara.DDX_DET.*syspara.DDY_DET;
dp = zeros(3,max(syspara.NDET_DET.*syspara.DM_DET),4);
dpc= zeros(3,max(syspara.NDET_DET)*max(syspara.DM_DET)*4,3);
% PET pixels location (positions are determined by adding shift of the PET center and the rotation of the PET plane)
for i = 1:4   
    
    dp(1,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i) = repmat(syspara.dxx0(:,i),syspara.DM_DET(i),1)';
    dp(2,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i) = repmat(syspara.dyy0(:,i),syspara.DM_DET(i),1)';
    dp(3,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i) = reshape(repmat(syspara.DDZ_DET(i)*(0.5:syspara.DM_DET(i)-0.5),syspara.NDET_DET(i),1),1,syspara.NDET_DET(i)*syspara.DM_DET(i));
    
    dp(1,dp(1,:,i)<0,i)= dp(1,dp(1,:,i)<0,i)-gap(1,i);
    dp(1,dp(1,:,i)>0,i)= dp(1,dp(1,:,i)>0,i)+gap(2,i);
    dp(2,dp(2,:,i)<0,i)= dp(2,dp(2,:,i)<0,i)-gap(3,i);
    dp(2,dp(2,:,i)>0,i)=  dp(2,dp(2,:,i)>0,i)+gap(4,i);
    if i < 4
        dpc(1,1:4:syspara.NDET_DET(i)*syspara.DM_DET(i)*4,i)=dp(1,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i)-syspara.DDX_DET(i)/2;
        dpc(1,2:4:syspara.NDET_DET(i)*syspara.DM_DET(i)*4,i)=dp(1,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i)+syspara.DDX_DET(i)/2;
        dpc(1,3:4:syspara.NDET_DET(i)*syspara.DM_DET(i)*4,i)=dp(1,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i)+syspara.DDX_DET(i)/2;
        dpc(1,4:4:syspara.NDET_DET(i)*syspara.DM_DET(i)*4,i)=dp(1,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i)-syspara.DDX_DET(i)/2;
        dpc(2,1:4:syspara.NDET_DET(i)*syspara.DM_DET(i)*4,i)=dp(2,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i)-syspara.DDY_DET(i)/2;
        dpc(2,2:4:syspara.NDET_DET(i)*syspara.DM_DET(i)*4,i)=dp(2,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i)-syspara.DDY_DET(i)/2;
        dpc(2,3:4:syspara.NDET_DET(i)*syspara.DM_DET(i)*4,i)=dp(2,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i)+syspara.DDY_DET(i)/2;
        dpc(2,4:4:syspara.NDET_DET(i)*syspara.DM_DET(i)*4,i)=dp(2,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i)+syspara.DDY_DET(i)/2;
        dpc(3,1:syspara.NDET_DET(i)*syspara.DM_DET(i)*4,i)=reshape(repmat(dp(3,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i),4,1),syspara.NDET_DET(i)*syspara.DM_DET(i)*4,1);
        % rotate and shift
        dpc(:,1:syspara.NDET_DET(i)*syspara.DM_DET(i)*4,i) = Mtr(:,:,i)*dpc(:,1:syspara.NDET_DET(i)*syspara.DM_DET(i)*4,i)+det_c(:,i);
    end
    dp(:,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i) = Mtr(:,:,i)*dp(:,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i)+det_c(:,i);
end
%cputime % end time
  
%'done2'


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



for na=0:syspara.NANG-1
    source_c = euler2TranMatrix(0,0,na*syspara.DANG)* source_c0;
    sp=zeros(syspara.NSX*syspara.NSY*syspara.NSZ,3);
    sp(:,1)=syspara.SDX*reshape(repmat(repmat(-(syspara.NSX-1)/2:(syspara.NSX-1)/2,1,syspara.NSY),1,syspara.NSZ),syspara.NSX*syspara.NSY*syspara.NSZ,1);
    sp(:,2)=syspara.SDY*reshape(repmat(repmat(-(syspara.NSY-1)/2:(syspara.NSY-1)/2,1,syspara.NSX),syspara.NSZ,1),syspara.NSX*syspara.NSY*syspara.NSZ,1);
    sp(:,3)=syspara.SDZ*reshape(repmat(repmat(-(syspara.NSZ-1)/2:(syspara.NSZ-1)/2,syspara.NSX,1),syspara.NSY,1,1),syspara.NSX*syspara.NSY*syspara.NSZ,1);
    sp(sp(:,1).^2+sp(:,2).^2+sp(:,3).^2>syspara.ptRad,:)=[];
    spc= exp(-2.5*(sp(:,1).^2+sp(:,2).^2+sp(:,3).^2));
    sp=sp+source_c';
    
    % start to calculate probability
    data=syspara.mdata{na+1};
%% 

    logpr = forward_proj_multi_angle_prob(data, size(data,2), size(sp,1), dp, dpc, sp(:,1),sp(:,2),sp(:,3),...
        Mtr, pt, reqpara, syspara.fn, spc);
    % toc;
    % [na logpr]
end
logpr=-logpr;
%'done4'

% (MLJ, 06252012) display the current geometry

end




