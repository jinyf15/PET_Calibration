function forward_proj_PET_multiang_display_new(x, syspara)

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
source_c0=[x(1);x(42);x(43)];

dp = zeros(3,max(syspara.NDET_DET.*syspara.DM_DET),4);

% PET pixels location (positions are determined by adding shift of the PET center and the rotation of the PET plane)
for i = 1:4   
    
    dp(1,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i) = repmat(syspara.dxx0(:,i),syspara.DM_DET(i),1)';
    dp(2,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i) = repmat(syspara.dyy0(:,i),syspara.DM_DET(i),1)';
    dp(3,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i) = reshape(repmat(syspara.DDZ_DET(i)*(0.5:syspara.DM_DET(i)-0.5),syspara.NDET_DET(i),1),1,syspara.NDET_DET(i)*syspara.DM_DET(i));
    
    dp(1,dp(1,:,i)<0,i)= dp(1,dp(1,:,i)<0,i)-gap(1,i);
    dp(1,dp(1,:,i)>0,i)= dp(1,dp(1,:,i)>0,i)+gap(2,i);
    dp(2,dp(2,:,i)<0,i)= dp(2,dp(2,:,i)<0,i)-gap(3,i);
    dp(2,dp(2,:,i)>0,i)=  dp(2,dp(2,:,i)>0,i)+gap(4,i);
    dp(:,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i) = Mtr(:,:,i)*dp(:,1:syspara.NDET_DET(i)*syspara.DM_DET(i),i)+det_c(:,i);
end
%cputime % end time
  
%'done2'
    for na=0:syspara.NANG-1
        source_c = euler2TranMatrix(0,0,na*syspara.DANG)* source_c0;
        sp=zeros(syspara.NSX*syspara.NSY*syspara.NSZ,3);
        sp(:,1)=syspara.SDX*reshape(repmat(repmat(-(syspara.NSX-1)/2:(syspara.NSX-1)/2,1,syspara.NSY),1,syspara.NSZ),syspara.NSX*syspara.NSY*syspara.NSZ,1);
        sp(:,2)=syspara.SDY*reshape(repmat(repmat(-(syspara.NSY-1)/2:(syspara.NSY-1)/2,1,syspara.NSX),syspara.NSZ,1),syspara.NSX*syspara.NSY*syspara.NSZ,1);
        sp(:,3)=syspara.SDZ*reshape(repmat(repmat(-(syspara.NSZ-1)/2:(syspara.NSZ-1)/2,syspara.NSX,1),syspara.NSY,1,1),syspara.NSX*syspara.NSY*syspara.NSZ,1);
        sp(sp(:,1).^2+sp(:,2).^2+sp(:,3).^2>syspara.ptRad^2,:)=[];
        sp=sp+source_c';


        % start to calculate probability
        msyspara.mdata=syspara.mdata{na+1};

        figure,clf;
        txt=['angle: '  num2str(na)];
        title(txt); hold on;

        for ii=1:size(msyspara.mdata,2)
            ind = find(msyspara.mdata(:,ii)~=0);
            dn1=msyspara.mdata(ind(1),ii);
            dn2=msyspara.mdata(ind(2),ii);
            plot3([dp(1,dn1,ind(1)), dp(1,dn2,ind(2))], [dp(2,dn1,ind(1)), dp(2,dn2,ind(2))], [dp(3,dn1,ind(1)), dp(3,dn2,ind(2))],'k');
        end

        hold on; % plot the locations of each pixels
        for kk=1:4
            plot3(dp(1,:,kk), dp(2,:,kk), dp(3,:,kk), '.');hold on;
        end
        hold on; % current source position
        plot3(sp(:,1), sp(:,2), sp(:,3), 'r.');
        axis equal
        %view(15, 30);
        %view([90 0]);
        %cputime
        %'finished'
        pause(1);
    end
end




