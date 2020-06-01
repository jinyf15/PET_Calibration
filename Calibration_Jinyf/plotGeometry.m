function stop = plotGeometry(x,optimValues,state,fullX,fitIndex,syspara,varargin)
global videoname
global resultfolder
stop = false;
iSave = 1;
fullX(fitIndex==1) = x;
x = fullX;
det_c = zeros(3,16);
PET_c = zeros(3,4);
eulerAng_PET = zeros(3,4);
eulerAng_det = zeros(3,16);
Mtr_PET = zeros(3,3,4);
Mtr_det = zeros(3,3,16);
cla;
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
%cputime % end time
  
%'done2'
    for na=0:syspara.NANG-1
        source_c = PointRotate(source_c0, intersection, RotAxisVector,na*syspara.DANG);
        sp=zeros(syspara.NSX*syspara.NSY*syspara.NSZ,3);
        sp(:,1)=syspara.SDX*reshape(repmat(repmat(-(syspara.NSX-1)/2:(syspara.NSX-1)/2,1,syspara.NSY),1,syspara.NSZ),syspara.NSX*syspara.NSY*syspara.NSZ,1);
        sp(:,2)=syspara.SDY*reshape(repmat(repmat(-(syspara.NSY-1)/2:(syspara.NSY-1)/2,1,syspara.NSX),syspara.NSZ,1),syspara.NSX*syspara.NSY*syspara.NSZ,1);
        sp(:,3)=syspara.SDZ*reshape(repmat(repmat(-(syspara.NSZ-1)/2:(syspara.NSZ-1)/2,syspara.NSX,1),syspara.NSY,1,1),syspara.NSX*syspara.NSY*syspara.NSZ,1);
        sp(sp(:,1).^2+sp(:,2).^2+sp(:,3).^2>syspara.ptRad^2,:)=[];
        sp=sp+source_c';


        % start to calculate probability
        msyspara.mdata=syspara.mdata{na+1};

        
        %txt=['angle: '  num2str(na)];
        txt = ['Iteration',num2str(optimValues.iteration),':fvalue=',num2str(optimValues.fval)];
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
        hold off
        if iSave
                drawnow
                
                frame = getframe(gcf);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                if optimValues.iteration == 0
                    videoname = [resultfolder,'results',num2str(syspara.iter),'iter',datestr(now,30),'.gif'];
                    imwrite(imind,cm,videoname,'gif', 'Loopcount',inf,'DelayTime', 0.1);
                else
                    imwrite(imind,cm,videoname,'gif','WriteMode','append');
                end
        end
        %view(15, 30);
        %view([90 0]);
        %cputime
        %'finished'
    end
end






