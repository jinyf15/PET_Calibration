function [outputArg] = preprocessing(foldername,filename,syspara)
%preprossing: 1)bin the continuous data to virtual voxels
%             2)pick 2 and energy around 511keV events
centerEnergy = 511;
energyWindow = 5;
range = [centerEnergy-energyWindow/2 centerEnergy+energyWindow/2]; 
fid = fopen([foldername,filename],'r');
formatSpec = ['%d' '%d' '%f' '%f' '%f' '%f\n'];
sizeA = [6 Inf];
A = fscanf(fid, formatSpec, sizeA);
fclose(fid);
A = A';
A = A(A(:,1)==2,:);
datalength = size(A,1)/2;
B = zeros(4,datalength*2);
nn=1;
M = zeros(3,3,4);
M(:,:,1)=euler2TranMatrix(0,-pi/2,-pi/2);
M(:,:,2)=euler2TranMatrix(0,-pi/2,pi/2);
M(:,:,3)=euler2TranMatrix(0,-pi/2,pi);
M(:,:,4)=euler2TranMatrix(0,-pi/2,0);
for i = 1:datalength
    if A(2*i-1,3)>range(1)&&A(2*i-1,3)<range(2)&&A(2*i,3)>range(1)&&A(2*i,3)<range(2)
        x1 = A(2*i-1,4);
        y1 = A(2*i-1,5);
        x2 = A(2*i,4);
        y2 = A(2*i,5);
        num(1) = uint8(absmax(x1,y1)>0)+uint8(abs(x1)>abs(y1))*2+1;
        num(2) = uint8(absmax(x2,y2)>0)+uint8(abs(x2)>abs(y2))*2+1;
        if num(1)==num(2) %2 events are on the same PET which is chance coincidnece
            continue;
        else %bin
            for j = 1:2
                x = A(2*i-2+j,4:6)*M(:,:,num(j));
                x(3)=x(3)+50;
                if abs(x(1))>=21.7 || abs(x(2))>=21.7||abs(x(1))<=0.8||abs(x(2))<=0.8  %exclude events outside detectors                    
                    break;
                end
                x(x(1:2)>0)=x(x(1:2)>0)-1.6;
                x(1:2)=x(1:2)+21.7;               
                B(num(j),nn) = floor(x(1)/syspara.DDX_DET(num(j)))+1+floor(x(2)/syspara.DDY_DET(num(j)))*syspara.DN_DET(num(j))+floor(x(3)/syspara.DDZ_DET(num(j)))*syspara.NDET_DET(num(j));                
            end            
            if length(find(B(:,nn)==0))~=2
                B(:,nn)=0;
                continue;
            end
            nn = nn + 1;
        end
    end
end
%outputArg = B(:,any(B));
outputArg = B(:,(B(1,:)&B(2,:))|(B(3,:)&B(4,:)));

