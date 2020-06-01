function convert_to_D1(foldername,filename1,filename2,syspara)

formatSpec = ['%d\t%d\t%f\t%f\t%f\t%f\n'];
sizeA = [6 1e6];

M = zeros(3,3,4);
M(:,:,1)=euler2TranMatrix(0,-pi/2,-pi/2);
M(:,:,2)=euler2TranMatrix(0,-pi/2,pi/2);
M(:,:,3)=euler2TranMatrix(0,-pi/2,pi);
M(:,:,4)=euler2TranMatrix(0,-pi/2,0);
tic;
fileID = fopen([foldername,filename1],'r');
fid2 = fopen([foldername,filename2],'w');
fclose(fid2);
fid2 = fopen([foldername,filename2],'a');
B = zeros(1,1e6*9);
N = 0;
nnn=0;
while ~feof(fileID)
    A = fscanf(fileID, formatSpec, sizeA);           
    if isempty(A)
        break;
    end
    sizeA2 = size(A,2);
    nn=1;
    mm=1;
    while nn<=sizeA2
        if nn+A(1,nn)-1 > sizeA2 
            A(:,end+1:nn+A(1,nn)-1) = fscanf(fileID, formatSpec, [6 -sizeA2+nn+A(1,nn)-1]);                                                    
        end
        iwrite = 1;
        temp = zeros(A(1,nn),8);
        for i = 0:A(1,nn)-1           
            x = A(4,nn+i);
            y = A(5,nn+i);
            ii = double(absmax(x,y)>0)+double(abs(x)>abs(y))*2+1;
            xx = A(4:6,nn+i)'*M(:,:,ii);
            xx(3)=xx(3)+50;
            if abs(xx(1))>=21.7 || abs(xx(2))>=21.7||abs(xx(1))<=0.8||abs(xx(2))<=0.8  %exclude events outside detectors                    
                iwrite = 0;
                break;
            end
            jj = double(xx(1)>0) + double(xx(2)>0)*2 + 1;
            temp(i+1,1) = N+1;
            temp(i+1,2) = i+1;
            temp(i+1,3) = (ii-1)*4 + jj;
            temp(i+1,4) = xx(1) - 11.25 * (double(xx(1)>0)*2-1);
            temp(i+1,5) = xx(2) - 11.25 * (double(xx(2)>0)*2-1);
            temp(i+1,6) = xx(3);
            temp(i+1,7) = floor((temp(i+1,4)+10.45)/syspara.DDX_DET(ii))+1+floor((temp(i+1,5)+10.45)/syspara.DDY_DET(ii))*syspara.DN_DET(ii)/2+floor(temp(i+1,6)/syspara.DDZ_DET(ii))*syspara.NDET_DET(ii)/4;                
            temp(i+1,8) = A(3,nn+i);
            temp(1:i+1,:) = check_NG(temp(1:i+1,:));
        end
        if iwrite
            N = N+1;
            B(mm) = A(1,nn);
            mm = mm+1;
            B(1,mm:mm+8*A(1,nn)-1)=reshape(temp',1,8*A(1,nn));
            mm=mm+A(1,nn)*8;
        end
        nn=nn+A(1,nn);
    end
    fwrite(fid2,B(1:mm-1),'double');
    nnn=nnn+1;
    disp(['complete processing ',num2str(nnn*1e6),' data ...']);
    time = toc

end
    fclose(fileID);
    fclose(fid2);
end



