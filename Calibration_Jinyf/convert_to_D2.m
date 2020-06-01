function convert_to_D2(foldername,filename1,filename2)
Energy = 511;
EnergyWindow = 5;
%range = [Energy-EnergyWindow/2 Energy+EnergyWindow/2];
%range = [461 515];
range = [300 515];
tic;
fileID = fopen([foldername,filename1],'r');
fid2 = fopen([foldername,filename2],'w');
fclose(fid2);
fid2 = fopen([foldername,filename2],'a');
N = 0;
nnn=0;
D2 = zeros(1,12);
while ~feof(fileID)
    A = fread(fileID,1,'double');          
    if isempty(A)
        break;
    end
    B = fread(fileID,A*8,'double');
    B=reshape(B(1:A*8),8,A);
    if (A==2) && (B(2,1)~=B(2,2)) 
        if (B(8,1)<=range(2)&&B(8,1)>=range(1)) && (B(8,2)>=range(1)&&B(8,2)<=range(2))
            if floor((B(3,1)-1)/4)+floor((B(3,2)-1)/4)==1 || floor((B(3,1)-1)/4)+floor((B(3,2)-1)/4)==5
                N=N+1;
                D2(1) = N;D2(2:6) = B(3:7,1)';
                D2(7) = N;D2(8:12) = B(3:7,2)';
                fwrite(fid2,D2,'double');
                if mod(N,10) ==0
                    disp(['complete converting ',num2str(N),' coincidence pairs ...']);
                    time = toc
                end
            end
        end
    end
    

end
    fclose(fileID);
    fclose(fid2);
end







