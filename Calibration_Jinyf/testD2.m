
foldername = 'D:\UIUC\RA\PET\Calibration\data\';
filename = 'D2test.bin';
% 
 fid2 = fopen([foldername,filename],'w');
 fclose(fid2);
fid2 = fopen([foldername,filename],'a');

D2 = zeros(1,12);
for N = 1:20
    ind1 = 1024+1024*(N-1);
    ind2 = 32+1024*(N-1);
    D2(1:6) = [N 1 0 0 0 ind1];
    D2(7:12) = [N 7 0 0 0 ind2];
    fwrite(fid2,D2,'double');
end
fclose(fid2);
          