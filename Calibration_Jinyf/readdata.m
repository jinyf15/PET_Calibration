function [output] = readdata(foldername,filename,syspara)


fid = fopen([foldername,filename],'r');
nn=1;
B = zeros(4,1);
while ~feof(fid)
    A = fread(fid,12,'double');          
    if isempty(A)
        break;
    end
    B(floor((A(2)-1)/4)+1,nn) = A(6) + mod(A(2)-1,4)*syspara.DM_DET(floor((A(2)-1)/4)+1)*syspara.NDET_DET(floor((A(2)-1)/4)+1)/4;
    B(floor((A(8)-1)/4)+1,nn) = A(12) + mod(A(8)-1,4)*syspara.DM_DET(floor((A(8)-1)/4)+1)*syspara.NDET_DET(floor((A(8)-1)/4)+1)/4;
    nn=nn+1;
end
fclose(fid);
output = B;

