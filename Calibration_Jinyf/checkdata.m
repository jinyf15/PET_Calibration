syspara.NANG=1;
syspara.DM_DET(1:4)=20;syspara.DN_DET(1:4)=64;
syspara.DDZ_DET(1:4)=0.5;syspara.DDX_DET(1:4)=2*20.9./syspara.DN_DET;syspara.DDY_DET(1:4)=2*20.9./syspara.DN_DET;
syspara.ptRad=1;
syspara.DANG=2*pi/syspara.NANG;
syspara.NDET_DET=syspara.DN_DET.*syspara.DN_DET;
syspara.SDX=0.1;syspara.SDY=0.1;syspara.SDZ=0.1;
syspara.NSX=20; syspara.NSY=20; syspara.NSZ=20;NIMG=syspara.NSX*syspara.NSY*syspara.NSZ;
syspara.mu_det(1:4)=0.58; syspara.fn=6;
folder = '..\data\';
file1 = 'D2.bin';
file2 = 'CombinedAllEvents_100ns.txt';
X1 = readdata(folder,file1,syspara);
X2 = preprocessing(folder,file2,syspara);
for i = 1:size(X1,2)
    ind1 = mod(X1(:,i)-1,1024);
    ind2 = floor((X1(:,i)-1)/(1024*20));
    depth = floor(mod(X1(:,i)-1,1024*20)/1024);
    X3(:,i) = (ind1+1024*ind2) + depth * 1024*4+1;
end
for i = 1:size(X3,2)
    ind1 = mod(X3(:,i)-1,1024);
    depth = floor((X3(:,i)-1)/4096);
    ind2 = floor(mod(X3(:,i)-1,4096)/1024);
    X4(:,i) = ind2*1024*20+depth*1024+ind1+1;
end