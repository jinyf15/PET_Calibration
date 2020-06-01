clc
close all

syspara.NANG=1;
syspara.DM_DET(1:4)=20;syspara.DN_DET(1:4)=64;
syspara.DDZ_DET(1:4)=0.5;syspara.DDX_DET(1:4)=2*20.9./syspara.DN_DET;syspara.DDY_DET(1:4)=2*20.9./syspara.DN_DET;
syspara.ptRad=1;
syspara.DANG=2*pi/syspara.NANG;
syspara.NDET_DET=syspara.DN_DET.*syspara.DN_DET;
syspara.SDX=0.05;syspara.SDY=0.05;syspara.SDZ=0.05;
syspara.NSX=40; syspara.NSY=40; syspara.NSZ=40;NIMG=syspara.NSX*syspara.NSY*syspara.NSZ;
syspara.mu_det(1:4)=0.058; syspara.fn=6;
%convert_to_D1('..\data\','CombinedAllEvents_100ns.txt','D1.bin',syspara);
%convert_to_D2('..\data\','D1.bin','D2.bin');
% read experiment data
subset = 8;
for na=0:syspara.NANG-1
    foldername = '..\data\';
    filename = 'D2test3.bin';  
    total_data{na+1}=readdata(foldername,filename,syspara);
    %syspara.mdata{na+1}=[0,0,1,1024*20+32;0,0,1024*20,1024*20*3+1]';
end

%system parameters
[fitX,fitIndex,fullX] = initX_v2;

[syspara.dxx0, syspara.dyy0] = initDetPixelPosi_v2(syspara.DN_DET,syspara.DDX_DET,syspara.DDY_DET);

%%%%% dispaly the initial guess
naa=0; 
syspara.na=0;
syspara.iter=0;
R = [8,5,4,3,2,1,0.75,0.5];
%R = [1,1,1,1,1,0.75,0.5,0.25];
% display the result
% load('100nsfitpara.mat');
% forward_proj_PET_multiang_display_new(fitX, syspara)
% fitX = x;
%load('../data/results8iter20200429T225029');
it = 1;
        syspara.ptRad=R(it);
        syspara.SDX=R(it)/20;syspara.SDY=R(it)/20;syspara.SDZ=R(it)/20;
        syspara.NSX=40; syspara.NSY=40; syspara.NSZ=40;
        x = fullX;
       forward_proj_PET_multiang_new_v3(x, syspara, fitIndex ,fullX,total_data);
  

%%% initial fitting using fmincon (end)

%parpool close




