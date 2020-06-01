clc
clear all
close all
iSave = 1;
global resultfolder
resultfolder = '../data/result10/';
if ~exist(resultfolder,'dir')
    mkdir(resultfolder);
end
%% unit of length: mm
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
convert_to_D2('../data/','D1.bin','D2_large.bin');
% read experiment data
subset = 4;
for na=0:syspara.NANG-1
    foldername = '..\data\';
    filename = 'D2.bin';  
    total_data{na+1}=readdata(foldername,filename,syspara);
    %syspara.mdata{na+1}=[0,0,1,1024*20+32;0,0,1024*20,1024*20*3+1]';
end

%system parameters
[fitX,fitIndex,fullX] = initX_v2;
[A,b,Aeq,beq,lb,ub] = constrains_v2(fullX,fitIndex);
nonlcon= @(x)ellipseparabola(x, fullX,fitIndex);

[syspara.dxx0, syspara.dyy0] = initDetPixelPosi_v2(syspara.DN_DET,syspara.DDX_DET,syspara.DDY_DET);

%%%%% dispaly the initial guess
naa=0; 
syspara.na=0;
syspara.iter=0;
%R = [6,5,4,3,2,1,0.75,0.5];
R = [4,3,2,1];
%R = [1,1,1,1,1,0.75,0.5,0.25];
% display the result
% load('100nsfitpara.mat');
% forward_proj_PET_multiang_display_new(fitX, syspara)
% fitX = x;
forward_proj_PET_multiang_display_new_v2(fullX, fitIndex,syspara,total_data,0);

%%% initial fitting using fmincon (start)

%options=optimset('LargeScale','on','TolFun',1,'TolX',1e-5,'MaxIter',50,'MaxFunEvals',length(fitX)*2000,'Display','iter','PlotFcns',{@optimplotx,@optimplotfval});
%options = optimset('Display','iter','MaxIter',10, 'DiffMinChange', 1e-8, 'TolFun',1e-8, 'TolCon', 1e-8, 'Algorithm', 'sqp', 'Useparallel', 'always');%, 'Outputfcn', @outfun);
start_t = datestr(now);
for it=1:subset
        if it > 0
            syspara.mdata{:} = total_data{:}(:,it:4:end);
        else
            syspara.mdata{:} = total_data{:}(:,it:8:end);
        end
        %syspara.mdata = total_data;
        
        syspara.na=naa;
        syspara.iter=it;
        syspara.ptRad=R(it);
        syspara.SDX=R(it)/20;syspara.SDY=R(it)/20;syspara.SDZ=R(it)/20;
        syspara.NSX=40; syspara.NSY=40; syspara.NSZ=40;
        myplot = @(x,optimValues,state)plotGeometry(x,optimValues,state,fullX,fitIndex,syspara);
        
        %% using standard fmincon
        if it>0
            options=optimset('LargeScale','on','TolFun',1e-10,'TolX',1e-5,'MaxIter',50,'MaxFunEvals',length(fitX)*2000,'Display','iter','Algorithm','sqp', 'Useparallel', 'always','PlotFcns',{@optimplotx,...
   @optimplotfval,myplot});
            [x,fval,exitflag,output]=fmincon(@(x)forward_proj_PET_multiang_new_v2(x, syspara, fitIndex ,fullX), fitX, [], [], [], [],lb, ub, nonlcon,options);
        %fitX=x;
        else
            options=optimset('LargeScale','on','TolFun',1,'TolX',1e-5,'MaxIter',20,'MaxFunEvals',length(fitX)*2000,'Display','iter','Algorithm','sqp', 'Useparallel', 'always','PlotFcns',{@optimplotx,...
   @optimplotfval,myplot});
            fwdproj_multiang=@(x)forward_proj_PET_multiang_new_v2(x, syspara, fitIndex ,fullX);
            %options = optimset('Display','iter', 'Algorithm', 'sqp', 'Useparallel', 'always');
            problem=createOptimProblem('fmincon', 'objective', fwdproj_multiang, 'x0', fitX, 'lb', lb, 'ub', ub,'nonlcon',nonlcon,'options', options); 
            %problem=createOptimProblem('fmincon', 'objective', fwdproj_multiang, 'x0', fitX, 'lb', lb, 'ub', ub,'options', options); 
            gs = GlobalSearch;
            [x,f] = run(gs,problem);
%             pts = repmat(fitX,8,1);
%             ms = MultiStart;
%             for i = 1:8
%                 pts(i,121:123)=[-10*rand 10*rand -10*rand];
%             end
%             tpoints = CustomStartPointSet(pts);
%             [x,fval] = run(ms,problem,tpoints);
        end
        fitX=x;
        %% display the result
        fullX(fitIndex==1)=x;
        forward_proj_PET_multiang_display_new_v2(fullX,fitIndex, syspara,total_data,iSave);
        if iSave
            filename=[resultfolder,'/results',num2str(it),'iter',datestr(now,30)];
            save(filename,'fitIndex','fullX');
            saveas(gcf,filename,'fig');
        end
        end_t(it,:)=datestr(now)
        
        
        pause(1);
end
[start_t; end_t]
fit_data_multiang_fmincon=x;
%%% initial fitting using fmincon (end)

%parpool close




