clc
clear all
close all
%% settings
iSave = 1;
iRead = 0;
subset = 8;
R0 = 4;
R = [3.5,3,2.5,2,1.5,1,0.75,0.5];
global resultfolder
resultfolder = '../data/result11_MT/';
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
%convert_to_D2('..\data\','D1.bin','D2.bin');
% read experiment data
for na=0:syspara.NANG-1
    foldername = '../data/';
    filename = 'D2_large.bin';  
    total_data{na+1}=readdata(foldername,filename,syspara);
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
forward_proj_PET_multiang_display_new_v2(fullX, fitIndex,syspara,total_data,0); % display

%% initial fitting using fmincon (start)
start_t = datestr(now);
syspara.na=naa;
syspara.iter=0;
syspara.ptRad=R0;
syspara.SDX=R0/20;syspara.SDY=R0/20;syspara.SDZ=R0/20;
syspara.NSX=40; syspara.NSY=40; syspara.NSZ=40;
syspara.mdata{:} = total_data{:}(:,1:4:end);
myplot = @(x,optimValues,state)plotGeometry(x,optimValues,state,fullX,fitIndex,syspara);     
options=optimset('LargeScale','on','TolFun',1e-10,'TolX',1e-10,'MaxIter',100,'MaxFunEvals',length(fitX)*2000,'Display','iter','Algorithm','sqp','PlotFcns',{@optimplotx,...
@optimplotfval,myplot});
[x,fval,exitflag,output]=fmincon(@(x)forward_proj_PET_multiang_new_v2(x, syspara, fitIndex ,fullX), fitX, [], [], [], [],lb, ub, nonlcon,options);
fitX=x;
%% display the result
fullX(fitIndex==1)=x;

forward_proj_PET_multiang_display_new_v2(fullX,fitIndex, syspara,total_data,0);
if iSave
    filename=[resultfolder,'/results0iter',datestr(now,30)];
    save(filename,'fitIndex','fullX');
    saveas(gcf,filename,'fig');
end

total_data = ExcludeChanceCoin(syspara, total_data, fullX);
end_t(1,:) = datestr(now);
for i = 1:subset
    subset_data{:,i} = total_data{:}(:,i:subset:end);
end

tolerance = 1e-10;
%% second fitting using fmincon
for it=1:subset
        syspara.na=naa;
        syspara.iter=it;
        syspara.ptRad=R(it);
        syspara.SDX=R(it)/20;syspara.SDY=R(it)/20;syspara.SDZ=R(it)/20;
        syspara.NSX=40; syspara.NSY=40; syspara.NSZ=40;
        syspara.mdata{:}=subset_data{:,it};
        myplot = @(x,optimValues,state)plotGeometry(x,optimValues,state,fullX,fitIndex,syspara);
        %tolerance = tolerance*sqrt(0.1)
        %% using standard fmincon
        options=optimset('LargeScale','on','TolFun',1e-10,'TolX',tolerance,'MaxIter',100,'MaxFunEvals',length(fitX)*2000,'Display','iter','Algorithm','sqp','PlotFcns',{@optimplotx,...
@optimplotfval,myplot});
        [x,fval,exitflag,output]=fmincon(@(x)forward_proj_PET_multiang_new_v2(x, syspara, fitIndex ,fullX), fitX, [], [], [], [],lb, ub, nonlcon,options);   
        fitX=x;
        %% display the result
        fullX(fitIndex==1)=x;
        forward_proj_PET_multiang_display_new_v2(fullX,fitIndex, syspara,total_data,iSave);
        if iSave
            filename=[resultfolder,'/results',num2str(it),'iter',datestr(now,30)];
            save(filename,'fitIndex','fullX');
            saveas(gcf,filename,'fig');
        end
        end_t(it+1,:)=datestr(now)
        
        
        pause(1);
end
[start_t; end_t]
fit_data_multiang_fmincon=x;
%%% initial fitting using fmincon (end)

%parpool close




