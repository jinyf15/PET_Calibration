clc
close all

syspara.NANG=1;
syspara.DM_DET(1:4)=20;syspara.DN_DET(1:4)=64;
syspara.DDZ_DET(1:4)=0.5;syspara.DDX_DET(1:4)=2*20.9./syspara.DN_DET;syspara.DDY_DET(1:4)=2*20.9./syspara.DN_DET;
syspara.ptRad=1;
syspara.DANG=2*pi/syspara.NANG;
syspara.NDET_DET=syspara.DN_DET.*syspara.DN_DET;
syspara.SDX=0.1;syspara.SDY=0.1;syspara.SDZ=0.1;
syspara.NSX=20; syspara.NSY=20; syspara.NSZ=20;NIMG=syspara.NSX*syspara.NSY*syspara.NSZ;
syspara.mu_det(1:4)=0.58; syspara.fn=6;

% read experiment data
for na=0:syspara.NANG-1
    foldername = '..\data\';
    filename = 'CombinedAllEvents_100ns.txt';  
    syspara.mdata{na+1}=preprocessing(foldername,filename,syspara);
    
end

%system parameters
[fitX,fitIndex,fullX] = initX;
[A,b,Aeq,beq,lb,ub] = constrains(fullX,fitIndex);
[syspara.dxx0, syspara.dyy0] = initDetPixelPosi(syspara.DN_DET,syspara.DDX_DET,syspara.DDY_DET);

%%%%% dispaly the initial guess
naa=0; 
syspara.na=0;
syspara.iter=0;
% display the result
% load('100nsfitpara.mat');
% forward_proj_PET_multiang_display_new(fitX, syspara)
% fitX = x;
forward_proj_PET_multiang_display_new(fullX, syspara);

%%% initial fitting using fmincon (start)

options=optimset('LargeScale','on','TolFun',1,'TolX',1e-5,'MaxIter',50,'MaxFunEvals',length(fitX)*2000,'Display','iter','Algorithm','sqp', 'Useparallel', 'always','PlotFcns',{@optimplotx,...
   @optimplotfval});
%options=optimset('LargeScale','on','TolFun',1,'TolX',1e-5,'MaxIter',50,'MaxFunEvals',length(fitX)*2000,'Display','iter','PlotFcns',{@optimplotx,@optimplotfval});
%options = optimset('Display','iter','MaxIter',10, 'DiffMinChange', 1e-8, 'TolFun',1e-8, 'TolCon', 1e-8, 'Algorithm', 'sqp', 'Useparallel', 'always');%, 'Outputfcn', @outfun);
for it=2
        start_t=datestr(now)
        syspara.na=naa;
        syspara.iter=it;           
        %% using standard fmincon
%         [x,fval,exitflag,output]=fmincon(@(x)forward_proj_PET_multiang_new(x, syspara), fitX, [], [], [], [],lb, ub, [],options);
%         fitX=x;
        
        %% save the current work space          
        %ms = MultiStart('UseParallel', 'always');
        fwdproj_multiang=@(x)forward_proj_PET_multiang_new(x, syspara, fitIndex ,fullX);
        %options = optimset('Display','iter', 'Algorithm', 'sqp', 'Useparallel', 'always');
        problem=createOptimProblem('fmincon', 'objective', fwdproj_multiang, 'x0', fitX, 'lb', lb, 'ub', ub,'options', options); 
        gs = GlobalSearch;
        [x,f] = run(gs,problem);
        fitX=x;
        %% display the result
        fullX(fitIndex==1)=x;
        forward_proj_PET_multiang_display_new(fullX, syspara);
        filename=[foldername,'/results',num2str(it),'iter',datestr(now,30)];
        save(filename,'fitIndex','fullX');
        saveas(gcf,filename,'fig');
        end_t=datestr(now)
        %[start_t end_t]
        
        pause(1);
end
fit_data_multiang_fmincon=x;
%%% initial fitting using fmincon (end)

%parpool close




