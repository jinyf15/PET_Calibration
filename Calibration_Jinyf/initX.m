function [fitX,fitIndex,fullX] = initX
fitX = zeros(1,43);
fitIndex = zeros(1,43);

fitX(1) = -2; %radius of rotation
fitIndex(1)=1;

fitX(2:4) = [0,-50,0]; % PET1 center position
fitX(5:7) = [0,-pi/2,-pi/2]; %Eular rotation angle
fitX(8:11) = [0.8,0.8,0.8,0.8]; % width of gap:x<0, x>0, y<0, y>0

fitX(12:14) = [0,50,0]; % PET2
fitX(15:17) = [0,-pi/2,pi/2];
fitX(18:21) = [0.8,0.8,0.8,0.8];

fitX(22:24) = [-50,0,0]; % PET3
fitX(25:27) = [0,-pi/2,pi];
fitX(28:31) = [0.8,0.8,0.8,0.8];

fitX(32:34) = [50,0,0]; % PET4
fitX(35:37) = [0,-pi/2,0];
fitX(38:41) = [0.8,0.8,0.8,0.8];
fitIndex(2:41) = 1;

fitX(42:43) = [2,-5];
fitIndex(42:43) = 1;
fullX = fitX;
fitX(fitIndex==0)=[];