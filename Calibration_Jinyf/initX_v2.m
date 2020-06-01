function [fitX,fitIndex,fullX] = initX_v2
fitX = zeros(1,127);
fitIndex = ones(1,127);


fitX(1:3) = [0,-50,0]; % PET1 center position
fitX(4:6) = [0,-pi/2,-pi/2]; %Eular rotation angle
for i = 0:3
    fitX(7+i*30:9+i*30) = [-11.25,-11.25,0]; % width of gap:x<0, x>0, y<0, y>0
    fitX(10+i*30:12+i*30) = [11.25,-11.25,0];
    fitX(13+i*30:15+i*30) = [-11.25,11.25,0];
    fitX(16+i*30:18+i*30) = [11.25,11.25,0];
    fitX(19+i*30:30+i*30) = zeros(1,12);
end
fitX(31:33) = [0,50,0]; % PET2
fitX(34:36) = [0,-pi/2,pi/2];

fitX(61:63) = [-50,0,0]; % PET3
fitX(64:66) = [0,-pi/2,pi];

fitX(91:93) = [50,0,0]; % PET4
fitX(94:96) = [0,-pi/2,0];

fitX(121:123) = [-4.25,0.75,-6.25];
fitX(124:125) = [0,0];
fitX(126:127) = [0,0]; %theta and phi

fullX = fitX;
fitX(fitIndex==0)=[];