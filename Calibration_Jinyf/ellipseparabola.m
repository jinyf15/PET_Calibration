function [c, ceq] = ellipseparabola(x,fullX,fitIndex)
fullX(fitIndex==1) = x;
x = fullX;

det_c = zeros(3,16);
PET_c = zeros(3,4);
eulerAng_PET = zeros(3,4);
eulerAng_det = zeros(3,16);
Mtr_PET = zeros(3,3,4);
Mtr_det = zeros(3,3,16);

for i = 0:3
    PET_c(:,i+1) = x(1+i*30:3+i*30);
    eulerAng_PET(:,i+1) = x(4+i*30:6+i*30);
    det_c(:,i*4+1:i*4+4) = reshape(x(7+i*30:18+i*30),3,4);
    eulerAng_det(:,i*4+1:i*4+4) = reshape(x(19+i*30:30+i*30),3,4);
	Mtr_PET(:,:,i+1) = euler2TranMatrix(eulerAng_PET(1,i+1),eulerAng_PET(2,i+1),eulerAng_PET(3,i+1));
    for j = 1:4
        Mtr_det(:,:,j+i*4) = euler2TranMatrix(eulerAng_det(1,i*4+j),eulerAng_det(2,i*4+j),eulerAng_det(3,i*4+j));
    end
end
NormalVector = zeros(3,4);
for i = 1:4
    NormalVector(:,i) = Mtr_PET(:,:,i) * [0;0;1];
end
costheta12 = VectorAngle(NormalVector(:,1),NormalVector(:,2));
costheta34 = VectorAngle(NormalVector(:,3),NormalVector(:,4));
cosVertical = VectorAngle(PET_c(:,2) - PET_c(:,1),PET_c(:,4) - PET_c(:,3));
c(1) = costheta12 + 0.995; %parallel
c(2) = costheta34 + 0.995; %parallel
c(3) = abs(cosVertical) - 0.005;  %vertical

t = (NormalVector(:,1)'*PET_c(:,1)-NormalVector(:,1)'*PET_c(:,2))/(NormalVector(:,1)'*NormalVector(:,1));
NewC = PET_c(:,2) + NormalVector(:,1)*t;
c(4) = sqrt((NewC-PET_c(:,1))'*(NewC-PET_c(:,1))) - 2;
t = (NormalVector(:,2)'*PET_c(:,2)-NormalVector(:,2)'*PET_c(:,1))/(NormalVector(:,2)'*NormalVector(:,2));
NewC = PET_c(:,1) + NormalVector(:,2)*t;
c(5) = sqrt((NewC-PET_c(:,2))'*(NewC-PET_c(:,2))) - 2;

t = (NormalVector(:,3)'*PET_c(:,3)-NormalVector(:,3)'*PET_c(:,4))/(NormalVector(:,3)'*NormalVector(:,3));
NewC = PET_c(:,4) + NormalVector(:,3)*t;
c(6) = sqrt((NewC-PET_c(:,3))'*(NewC-PET_c(:,3))) - 2;
t = (NormalVector(:,4)'*PET_c(:,4)-NormalVector(:,4)'*PET_c(:,3))/(NormalVector(:,4)'*NormalVector(:,4));
NewC = PET_c(:,3) + NormalVector(:,4)*t;
c(7) = sqrt((NewC-PET_c(:,4))'*(NewC-PET_c(:,4))) - 2;
 c = c*1e9;
ceq = [];

end

