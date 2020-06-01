function [outputArg1,outputArg2] = initDetPixelPosi(DN_DET, DDX_DET, DDY_DET)
%iniDetDixelPosi: init positions of  center of one-layer pixels
%   Detailed explanation goes here

X = zeros(max(DN_DET)^2,4);
Y = zeros(max(DN_DET)^2,4);
for k = 1:4
    x_ind=-(DN_DET(k)-1)/2:(DN_DET(k)-1)/2;
    x_ind=repmat(x_ind,1,DN_DET(k)) * DDX_DET(k);
    X(1:DN_DET(k)^2,k) = x_ind;
    y_ind=-(DN_DET(k)-1)/2:(DN_DET(k)-1)/2;
    y_ind=repmat(y_ind,DN_DET(k),1);
    y_ind=reshape(y_ind,1,DN_DET(k)*DN_DET(k)) * DDY_DET(k);
    Y(1:DN_DET(k)^2,k) = y_ind;
end


outputArg1 = X;
outputArg2 = Y;
end

