function output = check_NG(temp)

tmp = zeros(1,8);
for i = size(temp,1)-1:-1:1
    if temp(i,3)==temp(end,3)
        tmp = temp(end,:);
        tmp(2) = temp(i,2);
        temp(i+2:end,:) = temp(i+1:end-1,:);
        temp(i+1,:) = tmp;
        break;
    end
end
output = temp;

