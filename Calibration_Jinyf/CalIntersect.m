function [x,y,z] = CalIntersect(x1,x2,scan,index)

if scan == 'X'
    x = index;
    t = (x-x1(1))/(x2(1)-x1(1));
    y = t * (x2(2)-x1(2)) + x1(2);
    z = t * (x2(3)-x1(3)) + x1(3);
    x = x*2+50.5;
    y = round(y*2+50.5);
    z = round(z*2+50.5);
else
    y = index;
    t = (y-x1(2))/(x2(2)-x1(2));
    x = t * (x2(1)-x1(1)) + x1(1);
    z = t * (x2(3)-x1(3)) + x1(3);
    x = round(x*2+50.5);
    y = y*2+50.5;
    z = round(z*2+50.5);
end


end

