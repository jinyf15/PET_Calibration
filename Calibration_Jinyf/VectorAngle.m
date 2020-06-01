function output = VectorAngle(a,b)
A = sqrt(a'*a);
B = sqrt(b'*b);
output = a'*b/(A*B);

end

