clc
clear


x1 = [1 1]'
x2 = [3 4]'
x3 = [5 6]'

x = [x1 x2 x3]

L = [0 0 0;
    -1 1 0;
    0 -1 1]

n = length(L);
uk = zeros(length(x(:,1)),1);
u = uk
for i = 1:n
    uk = [0 0]';
    for j = 1:n
    uk = uk + L(i,j)*x(:,j);
    end
    if i == 1
        u = uk;
    else
        u = [u uk];
    end

end

u
