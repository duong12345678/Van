#include<stdio.h>
#include<stdlib.h>
#include<math.h>

function [x, k]=gs_hang(A, B, e)
%Phuong phap Gauss-Seidel cho ma tran cheo hang
%A, B la cac ma tran dau vao thoa man Ax=B
if nargin == 2
    e = 1e-6;
end
[m, n] = size(A);
[p, ~] = size(B);
if m ~= n
    fprintf("A khong la ma tran vuong\n");
elseif m ~= p
    fprintf("so chieu cua vector B khac voi co cua A\n");
elseif ktra(A, m) == 0
    fprintf("A khong la ma tran cheo troi hang\n");
else

    y = diag(A);
    A = -A;
    
    for i = 1:1:m
        A(i, i) = 0;
        A(i,:) = A(i,:)/y(i);
        B(i) = B(i)/y(i);
    end
    
    lambda = 0;
    for i = 1:1:m
        p = 0;
        q = 0;
        for j = i:1:m
            p = p + abs(A(i, j));
        end
        for j = 1:1:i-1
            q = q + abs(A(i, j));
        end
        if p/(1-q) > lambda
            lambda = p/(1-q);
        end
    end
    
    
    x = zeros(1,m)';
    tong=0;
    for i=1:1:m
        c=A(i,:)*x+B(i);
        if tong < abs(c-x(i))
            tong=abs(c-x(i));
        end
        x(i)=c;
    end
    u=tong;
    k=1;
    
    while u*lambda >= e*(1-lambda)
        k=k+1;
        tong=0;
        for i=1:1:m
            c=A(i,:)*x+B(i);
            if tong < abs(c-x(i))
                tong=abs(c-x(i));
            end
            x(i)=c;
        end
        u=tong;
    end
    
end

end


function [a]=ktra(A, m)
a=1;
for i=1:1:m
    tong=0;
    for j=1:1:m
        tong = tong + abs(A(i, j));
    end
    if 2*abs(A(i, i)) <= tong
        a=0;
        return;
    end
end
end

