function [P]=polynomialfit(x,y,n)

% Inputs:
% x is the Nx1 abscissa vector
% y is the Nx1 ordinate vector
% n is the polynomial degree

% Output:
% P is a nx1 vector whose components are the coefficients of the n-degree
% fitted polynomial, P(n+1)+P(n-2)x +...+ P(2)x^n-1+P(1)x^n.



x = x(:); y = y(:);

A = zeros(n+1,n+1);
Y = zeros(n+1,1);

for i=n+1:-1:1
    for j=n+1:-1:1
       k=((n+1)-j)+(n+1-i);
       A(i,j) = sum(x.^k);
    end    
    Y(i) = sum(y.*(x.^(n+1-i)));
end
P = A\Y; % polynomial coefficients

end