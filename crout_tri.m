function [x, lmin] = crout_tri(a, b, c, d)
% Resuelve A x = d con A tridiagonal mediante la factorización de Crout
% a: diagonal principal (n), b: subdiagonal (b(1) no se usa),
% c: superdiagonal (c(n) no se usa), d: término independiente
n = length(a);
l = zeros(n,1); u = zeros(n,1); z = zeros(n,1); x = zeros(n,1);

% Factorización A = L U y resolución de L z = d
l(1) = a(1);
u(1) = c(1)/l(1);
z(1) = d(1)/l(1);
for i = 2:n-1
    l(i) = a(i) - b(i)*u(i-1);
    u(i) = c(i)/l(i);
    z(i) = (d(i) - b(i)*z(i-1))/l(i);
end
l(n) = a(n) - b(n)*u(n-1);
z(n) = (d(n) - b(n)*z(n-1))/l(n);

% Resolución de U x = z
x(n) = z(n);
for i = n-1:-1:1
    x(i) = z(i) - u(i)*x(i+1);
end

lmin = min(abs(l));   % pivote mínimo (comprobación)
end