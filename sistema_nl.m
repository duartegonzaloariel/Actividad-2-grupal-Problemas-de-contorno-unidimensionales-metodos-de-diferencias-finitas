function [F, a, b, c] = sistema_nl(w, x, h)
% Evalúa F(w) y las diagonales de la jacobiana del problema
% y'' - x y y' + 3y^2 = 6x, 3y(0)-2y'(0)=0, 2y(1)-y'(1)/3=1
n = length(w);            % n = N+1 incógnitas w_0..w_N
F = zeros(n,1); a = zeros(n,1); b = zeros(n,1); c = zeros(n,1);

% Nodo x_0 = 0 (nodo ficticio eliminado)
F(1) = (2 + 3*h)*w(1) - 2*w(2) - 3*h^2*w(1)^2;
a(1) = 2 + 3*h - 6*h^2*w(1);
c(1) = -2;

% Nodos interiores
for k = 2:n-1
    yp   = (w(k+1) - w(k-1))/(2*h);
    F(k) = -w(k-1) + 2*w(k) - w(k+1) + h^2*(x(k)*w(k)*yp - 3*w(k)^2 + 6*x(k));
    b(k) = -1 - h/2*x(k)*w(k);
    a(k) =  2 + h^2*(x(k)*yp - 6*w(k));
    c(k) = -1 + h/2*x(k)*w(k);
end

% Nodo x_N = 1 (nodo ficticio eliminado)
F(n) = -2*w(n-1) + (2 - 12*h)*w(n) + 6*h + h^2*(3*w(n)^2 - 3*w(n) + 6);
b(n) = -2;
a(n) = 2 - 12*h + h^2*(6*w(n) - 3);
end