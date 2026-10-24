function [x, w, lmin] = df_lineal(N)
% Diferencias finitas de segundo orden para el problema lineal
%   y'' - x y' - 3y = 0,  3y(0)-2y'(0)=0,  2y(1)-y'(1)/3=1
% Nodos ficticios eliminados con las condiciones de Robin. Sistema
% tridiagonal A w = d resuelto con el algoritmo de Crout.
h = 1/N;
x = (0:N)'*h;
n = N + 1;

a = (2 + 3*h^2)*ones(n,1);   % diagonal principal
b = -1 - h/2*x;              % subdiagonal  (coef. de w_{i-1})
c = -1 + h/2*x;              % superdiagonal (coef. de w_{i+1})
d = zeros(n,1);

% Fila i = 0:  w_{-1} = w_1 - 3h w_0
a(1) = 2 + 3*h + 3*h^2;
c(1) = -2;

% Fila i = N:  w_{N+1} = w_{N-1} + 12h w_N - 6h
b(n) = -2;
a(n) = 2 - 12*h + 9*h^2;
d(n) = -6*h + 3*h^2;

[w, lmin] = crout_tri(a, b, c, d);
end