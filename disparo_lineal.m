function [x, y, yp, s, Y2fin] = disparo_lineal(N)
% Metodo de disparo lineal para
%   y'' = x y' + 3y,  3y(0)-2y'(0)=0,  2y(1)-y'(1)/3=1
% Con alpha_a=3, beta_a=-2, gamma_a=0 el primer PVI tiene solucion y1 = 0,
% asi que y = s*y2 con y2(0) = beta_a = -2, y2'(0) = -alpha_a = -3.
% El PVI se integra con Runge-Kutta de orden 4 y paso h = 1/N.
h = 1/N;
x = (0:N)'*h;
f = @(t, Y) [Y(2); t*Y(2) + 3*Y(1)];

Y = zeros(2, N+1);
Y(:,1) = [-2; -3];
for i = 1:N
    t  = x(i); Yi = Y(:,i);
    k1 = f(t,       Yi);
    k2 = f(t + h/2, Yi + h/2*k1);
    k3 = f(t + h/2, Yi + h/2*k2);
    k4 = f(t + h,   Yi + h*k3);
    Y(:,i+1) = Yi + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

Y2fin = Y(:,end);                          % [y2(1); y2'(1)]
s  = 1/(2*Y2fin(1) - Y2fin(2)/3);         % gamma_b/(alpha_b y2(1) + beta_b y2'(1))
y  = s*Y(1,:)';
yp = s*Y(2,:)';
end