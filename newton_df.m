function [x, w, iter, hist, lmin] = newton_df(N, metodo, tol, maxit)
% Diferencias finitas + Newton. metodo = 'crout' o 'gauss'
h = 1/N;
x = (0:N)'*h;
w = 2/9 + x/3;            % recta que cumple ambas condiciones de contorno
hist = zeros(maxit,1);
lmin = Inf;

for iter = 1:maxit
    [F, a, b, c] = sistema_nl(w, x, h);
    switch metodo
        case 'crout'
            [v, lm] = crout_tri(a, b, c, -F);
            lmin = min(lmin, lm);
        case 'gauss'
            J = diag(a) + diag(b(2:end), -1) + diag(c(1:end-1), 1);
            v = J \ (-F);
    end
    w = w + v;
    hist(iter) = norm(v, inf);
    if hist(iter) < tol
        hist = hist(1:iter);
        return
    end
end
warning('Newton no ha convergido en %d iteraciones', maxit);
end