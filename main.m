clear; clc; close all;
N = 80; tol = 1e-10; maxit = 50;

%% ---------------- APARTADO 2: Newton + Crout ----------------
[x, wC, itC, histC, lmin] = newton_df(N, 'crout', tol, maxit);
fprintf('=== APARTADO 2 (Crout) ===\n');
fprintf('Iteraciones: %d\n', itC);
for k = 1:itC
    fprintf('  k = %d   ||v||_inf = %.4e\n', k, histC(k));
end
fprintf('Pivote minimo |l_ii|: %.6f\n', lmin);

idx = N-8:N+1;             % nodos x_71 ... x_80
fprintf('\n  i      x_i        w_i (Crout)\n');
for i = idx
    fprintf('%3d   %.4f   %.10f\n', i-1, x(i), wC(i));
end

figure;
plot(x, wC, 'b-', 'LineWidth', 1.5); grid on;
xlabel('x'); ylabel('y(x)');
title('Solucion por diferencias finitas (N = 80, Crout)');
exportgraphics(gcf, 'fig_ap2_solucion.pdf', 'ContentType', 'vector');

%% ---------------- APARTADO 3: Gauss y bvp5c ----------------
[~, wG, itG] = newton_df(N, 'gauss', tol, maxit);

odefun = @(t, y) [y(2); t*y(1)*y(2) - 3*y(1)^2 + 6*t];
bcfun  = @(ya, yb) [3*ya(1) - 2*ya(2); 2*yb(1) - yb(2)/3 - 1];
solinit = bvpinit(x', @(t) [2/9 + t/3; 1/3]);
opts = bvpset('RelTol', 1e-10, 'AbsTol', 1e-12);
sol = bvp5c(odefun, bcfun, solinit, opts);
yB = deval(sol, x'); yB = yB(1,:)';

% Tiempos (timeit repite la ejecución y devuelve la mediana)
tC = timeit(@() newton_df(N, 'crout', tol, maxit));
tG = timeit(@() newton_df(N, 'gauss', tol, maxit));
tS = timeit(@() newton_df(N, 'sparse', tol, maxit));
tB = timeit(@() bvp5c(odefun, bcfun, solinit, opts));

fprintf('\n=== APARTADO 3 ===\n');
fprintf('Iteraciones Gauss: %d\n', itG);
fprintf('max|Crout - Gauss| = %.4e\n', max(abs(wC - wG)));
fprintf('max|Crout - bvp5c| = %.4e\n', max(abs(wC - yB)));
fprintf('Error max Crout respecto de x^3 = %.4e\n', max(abs(wC - x.^3)));
fprintf('Error max Gauss respecto de x^3 = %.4e\n', max(abs(wG - x.^3)));
fprintf('Error max bvp5c respecto de x^3 = %.4e\n', max(abs(yB - x.^3)));
fprintf('Nodos finales malla bvp5c: %d\n', numel(sol.x));
fprintf('Tiempo Crout: %.4e s\n', tC);
fprintf('Tiempo Gauss: %.4e s\n', tG);
fprintf('Tiempo Gauss sparse: %.4e s\n', tS);
fprintf('Tiempo bvp5c: %.4e s\n', tB);

fprintf('\n  i      x_i       Crout          Gauss          bvp5c        |Crout-bvp5c|\n');
for i = idx
    fprintf('%3d   %.4f   %.10f   %.10f   %.10f   %.3e\n', ...
        i-1, x(i), wC(i), wG(i), yB(i), abs(wC(i) - yB(i)));
end

figure;
subplot(1,2,1);
plot(x, wC, 'b-', 'LineWidth', 1.5); hold on;
plot(x, yB, 'r--', 'LineWidth', 1.5); grid on;
xlabel('x'); ylabel('y(x)'); legend('Dif. finitas (Crout)', 'bvp5c', 'Location', 'best');
subplot(1,2,2);
semilogy(x, abs(wC - yB) + eps, 'k-', 'LineWidth', 1.2); grid on;
xlabel('x'); ylabel('|w_i - y_{bvp5c}(x_i)|');
exportgraphics(gcf, 'fig_ap3_bvp5c.pdf', 'ContentType', 'vector');

%% ------- Comprobación del orden 2 (referencia: bvp5c) -------
Ns = [10 20 40 80 160 320];
err = zeros(size(Ns));
for j = 1:numel(Ns)
    [xj, wj] = newton_df(Ns(j), 'crout', tol, maxit);
    yj = deval(sol, xj'); 
    err(j) = max(abs(wj - yj(1,:)'));
end
fprintf('\n   N      error max      orden\n');
for j = 1:numel(Ns)
    if j == 1
        fprintf('%4d   %.4e      --\n', Ns(j), err(j));
    else
        fprintf('%4d   %.4e   %.3f\n', Ns(j), err(j), log2(err(j-1)/err(j)));
    end
end

figure;
loglog(1./Ns, err, 'bo-', 'LineWidth', 1.5); hold on;
loglog(1./Ns, err(end)*(Ns(end)./Ns).^2, 'k:', 'LineWidth', 1.2); grid on;
xlabel('h'); ylabel('Error maximo'); legend('Dif. finitas', 'O(h^2)', 'Location', 'best');
exportgraphics(gcf, 'fig_ap3_orden.pdf', 'ContentType', 'vector');

%% ------- Tiempos para mallas mas finas (Crout vs Gauss) -------
fprintf('\n    N      Crout        Gauss llena   Gauss sparse\n');
for NN = [80 320 1280 2560]
    t1 = timeit(@() newton_df(NN, 'crout',  tol, maxit));
    t2 = timeit(@() newton_df(NN, 'gauss',  tol, maxit));
    t3 = timeit(@() newton_df(NN, 'sparse', tol, maxit));
    fprintf('%5d   %.4e   %.4e   %.4e\n', NN, t1, t2, t3);
end

%% ---------------- APARTADOS 4 y 5: problema lineal ----------------
% y'' - x y' - 3y = 0 con las mismas condiciones de contorno
[xL, wL, lminL] = df_lineal(N);
[~, yS, ypS, s, Y2fin] = disparo_lineal(N);

fprintf('\n=== APARTADO 5 (problema lineal) ===\n');
fprintf('Disparo: y2(1) = %.10f, y2''(1) = %.10f, s = %.10f\n', Y2fin(1), Y2fin(2), s);
fprintf('CC disparo: 3y(0)-2y''(0) = %.2e, 2y(1)-y''(1)/3 = %.12f\n', ...
        3*yS(1) - 2*ypS(1), 2*yS(end) - ypS(end)/3);
fprintf('Pivote minimo Crout |l_ii|: %.6f\n', lminL);
fprintf('y(0): DF = %.10f, disparo = %.10f\n', wL(1), yS(1));

fprintf('\n  i      x_i       w_i (DF)        y_i (disparo)    |dif|\n');
for i = idx
    fprintf('%3d   %.4f   %.10f   %.10f   %.3e\n', ...
        i-1, xL(i), wL(i), yS(i), abs(wL(i) - yS(i)));
end
[dmax, imax] = max(abs(wL - yS));
fprintf('Diferencia maxima DF-disparo = %.4e en x = %.4f\n', dmax, xL(imax));

figure;
subplot(1,2,1);
plot(xL, wL, 'b-', 'LineWidth', 1.5); hold on;
plot(xL, yS, 'r--', 'LineWidth', 1.5); grid on;
xlabel('x'); ylabel('y(x)');
legend('Dif. finitas (Crout)', 'Disparo lineal (RK4)', 'Location', 'best');
subplot(1,2,2);
semilogy(xL, abs(wL - yS) + eps, 'k-', 'LineWidth', 1.2); grid on;
xlabel('x'); ylabel('|w_i - y_i^{disparo}|');
exportgraphics(gcf, 'fig_ap5_lineal.pdf', 'ContentType', 'vector');

% Orden de DF en el problema lineal (referencia: disparo con paso h/8)
errL = zeros(size(Ns));
for j = 1:numel(Ns)
    [~, wj] = df_lineal(Ns(j));
    [~, yr] = disparo_lineal(8*Ns(j));
    errL(j) = max(abs(wj - yr(1:8:end)));
end
fprintf('\n   N      error max      orden   (problema lineal)\n');
for j = 1:numel(Ns)
    if j == 1
        fprintf('%4d   %.4e      --\n', Ns(j), errL(j));
    else
        fprintf('%4d   %.4e   %.3f\n', Ns(j), errL(j), log2(errL(j-1)/errL(j)));
    end
end