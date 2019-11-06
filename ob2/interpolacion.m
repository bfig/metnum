% Funcion g1(x) = c/x^p
c = 0.96750; % con outlier
p = 2; % con outlier
% [c p] = [2 0.94613] sin outlier

% Funcion g2(x) = d/x^q
d = 1.9949;
q = 3;

a = 0.5;
b = 2;
y0 = 0;

x_exacta = linspace(a,b,100);
% sol analitica con p = 2 y q = 3:
k = (-d/c)*(1/c + 2)*exp(-2*c); % obtenido con la cond inicial.
y_exacta = d/c^2 + d./(c*x_exacta) + k*exp(c./x_exacta);

% Runge-Kuta con ode solver ode45()
function ans = ecdif(x,y)
    c = 0.96750;
    d = 1.9949;
    ans = -c*y/x^2 + d/x^3;
end
[x_RK, y_RK] = ode45(@ecdif, [a b], y0);

% interpolacion por splines cubicas
sp = linspace(a,b,100);
spp = spline(x_RK,y_RK,sp);


% Graficas comparativas de Euler hacia atras, Euler hacia adelante, R-K y sol analitica:
figure(1)
%Sol analitica
grosor = 1.5;
plot(x_exacta, y_exacta, "b",'LineWidth',grosor)
hold on
% Runge-Kuta
plot(x_RK, y_RK, "k")
% interpolacion por splines cubicas
plot(sp,spp,"--r",'LineWidth',grosor);

axis ([a b])
title ("Sol de la EDO");
xlabel ("x");
ylabel ("y");
hold off
legend_text = legend ("Sol Analitica", "Runge-Kuta", "Spline Cubica");
legend (legend_text, "location", "southeast");