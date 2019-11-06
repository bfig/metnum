% Funcion g1(x) = c/x^p
c = 0.96750; % con outlier
p = 2; % con outlier
% [c p] = [2 0.94613] sin outlier

% Funcion g2(x) = d/x^q
d = 1.9949;
q = 3;

% EDO:
% y' = g1(x)*y = g2(x)
% y(0.5) = 0
% x in [0.5 .. 2]

a = 0.5;
b = 2;
y0 = 0;
% calculo paso h dentro de region de estabilidad.
h = 0.015 ; % VALOR PROVISORIO %
N = ceil((b-a)/h);
% Recalculo h
h = (b-a)/N;

x_exacta = linspace(a,b,100);
% sol analitica con p = 2 y q = 3:
k = (-d/c)*(1/c + 2)*exp(-2*c); % obtenido con la cond inicial.
y_exacta = d/c^2 + d./(c*x_exacta) + k*exp(c./x_exacta);

% Euler hacia atras:
x = linspace(a,b,N);
y_EA = zeros(length(x),1);
y_EA(1) = y0;
for i = 2:length(x);
    y_EA(i) = (y_EA(i-1) + d*h/x(i)^3)/(1 + c*h/x(i)^2);
    %y_EA(i) = ( x(i)^2*y_EA(i-1) + d*h/x(i) )/(x(i)^2 + c*h); % otra forma
end

% Euler hacia adelante:
% usa el mismo x que euler hacia atras
y_E = zeros(length(x),1);
y_E(1) = y0;
for i = 2:length(x);
    y_E(i) = y_E(i-1) + h*(-c*y_E(i-1)/x(i-1)^2 + d/x(i-1)^3);
end

% Runge-Kuta con ode solver ode45()
function ans = ecdif(x,y)
    c = 0.96750;
    d = 1.9949;
    ans = -c*y/x^2 + d/x^3;
end

[x_RK, y_RK] = ode45(@ecdif, [a b], y0);

% Graficas comparativas de Euler hacia atras, Euler hacia adelante, R-K y sol analitica:
figure(1)
%Sol analitica
plot(x_exacta, y_exacta, "b")
hold on
% Euler hacia adelante
plot(x, y_E, "r")
% Euler hacia atras
plot(x, y_EA, "m")
% Runge-Kuta
plot(x_RK, y_RK, "k")

axis ([a b])
title ("Sol de la EDO");
xlabel ("x");
ylabel ("y");
hold off
legend_text = legend ("Sol Analitica", "Euler hacia adelante", "Euler hacia atras", "Runge-Kuta");
legend (legend_text, "location", "southeast");
% -----------------------------------------------------------------------------------