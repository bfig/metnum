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

% Grilla original
X = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0];
% Puntos intermedios de la grilla
X_medios = linspace(0.55,1.95,15);

% sol analitica con p = 2 y q = 3:
k = (-d/c)*(1/c + 2)*exp(-2*c); % obtenido con la cond inicial.
y_exacta = d/c^2 + d./(c*X_medios) + k*exp(c./X_medios);

% Runge-Kuta con ode solver ode45()
function ans = ecdif(x,y)
    c = 0.96750;
    d = 1.9949;
    ans = -c*y/x^2 + d/x^3;
end
[x_RK, y_RK] = ode45(@ecdif, [a b], y0);

% interpolacion por splines cubicas
spp = spline(x_RK,y_RK,X_medios);

% Graficas comparativas de Euler hacia atras, Euler hacia adelante, R-K y sol analitica:
figure(1)
%Sol analitica
grosor = 1.5;
plot(X_medios, y_exacta, "b",'LineWidth',grosor)
hold on
% Runge-Kuta
plot(x_RK, y_RK, "k")
% interpolacion por splines cubicas
plot(X_medios,spp,"--r",'LineWidth',grosor);

axis ([a b])
%title ("Sol de la EDO");
xlabel ("x");
ylabel ("y");
hold off
legend_text = legend ("Sol Analitica", "Lineal a trozos", "Spline Cubica");
legend (legend_text, "location", "southeast");

% Comparacion de interpolantes en puntos intermedios de la grilla original:
error_spp = log((spp - y_exacta).^2);
figure(2)
plot(X_medios,error_spp);
xlabel ("x");
ylabel ("Error^2");
legend_text = legend ("Error splines cubicos en escala logaritmica");
%legend (legend_text, "location", "southeast");