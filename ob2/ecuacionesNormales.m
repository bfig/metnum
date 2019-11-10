% Minimos cuadrados a los datos transformados de g1 y g2. PMCL resolviendo directamente las ecs normales.
X = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]';
%Y1
Y1 = [3.89,2.75,2.01,1.61,1.21,0.89,0.69,0.63,0.44,0.42,0.70,0.32,0.40,0.26,0.32,0.25]';

%Y2
%Y1 = [15.96,9.45,5.75,3.82,2.89,2.17,1.22,1.05,0.86,0.63,0.69,0.40,0.44,0.29,0.43,0.20]';

%   n cantidad de muestras 
n = length(X);

A = zeros(n,2);
Xtransf = log(X);
Ytransf = log(Y1);
C=zeros(2,1);
A = [Xtransf ones(n,1)];

C=(A' * A) \ (A' * Ytransf);

p= - round(C(1))
c= exp(C(2))

% Valor minimizado: ||AX - Y||^2
error_g = sum((c*X.^(-p) - Y1).^2);

%---------------------------------------
% Graficas para correccion del modelo: %
%---------------------------------------
% Funcion g1
% Grafica modelo lineal con transformacion logaritmica de los datos
figure (1)
plot(Xtransf, Ytransf, 'bs')
hold on
vect_Xcv = linspace(Xtransf(1),Xtransf(n),50);
plot(vect_Xcv, C(1)*vect_Xcv + C(2), 'r')
title ("Funcion g - lineal");
xlabel ("x");
ylabel ("y");
legend("datos transformados", "Sol Ecs Normales");

% Grafica sin transformacion logaritmica de los datos
figure (2)
plot(X, Y1, 'bs')
hold on
vect_X = linspace(X(1),X(n),50);
plot(vect_X, c*vect_X.^(-p), 'r')
title ("Funcion g - original");
xlabel ("x");
ylabel ("y");
legend("datos originales", "Sol Ecs Normales");
%----------------------------
