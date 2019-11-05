% Minimos cuadrados a los datos transformados de g1 y g2. PMCL usando descomposicion QR
clear all;
X = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]';
% Funcion g1
Y1 = [3.89,2.75,2.01,1.61,1.21,0.89,0.69,0.63,0.44,0.42,0.70,0.32,0.40,0.26,0.32,0.25]';
% Funcion g2
Y2 = [15.96,9.45,5.75,3.82,2.89,2.17,1.22,1.05,0.86,0.63,0.69,0.40,0.44,0.29,0.43,0.20]';

m = length(X);
n = 2; % modelo lineal

% Modelo de funcion g
% g = @(ci,p) @(x) ci/x^p;

X_cv = log(X);
Y1_cv = log(Y1);
Y2_cv = log(Y2);

A = [X_cv ones(m,1)];
% Descomposicion QR de A:
% [Q, R] = qr(A);
% Resulta mejor calcular Q_r y R_r directamente
% Q = [Q_r Q2]
Q_r = zeros(m,n);
Q_r(:,1) = A(:,1)/norm(A(:,1));
for i = 2:n;
    aux = A(:,i);
    for j = 1:(i-1);
        aux = aux - dot(A(:,i),Q_r(:,j))*Q_r(:,j);
    end
    Q_r(:,i) = aux/norm(aux);
end

% R = [   R_r  ]
%     [0 0 .. 0]
% R_r = R(1:n,:);% 1eras n filas de R
R_r = Q_r'*A;
% Los elementos fuera del triangulo quedan en x10^-17 etc, entonces los hago cero.
for i = 1:n;
    for j = 1:(i-1);
        R_r(i,j) = 0;
    end
end
%--------------
% Funcion g1: %
%--------------
% QtY1_1 = (Q'*Y1_cv)(1:n); % 1eros n elementos del vector Q'Y
QtY1_1 = Q_r'*Y1_cv;

% resuelvo el sist lineal R1X = (Q'Y)_1 con R1 matriz triangular superior, resolver con sustitucion hacia atras.
ec1_opt.UT = true;
Param_1 = linsolve(R_r,QtY1_1,ec1_opt);
% solucion es y = X(1).t + X(2)
% y = c.x^(-p) ; log(y) = -p.log(x) + log(c) ; y_cv = A.t + B
p = round(-Param_1(1));
c = exp(Param_1(2));

% Valor minimizado: ||AX - Y||^2
dist_g1 = sum((c*X.^(-p) - Y1).^2);

%--------------
% Funcion g2: %
%--------------
% QtY1_2 = (Q'*Y2_cv)(1:n); % 1eros n elementos del vector Q'Y
QtY1_2 = Q_r'*Y2_cv;

% resuelvo el sist lineal R1X = (Q'Y)_1 con R1 matriz triangular superior, resolver con sustitucion hacia atras.
ec2_opt.UT = true;
Param_2 = linsolve(R_r,QtY1_2,ec2_opt);
% solucion es y = X(1).t + X(2)
% y = d.x^(-q) ; log(y) = -q.log(x) + log(d) ; y_cv = A.t + B
q = round(-Param_2(1));
d = exp(Param_2(2));

% Valor minimizado: ||AX - Y||^2
dist_g2 = sum((d*X.^(-q) - Y2).^2);

%---------------------------------------
% Graficas para correccion del modelo: %
%---------------------------------------
% Funcion g1
% Grafica modelo lineal con transformacion logaritmica de los datos
figure (1)
plot(X_cv, Y1_cv, 'bs')
hold on
vect_Xcv = linspace(X_cv(1),X_cv(m),50);
plot(vect_Xcv, Param_1(1)*vect_Xcv + Param_1(2), 'r')
title ("Funcion g1 - lineal");
xlabel ("x");
ylabel ("y");
legend("datos transformados", "ajuste QR");

% Grafica sin transformacion logaritmica de los datos
figure (2)
plot(X, Y1, 'bs')
hold on
vect_X = linspace(X(1),X(m),50);
plot(vect_X, c*vect_X.^(-p), 'r')
title ("Funcion g1 - original");
xlabel ("x");
ylabel ("y");
legend("datos originales", "ajuste QR");
%----------------------------
% Funcion g2
% Grafica modelo lineal con transformacion logaritmica de los datos
figure (3)
plot(X_cv, Y2_cv, 'bs')
hold on
% vect_Xcv = linspace(X_cv(1),X_cv(m),50);
plot(vect_Xcv, Param_2(1)*vect_Xcv + Param_2(2), 'r')
title ("Funcion g2 - lineal");
xlabel ("x");
ylabel ("y");
legend("datos transformados", "ajuste QR");

% Grafica sin transformacion logaritmica de los datos
figure (4)
plot(X, Y2, 'bs')
hold on
% vect_X = linspace(X(1),X(m),50);
plot(vect_X, d*vect_X.^(-q), 'r')
title ("Funcion g2 - original");
xlabel ("x");
ylabel ("y");
legend("datos originales", "ajuste QR");

% dado x, y = 3*x.^2 + 2*x + 1; => polyfit(x,y,2) = [3 2 1]
