X = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]';
Y1 = [3.89,2.75,2.01,1.61,1.21,0.89,0.69,0.63,0.44,0.42,0.70,0.32,0.40,0.26,0.32,0.25]';
Y2 = [15.96,9.45,5.75,3.82,2.89,2.17,1.22,1.05,0.86,0.63,0.69,0.40,0.44,0.29,0.43,0.20]';

%   n cantidad de muestras 
n = length(X);

function C = PMCLQR(A,Y)
  % Resuelve el PMCL  min(|| AC - Y ||^2) usando descomp QR
  % Matriz A m x n ; vector Y m x 1
  % Calculo matriz Q reducida (1eras n filas de Q)
  % Q = [Q_r Q2]
  [m n] = size(A);
  Q_r = zeros(m,n);
  Q_r(:,1) = A(:,1)/norm(A(:,1));
  for i = 2:n;
    aux = A(:,i);
    for j = 1:(i-1);
      aux = aux - dot(A(:,i),Q_r(:,j))*Q_r(:,j);
    end
  Q_r(:,i) = aux/norm(aux);
  end
  % Calculo matriz R reducida (1eras n filas de R)
  % R = [   R_r  ]
  %     [0 0 .. 0]
  R_r = Q_r'*A;
  % Los elementos fuera del triangulo quedan en x10^-17 etc, entonces los hago cero.
  for i = 1:n;
    for j = 1:(i-1);
    R_r(i,j) = 0;
    end
  end
  % Resuelvo el sistema por sustitucion hacia atras
  ec_opt.UT = true;
  C = linsolve(R_r,Q_r'*Y,ec_opt);
end % PMCLQR

function PCk = PMCNL(PCinit,X,Y,iter=10)
  % Calcula la solucion al PMCNL con datos X e Y usando Gauss-Newton
  fun = @(PC) PC(2)*X.^(-PC(1)); % funcion f
  diffun = @(PC) [(-PC(2)*X.^(-PC(1)).*log(X)) (X.^(-PC(1)))]; % jacobiano de f
  PCk = PCinit;
  for i = 1:iter 
    Ak = diffun(PCk);
    Yk = Y - fun(PCk);
    % PCk = PCk + ols(Yk,Ak);
    PCk = PCk + PMCLQR(Ak,Yk);
  endfor
endfunction

% Calculo del error minimizado:
PC = PMCNL([1;1],X,Y1,10);
p = PC(1);
c = PC(2);
error_g1 = sum((c*X.^(-p) - Y1).^2);

QD = PMCNL([1;1],X,Y2,10);
q = QD(1);
d = QD(2);
error_g2 = sum((d*X.^(-q) - Y2).^2);

%---------------------------------------
% Graficas para correccion del modelo: %
%---------------------------------------
% Funcion g1
figure (2)
plot(X, Y1, 'bs')
hold on
vect_X = linspace(X(1),X(n),50);
plot(vect_X, c*vect_X.^(-p), 'r')
title ("Funcion g1");
xlabel ("x");
ylabel ("y");
legend("datos originales", "Sol Gauss-Newton");
%----------------------------
% Funcion g2
figure (4)
plot(X, Y2, 'bs')
hold on
% vect_X = linspace(X(1),X(n),50);
plot(vect_X, d*vect_X.^(-q), 'r')
title ("Funcion g2");
xlabel ("x");
ylabel ("y");
legend("datos originales", "Sol Gauss-Newton");