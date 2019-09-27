1; % Identifica al archivo como script

% Definicion de la funcion para R2 en formato matricial
function z = f(X)
  Q = [1 -1/3; -1/3 3];
  b = [ -1 ; 2 ] ;
	z = X' * Q * X - 2* b' * X + 2* e^(X(1)+X(2));
endfunction

% Funcion F(X) gradiente de f(X) (dividido 2) en R2
function z = df(X)
  Q = [1 -1/3; -1/3 3];
  b = [ -1 ; 2 ] ;
  z = Q * X - b + [1 ; 1] * e^(X(1)+X(2));
endfunction

function z = ej1_5()
	z = fsolve(@df,[1 ; 1]);
endfunction

% Definicion de la funcion a(X) en R2
function v = a(X)
  Q = [1 -1/3; -1/3 3];
  b = [ -1 ; 2 ] ;
  v = [1 ; 1] * sqrt(e^(X(1)+X(2)));
endfunction

% Funcion JF(X) jacobiana de F(X) en R2
function j = jF(X)
  Q = [1 -1/3; -1/3 3];
  b = [ -1 ; 2 ] ;
  a_p = a(X);
  j = Q + a_p * a_p';
endfunction

% Implementacion de Newton-Raphson
function sol = NR(X,iters = 10)
  tmp = X;
  for i = 1:iters
    j_k = jF(tmp);
    f_k = df(tmp);
    tmp = tmp - j_k \ f_k;
  endfor
  sol = tmp;
endfunction

% calculo del error absoluto
function error = ej2_2()
  res_1_5 = ej1_5();
  res_2_3 = NR([1;1]);
  error = abs(res_1_5 - res_2_3);
endfunction

% Calculo de la norma de F(z^(k)) para cada iteracion k
function reses = ej2_6()
    iters = 10
  x = zeros(iters,1);
  tmp = [1;1];
  for i = 1:iters;
    tmp = NR(tmp,1);
    x(i) = norm(df(tmp));
  endfor
  reses = x;
endfunction

% Calculo de Q para dimension n
function Q = Q_n(n)
  res = zeros(n,n);
  for i = 1:(n-1)
    res(i,i) = 2*i-1;
    res(i,i+1) = (-1)^i/(3*i);
    res(i+1,i) = res(i,i+1); % Fuerzo simetria
  end
  res(n,n) = 2*n-1;
  Q = res;
end

% Calculo de b para dimension n
function b = b_n(n)
  b = zeros(n,1);
  for k = 1:n
    b(k) = k*(-1)^k;
  end
end

% rank 1 Cholesky updating
% Calculo de matriz de descomp. de Cholesky para
% actualizacion de rango 1 de Q
% L factor de Cholesky Triagular Superior de Q
function [L, a] = CholeskyModifyB(L_dato, a_dato)
  L = L_dato;
  a = a_dato;
  n = length(a);
  c = zeros(n);
  s = zeros(n);
  for i = 1:n
    % Compute
    w = sqrt(L(i,i)^2 + a(i)^2);
    c(i) = w/L(i,i);
    s(i) = a(i)/L(i,i);
    L(i,i) = w;
    % fin Compute
    for j = (i+1):n
      % Apply
      L(i,j) = (L(i,j) + s(i)*a(j))/c(i);
      a(j) = c(i)*a(j) - s(i)*L(i,j);
      %fin Apply
    end % for
  end % for
end % CholeskyModifyB

% Calculo de matriz de descomp. de Cholesky para
% matriz tridiagonal simetrica
% Devuelve matriz triangular inferior
function L = cholesky_tridiagonal(Q)
  delta = Q(1,1); % delta_1 = alpha_1
  n = size(Q)(1);
  L = zeros(n,n);
  L(1,1) = sqrt(delta);
  for i = 2:n
    L(i,i-1) = Q(i,i-1)/sqrt(delta); % use delta_(i-1)
    delta = Q(i,i) - (Q(i,i-1)^2)/delta; % calculo delta_(i)
    L(i,i) = sqrt(delta);
  end
end

% Calculo de matriz de descomp. de Cholesky para
% matriz tridiagonal simetrica
% Devuelve matriz triangular superior
function L = cholesky_tridiagonal_2(Q)
  delta = Q(1,1); % delta_1 = alpha_1
  n = size(Q)(1);
  L = zeros(n,n);
  for i = 1:(n-1)
    L(i,i) = sqrt(delta);
    L(i,i+1) = Q(i,i+1)/sqrt(delta); % use delta_(i-1)
    delta = Q(i+1,i+1) - (Q(i,i+1)^2)/delta; % calculo delta_(i)
  end
  L(n,n) = sqrt(delta);
end

function res = az(z)
  res = (exp(sum(z)))*ones(length(z),1); % [1 1 ... 1]'*exp(sum(x1,x2,...,xn))
end

function z = NR2(sQ, iters=1)
  z_tmp = (ones(sQ,1)/sQ); % [1/n 1/n ... 1/n]'
  b = b_n(sQ);
  Q = Q_n(sQ);
  L = cholesky_tridiagonal_2(Q);
  ec1_opt.LT = true;
  ec2_opt.UT = true;
  for i = 1:iters
    a_temp = az(z_tmp);
    a = sqrt(a_temp);
    %[Lm,info] = cholupdate(L',a);
    Lm = CholeskyModifyB(L,a);
    %NR: F(z_new-z_tmp) = F(z_tmp) + dF(z_tmp) * (z_new-z_tmp)
    % igualando a 0
    % z_new = z_tmp - dF^-1(z_tmp)(F(z_tmp))
    % 
    %(Q+b)x = -F(z_tmp)
    % LL'x = -F(z_tmp)
    % Ly = -F(z_tmp)
    % L' z_newtmp = y
    Fz = Q*z_tmp - b + a_temp;
    y = linsolve(Lm',-Fz,ec1_opt);
    x = linsolve(Lm,y,ec2_opt);
    z_tmp = z_tmp + x;
  endfor
  z = z_tmp;
endfunction

%function res1,res2 = ej5_1_2()
%  enes = [1,2,3,4,5,10,20,30,40,50]*100
%  valores = zeros(size(enes),1)
%  normas = zeros(size(enes),1)
%  tmp = ones(1,size(Q)(1))/size(Q)(1)
%  for i = 1:size(enes);
%    %
%    tmp = NR(tmp,1);
%    x(i) = norm(df(tmp));
%  endfor
%  res1 = x
%endfunction

% Calculo de los tiempos de la ejecucion de NR2 y calculo de norma de F(z_NR)
function [norma_Fz,tiempos] = ej5_2_3()
  enes = [1 2 3 4 5 10 20 30 40 50]*100;
  norma_Fz = zeros(1,length(enes));
  tiempos = zeros(1,length(enes));
  for (i = 1:length(enes))
    tic();
    z = NR2(enes(i),20);
    tiempos(i) = toc();

    Q = Q_n(enes(i));
    b = b_n(enes(i));
    a_temp = az(z);
    norma_Fz(i) = norm(Q*z - b + a_temp);
  end
end

function tiempo = prueba_tiempo(enes)
    tic();
    NR2(enes,20);
    tiempo = toc();
end
