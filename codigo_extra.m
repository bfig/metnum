clear

% Funcion f(x,y)
function z = f_11(x,y)
  z = x^2 + 3*y^2 - (2/3)*x*y + 2*x - 4*y + 2*exp(x+y);
end

% Funcion f(X) usando Q y b
function z = f_12(X)
  % X vector columna (x ; y)
  Q = [1 -1/3 ; -1/3 3];
  b = [-1 ; 2];
  z = (X'*Q*X - 2*b'*X) + 2*exp(X(1)+X(2));
end

% Funcion F(X) gradiente de f(X)
function Z = F(X)
  % X vector columna (x ; y)
  %Z = zeros(2,1); % vector columna
  Q = [1 -1/3 ; -1/3 3];
  b = [-1 ; 2];
  Z = (Q*X - b) + [1 ; 1]*exp(X(1)+X(2));
end

% Funcion JF(X) jacobiana de F(X)
function M = JF(X)
  M = zeros(2,2); % matriz cuadrada 2x2
  Q = [1 -1/3 ; -1/3 3];
  M = Q + [1 1 ; 1 1]*exp(X(1)+X(2))
end

% Funcion SL(X_k, s_k) solucion al sist lineal JF(X_k)*s_k = -F(X_k)
function Z = SL(X_k, s_k)
  Z = JF(X_k)*s_k + F(X_k);
end

function Q = Q_n(n)
  res = zeros(n,n);
  for i = 1:(n-1)
    res(i,i) = 2*i-1;
    res(i,i+1) = (-1)^i/(3*i);
    res(i+1,i) = res(i,i+1); %Fuerzo simetria
  end
  res(n,n) = 2*n-1;
  Q = res;
end

% rank 1 Cholesky updating
% Funciones auxiliares

% No vale la pena, se incluye su codigo directamente en CholeskyModifyA.
function [c, s, L_ans] = Compute(L_dato, a, i)
  w = sqrt(L_dato(i,i)^2 + a(i)^2);
  c(i) = w/L_dato(i,i);
  s(i) = a(i)/L_dato(i,i);
  L_ans(i,i) = w;
end

function [L_ans, a_ans] = Apply(c, s, L_dato, a_dato, i, j)
  L_ans(i,j) = (L_dato(i,j) + s(i)*a_dato(j))/c(i);
  a_ans(j) = c(i)*a_dato(j) - s(i)*L_ans(i,j);
end

function [L, a] = CholeskyModifyA(L_dato, a_dato)
  L = L_dato;
  a = a_dato;
  n = length(a);
  c = zeros(n);
  s = zeros(n);
  for i = 1:n
    % Compute
    i
    w = sqrt(L(i,i)^2 + a(i)^2);
    c(i) = w/L(i,i);
    s(i) = a(i)/L(i,i);
    L(i,i) = w
    input('Press any key for next...');
    % fin Compute
    for j = 1:(i-1)
      % Apply
      texto = sprintf('(i,j) = (%d,%d)',i,j);
      disp(texto)
      L(j,i) = (L(j,i) + s(j)*a(i))/c(j)
      a(i) = c(j)*a(i) - s(j)*L(j,i)
      input('Press any key for next...');
      %fin Apply
    end % for
  end % for
end % CholeskyModifyA

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

%---------------------------%
% Programa Principal
%---------------------------%

% [minimo, fval, info] = fsolve(@F, [1;2]);
% minimo
% fval
% info

n = input('Valor de n: ');
a = ones(n,1)
Q = Q_n(n)
L = chol(Q)
disp('Calculo L para Q + aat')
[L_nueva,a_nueva] = CholeskyModifyA(L,a);