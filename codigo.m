1;

function z = f(X)
  Q = [1 -1/3; -1/3 3];
  b = [ -1 ; 2 ] ;
	z = X' * Q * X - 2* b' * X + 2* e^(X(1)+X(2));
endfunction

function z = df(X)
  Q = [1 -1/3; -1/3 3];
  b = [ -1 ; 2 ] ;
  z = Q * X - [ -1 ; 2 ] + [1 ; 1] * e^(X(1)+X(2));
endfunction

function z = ej1_5()
	z = fsolve(@df,[1 ; 1]);
endfunction

function v = a(X)
  Q = [1 -1/3; -1/3 3];
  b = [ -1 ; 2 ] ;
  v = [1 ; 1] * sqrt(e^(X(1)+X(2)));
endfunction

function j = jF(X)
  Q = [1 -1/3; -1/3 3];
  b = [ -1 ; 2 ] ;
  a_p = a(X);
  j = Q + a_p * a_p';
endfunction

function sol = NR(X,iters = 10)
  tmp = X;
  for i = 1:iters
    j_k = jF(tmp);
    f_k = df(tmp);
    tmp = tmp - j_k \ f_k;
  endfor
  sol = tmp;
endfunction

function error = ej2_2()
  res_1_5 = ej1_5();
  res_2_3 = NR([1;1]);
  error = abs(res_1_5 - res_2_3);
endfunction

function reses = ej2_6()
    iters = 10
  x = zeros(iters,1)
  tmp = [1;1];
  for i = 1:iters;
    tmp = NR(tmp,1);
    x(i) = norm(df(tmp));
  endfor
  reses = x
endfunction

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

function res = az(z)
  res = sqrt(exp(sum(z)))*ones(size(z))';
endfunction
function res = afz(z)
  res = exp(sum(z))*ones(size(z))';
endfunction
function res = b_n(n)
  res = ones(1,n)';
  for i = 1:n
    res(i) = (-1)^i * i;
  endfor
endfunction

function z = NR2(sQ, iters=1)
  
  z_tmp = (ones(1,sQ)/sQ);
  b = b_n(sQ)';
  Q = Q_n(sQ)
  L = cholesky_tridiagonal(Q);
  lowerlinsolve = opts.LT = true;
  upperlinsolve = opts.UT = true;
  for i = 1:iters
    a = az(z_tmp);
    aF = afz(z_tmp);
    [Lm,info] = cholupdate(L',a);
    %NR: F(z_new-z_tmp) = F(z_tmp) + dF(z_tmp) * (z_new-z_tmp)
    % igualando a 0
    % z_new = z_tmp - dF^-1(z_tmp)(F(z_tmp))
    % 
    %(Q+b)x = -F(z_tmp)
    % LL'x = -F(z_tmp)
    % Ly = -F(z_tmp)
    % L' z_newtmp = y
    Fz = Q*z_tmp' - b' + aF;
    y = linsolve(Lm,Fz);
    x = linsolve(Lm',y);
    z_tmp = z_tmp - x';
  endfor
  z = z_tmp';
endfunction

function res1,res2 = ej5_1_2()
  enes = [1,2,3,4,5,10,20,30,40,50]*100
  valores = zeros(size(enes),1)
  normas = zeros(size(enes),1)
  tmp = ones(1,size(Q)(1))/size(Q)(1)
  for i = 1:size(enes);
    %
    tmp = NR(tmp,1);
    x(i) = norm(df(tmp));
  endfor
  res1 = x
endfunction
