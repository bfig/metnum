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
function error = ej2_5()
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

function [ci,si,L] = compute(L,V,i)
  w = sqrt(L(i,i)^2 + V(i)^2)
  ci = w/L(i,i)
  si = V(i)/L(i,i)
  L(i,i) = w
end

function [L,V] = apply(ci,si,L,V,i,j)
  L(i,j) = (L(i,j)+si*V(j))/ci
  V(j) = ci * V(j) - si * L(i,j)
endfunction


function z = NR2(Q, iters = 10)
    tmp = X;
  for i = 1:iters
    chol
    j_k = jF(tmp);
    f_k = df(tmp);
    tmp = tmp - j_k \ f_k;
  endfor
  z = tmp;
endfunction

