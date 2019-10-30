X = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0];
%Y1
Y1 = [3.89,2.75,2.01,1.61,1.21,0.89,0.69,0.63,0.44,0.42,0.70,0.32,0.40,0.26,0.32,0.25];

%Y2
%Y1 = [15.96,9.45,5.75,3.82,2.89,2.17,1.22,1.05,0.86,0.63,0.69,0.40,0.44,0.29,0.43,0.20];

%   n cantidad de muestras 
n = length(X);



% Implementacion de Euler hacia adelante
function sol = EulerAdelante(h=0.1,iters = 16)
  tmpX = 1/2
  tmpY = 0
  d=1.95;
  c=0.97;
  for i = 1:iters
    tmpY = tmpY + h * (d/tmpX^3 - (c/tmpX^2)* tmpY )
    tmpX = tmpX + h
  endfor
  %sol = tmpY;
endfunction
