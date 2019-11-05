X = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0];
Y1 = [3.89,2.75,2.01,1.61,1.21,0.89,0.69,0.63,0.44,0.42,0.70,0.32,0.40,0.26,0.32,0.25];
Y2 = [15.96,9.45,5.75,3.82,2.89,2.17,1.22,1.05,0.86,0.63,0.69,0.40,0.44,0.29,0.43,0.20];
X2 = X;
X2(11) = [];
Y12 = Y1;
Y12(11)= [];

g = @(ci,p) @(x) ci/x^p;
lsval = @(X,Y) @(PC) norm(Y-PC(2)*X.^(-PC(1)));

function res = Dhat(ci,p,X,Y)
  res = 0;
  g = @(ci,p) @(x) ci/x^p;
  for x = 1:length(X)
    res += (g(ci,p)(X(x)) - Y(x))^2;
  endfor
endfunction


function res = DifDhat(ci,p,X,Y)
  res = -2*ci;
  acc = 0;
  g = @(ci,p) @(x) ci/x^p;
  for x = 1:length(X)
    acc += (g(ci,p)(X(x)) - Y(x))* log(X(x))/X(x)^p;
  endfor
  res *= acc;
endfunction

function pn = NRstep(ci,oldp,X,Y)
  pn = oldp - Dhat(ci,oldp,X,Y)/DifDhat(ci,oldp,X,Y);
endfunction

function PCk = PMCNL(PCinit,X,Y,iter=10)
  fun = @(PC) PC(2)*X.^(-PC(1))
  diffun = @(PC) [(-PC(2)*X.^(-PC(1)).*log(X))' (X.^(-PC(1)))']
  PCk = PCinit;
  for i = 1:iter 
    Ak = diffun(PCk);
    Yk = Y' - fun(PCk)';
    PCk = PCk + ols(Yk,Ak);
  endfor
endfunction

function plotbasic(ci,p,X,Y)
  g = @(ci,p) @(x) ci/x^p;
  newplot
  hold on
  scatter(X,Y)
  plot(X,arrayfun(g(ci,p),X))
  hold off
endfunction

function plotiterNR(X,Y,p0,iter=10)
  g = @(ci,p) @(x) ci/x^p;
  newplot
  hold on
  #scatter(X,Y);
  pn = p0;
  color = [1,1,1];
  for i = 1:iter
    #color = [1,1,1]*log(i)/log(iter);
    ci = ols(Y',arrayfun(g(1,pn),X'))
    pn = NRstep(ci,pn,X,Y)
    scatter(i,Dhat(ci,pn,X,Y))
    #plot(X,arrayfun(g(ci,pn),X),'color',color)
    text(i,Dhat(ci,pn,X,Y),mat2str(i))
  endfor
  hold off
endfunction



#todo mal
function errorsfc(X,Y,P,C)
  err = zeros(length(P),length(C))
  for i=P
    for j=C
      err(i,j) = Dhat(i,j,X,Y)
    endfor
  endfor
  p = linspace (0.1, 3, 41)';
  c = linspace(-10,10,41)';
  [pp, cc] = meshgrid (p,c);
  r = sqrt ( (cc/x^pp - y)^2);
  tz = sin (r) ./ r;
  mesh (tx, ty, tz);
  xlabel ("tx");
  ylabel ("ty");
  zlabel ("tz");
  title ("3-D Sombrero plot");
endfunction

function plotspline(X,Y)
  
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

  ecdif = @(x,y) -c*y/x^2 + d/x^3;
  
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
  
  
  [x_RK, y_RK] = ode45(ecdif, [a b], y0);
  plot(x_RK, y_RK, "k")
  sp = linspace(X(1),X(end));
  spp = spline(x_RK,y_RK,sp);
  hold on
  plot(sp,spp);
  plot(x_RK,y_RK);
  hold off
endfunction

function plotinterp(X,Y)
  
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

  ecdif = @(x,y) -c*y/x^2 + d/x^3;
  
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
  
  
  [x_RK, y_RK] = ode45(ecdif, [a b], y0);
  sp = linspace(X(1),X(end));
  spp = interp1(x_RK,y_RK,sp);
  hold on
  plot(sp,spp);
  plot(x_RK, y_RK, "k")
  hold off
endfunction