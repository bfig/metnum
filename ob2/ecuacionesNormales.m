X = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0];
%Y1
Y1 = [3.89,2.75,2.01,1.61,1.21,0.89,0.69,0.63,0.44,0.42,0.70,0.32,0.40,0.26,0.32,0.25];

%Y2
%Y1 = [15.96,9.45,5.75,3.82,2.89,2.17,1.22,1.05,0.86,0.63,0.69,0.40,0.44,0.29,0.43,0.20];

%   n cantidad de muestras 
n = length(X);

log_X = zeros(1,length(X));
log_Y1 = zeros(1,length(Y1));
log_Xxlog_Y1 = zeros(1,length(X));
log_Xxlog_X = zeros(1,length(X));
Sumlog_Y1 = 0;
Sumlog_X = 0;
Sumlog_Xxlog_Y1 = 0;
Sumlog_Xxlog_X = 0;

  for i = 1:length(X)
    log_X(i) = log(X(i));
    log_Y1(i) = log(Y1(i));
    log_Xxlog_Y1(i) = log(X(i)) * log(Y1(i));
    log_Xxlog_X(i) = log(X(i)) * log(X(i));
    
    Sumlog_Y1 = Sumlog_Y1 + log_Y1(i);
    Sumlog_X = Sumlog_X + log_X(i);
    Sumlog_Xxlog_Y1 = Sumlog_Xxlog_Y1 + log_Xxlog_Y1(i);
    Sumlog_Xxlog_X = Sumlog_Xxlog_X + log_Xxlog_X(i);
  endfor

%ecuaciones normales
% Sumlog_Y1 = n * log c + p * Sumlog_X    (1)
% Sumlog_Xxlog_Y1 = log c * Sumlog_X  + p * Sumlog_Xxlog_X (2)

p=round((Sumlog_Y1 * Sumlog_X - n * Sumlog_Xxlog_Y1) /( n * Sumlog_Xxlog_X - Sumlog_X * Sumlog_X ))
c= exp((Sumlog_Y1 + p * Sumlog_X ) /n)

