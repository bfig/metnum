chol(Q_n(2))

% Funcion
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

