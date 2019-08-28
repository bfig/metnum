function z = df(X)
  z = X' * Q() * X + [1 1]' * e^(X(1)+X(2))
end




