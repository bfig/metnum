function z = df(X)
	z = X' * [1 -1/3; -1/3 3] * X + [1 1]' * e^(X(1)+X(2))
end




