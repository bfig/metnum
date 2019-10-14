function z = f(X)
	z = X' * [1 -1/3; -1/3 3] * X - 2* [ -1 2 ] * X + 2* e^(X(1)+X(2))
end