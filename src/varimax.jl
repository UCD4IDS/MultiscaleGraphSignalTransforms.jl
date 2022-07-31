"""
	varimax(A; gamma = 1.0, minit = 20, maxit = 1000, reltol = 1e-12)

VARIMAX perform varimax (or quartimax, equamax, parsimax) rotation to the column vectors of the input matrix.

# Input Arguments
- `A::Matrix{Float64}`: input matrix, whose column vectors are to be rotated. d, m = size(A).
- `gamma::Float64`: default is 1. gamma = 0, 1, m/2, and d(m - 1)/(d + m - 2), corresponding to quartimax, varimax, equamax, and parsimax.
- `minit::Int`: default is 20. Minimum number of iterations, in case of the stopping criteria fails initially.
- `maxit::Int`: default is 1000. Maximum number of iterations.
- `reltol::Float64`: default is 1e-12. Relative tolerance for stopping criteria.

# Output Argument
- `B::Matrix{Float64}`: output matrix, whose columns are already been rotated.

Implemented by Haotian Li, Aug. 20, 2019
"""
function varimax(A; gamma = 1.0, minit = 20, maxit = 1000, reltol = 1e-12)
	# Get the sizes of input matrix
	d, m = size(A)

	# If there is only one vector, then do nothing.
	if m == 1
		return A
	end

	if d == m && rank(A) == d
		return Matrix{Float64}(I, d, m)
	end

	# Warm up step: start with a good initial orthogonal matrix T by SVD and QR
	T = Matrix{Float64}(I, m, m)
	B = A * T
	L,_,M = svd(A' * (d*B.^3 - gamma*B * Diagonal(sum(B.^2, dims = 1)[:])))
	T = L * M'
	if norm(T-Matrix{Float64}(I, m, m)) < reltol
		T = qr(randn(m,m)).Q
		B = A * T
	end

	# Iteration step: get better T to maximize the objective (as described in Factor Analysis book)
	D = 0
	for k in 1:maxit
		Dold = D
		L,s,M = svd(A' * (d*B.^3 - gamma*B * Diagonal(sum(B.^2, dims = 1)[:])))
		T = L * M'
		D = sum(s)
		B = A * T
		if (abs(D - Dold)/D < reltol) && k >= minit
			break
		end
	end

	# Adjust the sign of each rotated vector such that the maximum absolute value is positive.
	for i in 1:m
		if abs(maximum(B[:,i])) < abs(minimum(B[:,i]))
			B[:,i] .= - B[:,i]
		end
	end

	return B
end
