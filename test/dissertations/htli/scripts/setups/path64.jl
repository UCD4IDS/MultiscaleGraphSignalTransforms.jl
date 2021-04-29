N = 64; G = path_graph(N)
L = Matrix(laplacian_matrix(G))
𝛌, 𝚽 = eigen(L); 𝚽 = 𝚽 .* sign.(𝚽[1,:])'
Q = incidence_matrix(G; oriented = true)
∇𝚽 = Q' * 𝚽
