function X = project_psd(X)

[V, D] = eig(X);
X = V * diag(max(0, diag(D))) * V';

end
