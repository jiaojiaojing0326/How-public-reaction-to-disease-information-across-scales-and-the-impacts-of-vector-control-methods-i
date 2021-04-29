function D=distance(A)
% calculate de distance matrix from a given adjancency matrix A

N = max(size(A));
D = NaN(N);
B = A;
k = 1;
%while any(isnan(D(:)))
while k <= N
    % Check for new walks, and assign distance
    D(B > 0 & isnan(D)) = k;
    % Iteration
    k = k + 1;
    B = B * A;
end
D(1:N+1:end) = 0; % diagonal elements are always 0
% Now D contains the distance matrix
end 