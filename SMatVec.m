function A = SMatVec(v)
% Determines the n*n symmetric matrix from an n*(n+1)/2 column vector
% containing the independent elements of the symmetric matrix.

    % Determine the size of the matrix
    L=length(v);
    n=round((sqrt(1+8*L)-1)/2);         % Solve n(n+1)/2 = L
    
    % Initialize the symmetric matrix
    A = zeros(n, n);
    
    % Index to track the position in the vector
    idx = 1;
    
    % Fill the diagonal elements
    for i = 1:n
        A(i, i) = v(idx);
        idx = idx + 1;
    end
    
    % Fill the upper triangular elements
    for i = 1:n-1
        for j = i+1:n
            A(i, j) = v(idx);
            A(j, i) = v(idx); % Ensure symmetry
            idx = idx + 1;
        end
    end
end