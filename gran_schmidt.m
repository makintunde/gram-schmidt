function[result] = gran_schmidt(varargin)
    A_raw = cellfun(@transpose, varargin, 'UniformOutput', false);  
    A = cell2mat(A_raw);
    V(:,1) = A(:, 1);
    Q(:,1) = V(:,1) / norm(V(:,1));
    %args = zeros(size(A, 1),1);
    for i=2:nargin,
        V(:,i) = A(:, i) - sum(A, Q, i, nargin);
        Q(:,i) = V(:, i) / norm(V(:,i));
    end
    result = (Q' * Q) - eye(nargin);
end

function[result] = sum(A, Q, i, n)
    R = zeros(size(A, 1), 1);
    for j=1:i-1,
        R = R + dot(A(:,i), Q(:,j)) * Q(:,j);
    end    
    result = R;
end
