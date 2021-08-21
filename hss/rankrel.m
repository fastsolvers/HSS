function r = rankrel(A,tol)

if nargin == 1
    tol = 1e-6;
end

s = svd(full(A));

if 0%s(1) < tol
   r = 1;
   warning(['s(1) = ' num2str(s(1)) ' < tol = ' num2str(tol)])
else
   r = length(find(s/s(1) >= tol));
end

% using s(1)*tol may yield wrong results for zero matrices