function [r,rr] = hssrank(A,tol,lmax)
% rr: rank pattern

n = size(A,1);

if nargin < 3
    lmax = 4;
end
if nargin < 2
    tol = 1e-6;
end

m1 = floor(n/2);
m = [m1 n-m1];
r0 = min(m);

fprintf('size A: %i x %i\n',size(A));
l = 1;
rr = [];
while min(m) >= r0/2 && l <= lmax
    r1 = zeros(length(m),1);
    m1 = 0;
    for i = 1:length(m)
        T = [A(m1+1:m1+m(i),1:m1) A(m1+1:m1+m(i),m1+m(i)+1:end)];
        r1(i) = rankrel(T,tol);
        if nargout ~= 1
            fprintf('%i ',r1(i));
        end
        m1 = m1+m(i);
    end
    r0 = min(r0,max(r1));
    rr = [rr max(r1)];
    if nargout ~= 1
        fprintf('\n---- level %i, block size %i, rank %i ----\n',l,min(m),max(r1));
    end
    m2 = [];
    for i = 1:length(m)
        m1 = floor(m(i)/2);
        m2 = [m2 m1 m(i)-m1];
    end
    m = m2;
    l = l+1;
end

r = max(rr);
fprintf('HSS rank: %i\n',r);
