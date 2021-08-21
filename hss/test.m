%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% test/driver routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

n = 2^10;
v = randn(n,1);
T = toeplitz(v); % symmetric random Toeplitz matrix
A = ifft((ifft(T)*sqrt(n))')*sqrt(n);

n = size(A,1);
m0 = 64; % HSS block row size
tol = 1e-6;
[tr,m] = npart(n,m0);  % binary tree; partition sequence of n

%%%%%%%%%% Generation of HSS approximation
[D,U,R,B,W,V,nflops1] = mat2hss(A,tr,m,'tol',tol);

%%%%%%%%%% HSS system solution test
x0 = randn(n,1); b = A*x0;
[x,nflops2] = hssulvsol(tr,D,U,R,B,W,V,length(tr),b);
fprintf('HSS solution      flops: %e\n',nflops2);
fprintf('Regular solution  flops: %e\n',2/3*n^3);
fprintf('relative       residual: %e\n',norm(A*x-b)/norm(b));

%%%%%%%%%% Test the accuracy of the HSS matrix approximation
A1 = hss2mat(D,U,R,B,W,V,tr);
fprintf('||A-tilde A||_2/||A||_2: %e\n',norm(A-A1)/norm(A));