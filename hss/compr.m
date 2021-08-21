function [Q,R,nflops,rk] = compr(A,type,par,opt,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compression routine (rank-revealing QR) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ATTN: ONLY opt 0 supplied in this package for simplicity!
% opt = 0: Default: truncated SVD; may be slow in large matlab test
% opt = 1: A C mex rank-revealing QR with column pivoting
% opt = 2: A matlab rank-revealing QR with column pivoting; may be slow in large matlab test
% opt = 3: Randomized SVD

if nargin < 3
    error('3 or more inputs are needed')
elseif ~strcmp(type,'tol') && ~strcmp(type,'rank')
    error('The 2nd input must be ''tol'' or ''rank''')
end

if nargin < 4
    opt = 0;
end

if nnz(A) == 0 % full zero matrix; need to address this in mgsclpv.cpp
    Q = zeros(size(A,1),1);
    R = zeros(1,size(A,2));
    nflops = 0;
    rk = 1;
    return
end

%opt = 4;
switch opt
    case 0
        nflops = NaN; % No flop count
        [Q,R,V] = svd(A);
        if strcmp(type,'tol')
            rk = length(find(diag(R) > par*R(1,1)));
        else
            rk = par;
        end
        Q = Q(:,1:rk); R = R(1:rk,1:rk); V = V(:,1:rk);
        R = R*V';
    case 1
        if strcmp(type,'tol')
            [Q,R,nflops,rk] = mgsclpv(A,par);
        else
            [Q,R,nflops,rk] = mgsclpvr(A,par);
        end
    case 2
        [Q,R,nflops,rk] = mgsclpvm(A,typ,par);
        rk = min([rk,size(R)]);
        Q = Q(:,1:rk); R = R(1:rk,:);
    case 3
        if strcmp(type,'rank')
            [Q,R,nflops] = rsvd(A,r,p,k,varargin);
        else
            error('Type ''tol'' for rsvd method unavailable')
        end
end
% fprintf('compression flops: %8.2e\n',nflops);
% maxerr(Q*R,A)