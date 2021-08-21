function [D,U,R,B,nflops] = mat2hsssym(A,tr,m,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Symmetric HSS construction %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input parameters
% A   --- symmetric matrix
% m   --- partition information
% tr  --- HSS tree information, tr(i) is the parent of i
% tol --- tolerance
%
%%% Output parameters
% D,U,R,B --- generator arrays (cell arrays)
%
%%% Example
% A = rand(128,128); A = A+A';
% m = [32 32 32 32];
% tr = [3 3 7 6 6 7 0];
% or [tr,m] = npart(128,32);
% tol = 1e-6;
% [D,U,R,B] = mat2hsssym(A,tr,m,'tol',1e-6);

nflops = 0;
N = size(A,1);

if nargin == 1
    error('Input tr,m needed')
elseif nargin == 3
    ctype = 'tol';
    par = 1e-6;
elseif length(varargin) == 2
    ctype = varargin{1};
    par = varargin{2};
else
    error('Not enough input')
end

D = {[]}; U = {[]}; R = {[]}; B = {[]};

if length(m) == 1
    D{1} = A;
    return
end

n = length(tr);
ch = child(tr);

l = zeros(n,2);
l(1,:) = [1 m(1)];
lt = 1; it = 1;
for i = 1:n
    if isempty(ch{i})
        l(i,:) = [lt lt+m(it)-1];
        lt = l(i,2)+1; it = it+1;
    else
        l(i,:) = [l(ch{i}(1),1) l(ch{i}(2),2)];
    end
end

st = []; % stack
ns = 0;

%%% D, U
for i = 1:n-1
    %fprintf('-------- %d -------\n',i);
    if isempty(ch{i})
        D{i} = A(l(i,1):l(i,2),l(i,1):l(i,2));

        T{i} = [];
        if i > 1
            st = [st i-1]; ns = ns+1;
            T{i} = [];
            for j = 1:ns
                T{i} = [T{i} T{st(j)}(:,end-(N-l(i,1)):end-(N-l(i,2)))'];
            end
        end
        
        T{i} = [T{i} A(l(i,1):l(i,2),l(i,2)+1:N)];
        [U{i},T{i},nflops1] = compr(T{i},ctype,par);
        nflops = nflops+nflops1;
    else
        st(end) = []; ns = ns-1;
        
        i1 = ch{i}(1); i2 = ch{i}(2);
        [r1,n1] = size(T{i1}); [r2,n2] = size(T{i2});

        T{i} = [T{i1}(:,1:n1-(N-l(i1,2))) T{i1}(:,n1-(N-l(i,2))+1:end);
                T{i2}(:,1:n1-(N-l(i1,2))) T{i2}(:,n2-(N-l(i,2))+1:end)];
        B{i1} = T{i2}(:,n1-(N-l(i1,2))+1:n2-(N-l(i,2)))';
        T{i1} = []; T{i2} = [];
        [U{i},T{i},nflops1] = compr(T{i},ctype,par);
        nflops = nflops+nflops1;
    end
end

i = n;
% fprintf('-------- %d -------\n',i);
i1 = ch{i}(1); i2 = ch{i}(2);
[r1,n1] = size(T{i1}); [r2,n2] = size(T{i2});
B{i1} = T{i2}(:,n1-(N-l(i1,2))+1:n2-(N-l(i,2)))';

%%% R
for i = n-1:-1:1
    if ~isempty(ch{i})
        i1 = ch{i}(1); i2 = ch{i}(2);
        sz = size(U{i1},2);
        R{i1} = U{i}(1:sz,:);
        R{i2} = U{i}(sz+1:end,:);
        U{i} = [];
    end
end

% fprintf('mat2hsssym nflops: %8.2e\n',nflops);
return;