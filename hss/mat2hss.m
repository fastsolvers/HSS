function [D,U,R,B,W,V,nflops] = mat2hss(A,tr,m,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  General HSS construction  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input:
% A:           dense matrix
% tr:          binary tree
% m:           partition information; example: [tr,m] = npart(1024,64);
% varargin:    compression info; currently default version in the test
%%% Output:
% D,U,R,B,W,V: HSS generators

nflops = 0;

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

D = {[]}; U = {[]}; V = {[]}; R = {[]}; W = {[]}; B = {[]};

n = length(tr);

if n == 1
    D{1} = A;
    return
end

N = size(A,1);
ch = child(tr);

% l(i,1:2): range of block i
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

ns1 = 1; % (# of blocks in stack)+1
ls1(ns1) = 1;
ws1(ns1) = l(1,2)+1;

ns2 = 1; % (# of blocks in stack)+1
ls2(ns2) = 1;
ws2(ns2) = l(1,2)+1;

for i = 1:n
    %fprintf('-------- %d -------\n',i);
    ch1 = ch{i};
    
    if isempty(ch1)
        D{i} = A(l(i,1):l(i,2),l(i,1):l(i,2));
        
        % off-diag row compression
        if ns1 == 1
            T1 = A(l(i,1):l(i,2),l(i,2)+1:N);
        else
            T1 = [A(ws2(ns2-1):ws2(ns2-1)+(l(i,2)-l(i,1)),1:ls2(ns2)-1) A(l(i,1):l(i,2),l(i,2)+1:N)];
        end
        r1 = size(T1,1);
        % compress
        [U{i},T1,nflops1] = compr(T1,ctype,par);
        nflops = nflops+nflops1;
        r1 = r1-size(T1,1);
        
        ls1(ns1+1) = ls1(ns1)+size(T1,1); % expand ls
        ws1(ns1) = l(i,2)+1;
        A(ls1(ns1):ls1(ns1+1)-1,l(i,2)+1:N) = T1(:,end-(N-l(i,2))+1:end);
  
        % off-diag col compression
        if ns2 == 1
            T2 = A(l(i,2)+1:N,l(i,1):l(i,2));
        else
            T2 = [A(1:ls1(ns1)-1,ws1(ns1-1):ws1(ns1-1)+(l(i,2)-l(i,1))); A(l(i,2)+1:N,l(i,1):l(i,2))];
        end
        r2 = size(T2,2);
        % compress
        [V{i},T2,nflops1] = compr(T2',ctype,par); T2 = T2';
        nflops = nflops+nflops1;
        r2 = r2-size(T2,2);
        
        ls2(ns2+1) = ls2(ns2)+size(T2,2); % expand ls
        ws2(ns2) = l(i,2)+1;
        A(l(i,2)+1:N,ls2(ns2):ls2(ns2+1)-1) = T2(end-(N-l(i,2))+1:end,:);

        ns1 = ns1+1;
        ns2 = ns2+1;
    else
        B{ch1(1)} = A(ls1(ns1-2):ls1(ns1-1)-1,ws1(ns1-2):ws1(ns1-1)-1);
        B{ch1(2)} = A(ws2(ns2-2):ws2(ns2-1)-1,ls2(ns2-2):ls2(ns2-1)-1);
        if i == n; break; end
        
        % off-diag row compression
        if ns1 == 3
            T1 = A(ls1(ns1-2):ls1(ns1)-1,l(i,2)+1:N);
        else
            T1 = [A(ws2(ns2-3):ws2(ns2-1)-1,1:ls2(ns2-2)-1) A(ls1(ns1-2):ls1(ns1)-1,l(i,2)+1:N)];
        end
        r1 = size(T1,1);
        % compress
        [U{i},T1,nflops1] = compr(T1,ctype,par);
        nflops = nflops+nflops1;
        r1 = r1-size(T1,1);
        ls1(ns1-1) = ls1(ns1-2)+size(T1,1);
        A(ls1(ns1-2):ls1(ns1-1)-1,l(i,2)+1:N) = T1(:,end-(N-l(i,2))+1:end);

        % off-diag col compression
        if ns2 == 3
            T2 = A(l(i,2)+1:N,ls2(ns2-2):ls2(ns2)-1);
        else
            T2 = [A(1:ls1(ns1-2)-1,ws1(ns1-3):ws1(ns1-1)-1); A(l(i,2)+1:N,ls2(ns2-2):ls2(ns2)-1)];
        end
        r2 = size(T2,2);
        % compress
        [V{i},T2,nflops1] = compr(T2',ctype,par); T2 = T2';
        nflops = nflops+nflops1;
        r2 = r2-size(T2,2);
        ls2(ns2-1) = ls2(ns2-2)+size(T2,2);
        A(l(i,2)+1:N,ls2(ns2-2):ls2(ns2-1)-1) = T2(end-(N-l(i,2))+1:end,:);

        ns1 = ns1-1; % shrink ls
        ns2 = ns2-1; % shrink ls
    end
   
    for j = 1:ns2-3
        A(ws2(j)+r1:ws2(ns1-2)-1+r1,ls2(j):ls2(j+1)-1) = A(ws2(j):ws2(ns2-2)-1,ls2(j):ls2(j+1)-1);
        ws2(j) = ws2(j)+r1;
    end
    ws2(ns2-1) = l(i,2)+1;
    if ns2 > 2
        A(ws2(ns2-2)+r1:ws2(ns2-1)-1,1:ls2(ns2-1)-1) = T1(:,1:ls2(ns2-1)-1);
        ws2(ns2-2) = ws2(ns2-2)+r1;
    end
    
    for j = 1:ns1-3
        A(ls1(j):ls1(j+1)-1,ws1(j)+r2:ws1(ns1-2)-1+r2) = A(ls1(j):ls1(j+1)-1,ws1(j):ws1(ns1-2)-1);
        ws1(j) = ws1(j)+r2;
    end
    ws1(ns1-1) = l(i,2)+1;
    if ns1 > 2
        A(1:ls1(ns1-1)-1,ws1(ns1-2)+r2:ws1(ns1-1)-1) = T2(1:ls1(ns1-1)-1,:);
        ws1(ns1-2) = ws1(ns1-2)+r2;
    end
end

for i = n-1:-1:1
    ch1 = ch{i};
    if ~isempty(ch1)
        sz = size(U{ch1(1)},2);
        R{ch1(1)} = U{i}(1:sz,:);
        R{ch1(2)} = U{i}(sz+1:end,:);
        U{i} = [];
        
        sz = size(V{ch1(1)},2);
        W{ch1(1)} = V{i}(1:sz,:);
        W{ch1(2)} = V{i}(sz+1:end,:);
        V{i} = [];
    end
end