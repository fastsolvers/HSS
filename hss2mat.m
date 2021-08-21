function [A,nflops] = hss2mat(D,U,R,B,W,V,tr)

n = length(tr);
lv = hsslevel(tr);
ch = child(tr);
nflops = 0;

cl = [];
k = 0;
l(1) = 1;
for i = 1:n
    if isempty(ch{i})
        k = k+1;
        cl(k) = i;
        m(k) = size(D{i},1);
        l(k+1) = l(k)+m(k);
    end
end

A = zeros(l(end)-1,l(end)-1);
for i = 1:k
    ii = cl(i);
    A(l(i):l(i+1)-1,l(i):l(i+1)-1) = D{ii};
    
    for j = [1:i-1 i+1:k]
        ii = cl(i); jj = cl(j);
        T1 = U{ii}; T2 = V{jj};
        while lv(ii) > lv(jj)
            sz1 = size(T1); sz2 = size(R{ii});
            nflops = nflops + sz1(1)*(2*sz1(2)-1)*sz2(2);
            T1 = T1*R{ii}; ii = tr(ii);
        end
        while lv(ii) < lv(jj)
            sz1 = size(T2); sz2 = size(W{jj});
            nflops = nflops + sz1(1)*(2*sz1(2)-1)*sz2(2);
            T2 = T2*W{jj}; jj = tr(jj);
        end
        while tr(ii) ~= tr(jj)
            sz1 = size(T1); sz2 = size(R{ii});
            nflops = nflops + sz1(1)*(2*sz1(2)-1)*sz2(2);
            T1 = T1*R{ii}; ii = tr(ii);
            sz1 = size(T2); sz2 = size(W{jj});
            nflops = nflops + sz1(1)*(2*sz1(2)-1)*sz2(2);
            T2 = T2*W{jj}; jj = tr(jj);
        end
        sz1 = size(T1); sz2 = size(B{ii}); sz3 = size(T2');
        nflops = nflops + sz1(1)*(2*sz1(2)-1)*sz2(2) + sz1(1)*(2*sz2(2)-1)*sz3(2);
        A(l(i):l(i+1)-1,l(j):l(j+1)-1) = T1*B{ii}*T2';
    end
end