function [A,nflops] = hss2matsym(D,U,R,B,tr)

nflops = 0;
n = length(tr);

if n == 1
    A = D{1};
    return
end

lv = hsslevel(tr);

ch = child(tr);

cl = [];
for i = 1:n
    if isempty(ch{i})
        cl = [cl i];
    end
end

n = length(cl);
for i = 1:n
    m(i) = size(U{cl(i)},1);
end
l(1) = 1;
for j = 1:length(m)
    l(j+1) = l(j)+m(j);
end

A = zeros(l(end)-1,l(end)-1);
for i = 1:n-1
    ii = cl(i);
    A(l(i):l(i+1)-1,l(i):l(i+1)-1) = D{ii};
    for j = i+1:n
        ii = cl(i);    jj = cl(j);
        % fprintf('[i,j]: [%d,%d]; [ii,jj]: [%d,%d]\n',i,j,ii,jj);
            T1 = U{ii}; T2 = U{jj};
            while lv(ii) > lv(jj)
                nflops = nflops+flops('prod',T1,'n',R{ii},'n');
                T1 = T1*R{ii}; ii = tr(ii);
            end
            while lv(ii) < lv(jj)
                nflops = nflops+flops('prod',T2,'n',R{jj},'n');
                T2 = T2*R{jj}; jj = tr(jj);
            end
            while tr(ii) ~= tr(jj)
                nflops = nflops+flops('prod',T1,'n',R{ii},'n')+flops('prod',T2,'n',R{jj},'n');
                T1 = T1*R{ii}; ii = tr(ii);
                T2 = T2*R{jj}; jj = tr(jj);
            end
            nflops = nflops+flops('prod',T1,'n',B{ii},'n',T2,'t');
            A(l(i):l(i+1)-1,l(j):l(j+1)-1) = T1*B{ii}*T2';
            A(l(j):l(j+1)-1,l(i):l(i+1)-1) = A(l(i):l(i+1)-1,l(j):l(j+1)-1)';
    end
end
i = n;
ii = cl(i);
A(l(i):l(i+1)-1,l(i):l(i+1)-1) = D{ii};