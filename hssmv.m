function [y,nflops] = hssmv(D,U,R,B,W,V,tr,x,trans)
% trans (optional): with any input trans, it is to compute A^T*x

k = length(tr);

if k == 1
    y = D{1}*x;
    nflops = flops('mv',D{1},'n');
    return
end

ch = child(tr);

nflops = 0;

if nargin == 8
    l = 1;
    for i = 1:k-1
        if isempty(ch{i})
            g{i} = V{i}'*x(l:l+size(D{i},1)-1,:);
            l = l+size(D{i},1);
            nflops = nflops+flops('mv',V{i},'t');
        else
            c1 = ch{i}(1); c2 = ch{i}(2);
            g{i} = W{c1}'*g{c1}+W{c2}'*g{c2};
            nflops = nflops+flops('mv',W{c1},'t')+flops('mv',W{c2},'t')+length(g{i});
        end
    end
    
    y = zeros(l-1,1);
    for i = k-1:-1:1
        s = sib(tr,ch,i);
        if tr(i) == k
            f{i} = B{i}*g{s};
            nflops = nflops+flops('mv',B{i},'n');
        else
            f{i} = B{i}*g{s}+R{i}*f{tr(i)};
            nflops = nflops+flops('mv',B{i},'n')+flops('mv',R{i},'n')+length(f{i});
        end
        
        if isempty(ch{i})
            y(l-size(D{i},1):l-1,1:size(x,2)) = D{i}*x(l-size(D{i},1):l-1,:)+U{i}*f{i};
            l = l-size(D{i},1);
            nflops = nflops+flops('mv',D{i},'n')+flops('mv',U{i},'n')+size(D{i},1);
        end
    end    
else
    l = 1;
    for i = 1:k-1
        if isempty(ch{i})
            g{i} = U{i}'*x(l:l+size(D{i},1)-1,:);
            l = l+size(D{i},1);
            nflops = nflops+flops('mv',U{i},'t');
        else
            c1 = ch{i}(1); c2 = ch{i}(2);
            g{i} = R{c1}'*g{c1}+R{c2}'*g{c2};
            nflops = nflops+flops('mv',R{c1},'t')+flops('mv',R{c2},'t')+length(g{i});
        end
    end
    
    y = zeros(l-1,1);
    for i = k-1:-1:1
        s = sib(tr,ch,i);
        if tr(i) == k
            f{i} = B{s}'*g{s};
            nflops = nflops+flops('mv',B{s},'t');
        else
            f{i} = B{s}'*g{s}+W{i}*f{tr(i)};
            nflops = nflops+flops('mv',B{s},'t')+flops('mv',W{i},'n')+length(f{i});
        end
        
        if isempty(ch{i})
            y(l-size(D{i},1):l-1,1:size(x,2)) = D{i}'*x(l-size(D{i},1):l-1,:)+V{i}*f{i};
            l = l-size(D{i},1);
            nflops = nflops+flops('mv',D{i},'t')+flops('mv',V{i},'n')+size(D{i},1);
        end
    end
end