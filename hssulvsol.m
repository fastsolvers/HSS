function [x,nflops] = hssulvsol(tr,D,U,R,B,W,V,npv,b)
%%%%%%%%%%%%%%%%%%%%%%
%%%  HSS solution  %%%
%%%%%%%%%%%%%%%%%%%%%%

n = length(tr);
ch = child(tr);

nflops = 0;

xt(1:n) = {[]};
S = [];
l = 1;
for i = 1:npv
%     fprintf('-------- %d -------\n',i);
    if isempty(ch{i})
        Di = D{i}; Ui = U{i}; Vi = V{i};
        xt{i} = b(l:l+size(D{i},1)-1);  l = l+size(D{i},1);
    else
        c1 = ch{i}(1); c2 = ch{i}(2);
        [S,D1] = pop(S); [S,U1] = pop(S); [S,V1] = pop(S);
        % the previous node is a right child, no need to push/pop
        T1 = U1*B{c1}*Vi'; T2 = Ui*B{c2}*V1';
        nflops = nflops+flops('prod',U1,'n',B{c1},'n')+2*size(U1,1)*numel(Vi);
        nflops = nflops+flops('prod',Ui,'n',B{c2},'n')+2*size(Ui,1)*numel(V1);
        Di = [D1 T1; T2 Di];
        if i < n
            Ui = [U1*R{c1}; Ui*R{c2}];
            Vi = [V1*W{c1}; Vi*W{c2}];
            nflops = nflops+flops('prod',U1,'n',R{c1},'n')+flops('prod',Ui,'n',R{c2},'n');
            nflops = nflops+flops('prod',V1,'n',W{c1},'n')+flops('prod',Vi,'n',W{c2},'n');
        end
    end

    if i == n
        xt{i} = Di\xt{i};
        nflops = nflops+flops('chol',Di)+2*flops('rdiv',xt{i},'t',Di,'tri');
        break;
    end
    
    sz(i,1:2) = size(Ui);
    [Qi,Ui] = qr(Ui); Qi = Qi(:,[sz(i,2)+1:sz(i,1) 1:sz(i,2)]); Ui = Ui(1:sz(i,2),:);
    nflops = nflops+2*sz(i,1)*sz(i,2)^2-2/3*sz(i,2)^3;
    Di = Qi'*Di;
    xt{i} = Qi'*xt{i};
    nflops = nflops+2*(size(Di,1)+1)*sz(i,2)*(2*sz(i,1)-sz(i,2));
    [Q{i},Li] = qr(Di(1:sz(i,1)-sz(i,2),:)'); Li = Li(1:sz(i,1)-sz(i,2),:)';
    nflops = nflops+2*size(Di,1)*(sz(i,1)-sz(i,2))^2-2/3*(sz(i,1)-sz(i,2))^3;
    Di(sz(i,1)-sz(i,2)+1:end,:) = Di(sz(i,1)-sz(i,2)+1:end,:)*Q{i};
    nflops = nflops+2*sz(i,1)*sz(i,2)*(2*sz(i,1)-sz(i,2));
   
    Vi = Q{i}'*Vi;
    nflops = nflops+2*size(Vi,2)*sz(i,2)*(2*sz(i,1)-sz(i,2));
    
    xt{i}(1:sz(i,1)-sz(i,2)) = Li\xt{i}(1:sz(i,1)-sz(i,2));
    nflops = nflops+flops('chol',Li);
    
    if isempty(ch{i})
        Vx{i} = Vi(1:sz(i,1)-sz(i,2),:)'*xt{i}(1:sz(i,1)-sz(i,2));
        nflops = nflops+2*numel(Vi(1:sz(i,1)-sz(i,2),:));
    else
        Vx{i} = Vx{i}+Vi(1:sz(i,1)-sz(i,2),:)'*xt{i}(1:sz(i,1)-sz(i,2));
        nflops = nflops+2*numel(Vi(1:sz(i,1)-sz(i,2),:))+length(Vx{i});
    end
    if i == ch{tr(i)}(1)
        xt{tr(i)} = xt{i}(sz(i,1)-sz(i,2)+1:end)-Di(sz(i,1)-sz(i,2)+1:end,1:sz(i,1)-sz(i,2))*xt{i}(1:sz(i,1)-sz(i,2),:);
        nflops = nflops+2*numel(Di(sz(i,1)-sz(i,2)+1:end,1:sz(i,1)-sz(i,2)))+length(xt{tr(i)});
    else
        Stmp = pop(S); [Stmp,Utmp] = pop(Stmp);
        xt{tr(i)} = xt{tr(i)}-Utmp*(B{ch{tr(i)}(1)}*Vx{i});
        nflops = nflops+2*numel(B{ch{tr(i)}(1)})+2*numel(Utmp)+length(xt{tr(i)});
        xt{tr(i)} = [xt{tr(i)}; 
            xt{i}(sz(i,1)-sz(i,2)+1:end)-Di(sz(i,1)-sz(i,2)+1:end,1:sz(i,1)-sz(i,2))*xt{i}(1:sz(i,1)-sz(i,2),:)-Ui*B{ch{tr(i)}(2)}*Vx{ch{tr(i)}(1)}];
        nflops = nflops+2*sz(i,2)*sz(i,1)+2*numel(Ui)+2*sz(i,2); 
        if tr(i) < n
            Vx{tr(i)} = W{ch{tr(i)}(1)}'*Vx{ch{tr(i)}(1)}+W{i}'*Vx{i};
            nflops = nflops+2*numel(W{ch{tr(i)}(1)})+2*numel(W{i});
        end
    end
    xt{i} = xt{i}(1:sz(i,1)-sz(i,2));
    Vi = Vi(sz(i,1)-sz(i,2)+1:end,:);
    Di = Di(sz(i,1)-sz(i,2)+1:end,sz(i,1)-sz(i,2)+1:end);
    if i == ch{tr(i)}(1)
        S = push(S,Vi); S = push(S,Ui); S = push(S,Di);
    end
end

for i = n-1:-1:1
    if i == ch{tr(i)}(1)
        xt{i} = Q{i}*[xt{i}; xt{tr(i)}(1:sz(i,2))];
    else
        xt{i} = Q{i}*[xt{i}; xt{tr(i)}(end-sz(i,2)+1:end)];
    end
    nflops = nflops+4*sz(i,1)*sz(i,2);
end

x = [];
for i = 1:n
    if isempty(ch{i})
        x = [x; xt{i}];
    end
end

% fprintf('hsssol nflops: %e\n',nflops);