function [P,L,Q,DD,UU,VV,sz,nflops] = hssulv(tr,D,U,R,B,W,V,npv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HSS ULV factorization %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(tr);
ch = child(tr);

if nargin < 8
    npv = length(tr);
end

nflops = 0;

S = [];
l = 1;
for i = 1:npv
%     fprintf('-------- %d -------\n',i);
    if isempty(ch{i})
        Di = D{i}; Ui = U{i}; Vi = V{i};
        l = l+size(D{i},1);
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
        L{i} = Di;
        break;
    end
    
    sz(i,1:2) = size(Ui);
    [P{i},Ui] = qr(Ui); P{i} = P{i}(:,[sz(i,2)+1:sz(i,1) 1:sz(i,2)]); Ui = Ui(1:sz(i,2),:);
    UU{i} = Ui;
    nflops = nflops+2*sz(i,1)*sz(i,2)^2-2/3*sz(i,2)^3;
    Di = P{i}'*Di;
    nflops = nflops+flops('prod',P{i},'t',Di,'n');
    
    % LQ factorization of partial diag
    [Q{i},L{i}] = qr(Di(1:sz(i,1)-sz(i,2),:)'); L{i} = L{i}(1:sz(i,1)-sz(i,2),:)';
    nflops = nflops+2*size(Di,1)*(sz(i,1)-sz(i,2))^2-2/3*(sz(i,1)-sz(i,2))^3;
    Di(sz(i,1)-sz(i,2)+1:end,:) = Di(sz(i,1)-sz(i,2)+1:end,:)*Q{i};
    nflops = nflops+2*sz(i,1)*sz(i,2)*(2*sz(i,1)-sz(i,2));
    
    DD{i} = Di;
   
    Vi = Q{i}'*Vi;
    VV{i} = Vi;
    nflops = nflops+2*size(Vi,2)*sz(i,2)*(2*sz(i,1)-sz(i,2));

    Vi = Vi(sz(i,1)-sz(i,2)+1:end,:);
    Di = Di(sz(i,1)-sz(i,2)+1:end,sz(i,1)-sz(i,2)+1:end);
    if i == ch{tr(i)}(1)
        S = push(S,Vi); S = push(S,Ui); S = push(S,Di);
    end
end

% Schur complement
if npv == n-2
    DD{n-1} = D{n-1}-U{n-1}*B{n-1}*Vi'/Di*Ui*B{n-2}*V{n-1}';
end
% fprintf('hssulv nflops: %e\n',nflops);