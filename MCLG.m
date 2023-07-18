function  H = MCLG(dataname,X,Y,Gv,S,Fv,beta,mu,alpha,lambda,gamma,theta) %
         % Output: Z^T

V = size(X,2);      % view number
c = max(Y);         % class number
n = length(Y);      % object number

%% =========================Initialization==========================
% Mv Wv Qv ...
for iv=1:V
    dv(iv) = size(X{iv},1); % dimension of each view
    Wv{iv} = ones(dv(iv), c);
    Mv{iv} = ones(dv(iv), c);
    Qv{iv} = ones(c, c);
    
    pv(iv) = 1.0/V;
    qv(iv) = 1.0/V;
    mv(iv) = 1.0/V;
    nv(iv) = 1.0/V;
    rv(iv) = 1.0/V;
    
    D = diag(sum(Gv{iv}));
    LGv{iv} = D-Gv{iv}; 
end

% U A
A = ones(c,c);
U = ones(n,n);

% F Z
F = ones(n, c);
%F = orth(F);

Z = ones(c, n);
%Z = orth(F)';

%% ==========================Optimization===========================
iter = 1;
maxIter = 5;
obj = zeros(1,maxIter);
objfile = strcat('result/',dataname,'Obj.csv');% record obj value
fp = fopen(objfile, 'a', 'n', 'utf8');
while iter<=maxIter
    
    % Updating Wv
    for iv = 1:V
        [UU,~,VV] = svds(X{iv}*Z'*A');
        Wv{iv} = UU*VV';
    end
    
    % Updating A
    sumWXZ = zeros(c,c);
    for iv = 1:V
        sumWXZ = sumWXZ + qv(iv)*Wv{iv}'*X{iv}*Z';
    end
    [UU,~,VV] = svds(sumWXZ);
     A = UU*VV';
       
    % Updating U
    for i=1:n
        for j=1:n
            Q(i,j)=0.5*norm(F(i,:)-F(j,:))^2;
        end
    end
    sumGv = zeros(n,n);
    for iv = 1:V
        sumGv = sumGv + mv(iv)* Gv{iv};
    end
    U = (-0.5*theta*Q+alpha*S+sumGv)/alpha;
    U = (abs(U)+abs(U'))/2.0;  
    
    % Updating Mv
    for iv = 1:V
        I = eye(dv(iv));
        Mv{iv} = 1.0/pv(iv) * inv(1.0/pv(iv)*X{iv}*X{iv}'+beta*I)*X{iv}*F;
        
        %Mv{iv}(Mv{iv}<0) = 0;
    end
    
    % Updating F
    %%%L 
    D = diag(sum(U));
    Lu = D-U; 
    
    L = theta*Lu;
    
    [~,y] = eig(L);
    m = diag(y);
    lamdaMax = max(m); % 
    
    %%%C
    C1 = zeros(n, c);
    for iv = 1:V
        C1 = C1 + 1.0/pv(iv)*X{iv}'*Mv{iv};
    end
    C2 = zeros(n, c);
    for iv = 1:V
        C2 = C2 + nv(iv)*(Z'-Fv{iv}*Qv{iv});
    end
    
    C = C1+lambda*C2;

    %%%calculate
    F1 = ones(n,c);
    I=eye(n);
    %sign=0;
    while 1
        M = 2*(lamdaMax*I-L)*F+2*C;
        [UU,~,VV] = svds(M);
        F = UU*VV';

        if (norm(F-F1))<1e-3
            break;
        end
        F1 = F;
    end
    
    % Updating Z
    H = Z';
    %%%L 
    sumLGv = zeros(n,n);
    for iv=1:V
        sumLGv = sumLGv+ rv(iv)*LGv{iv};
    end
    L = gamma * sumLGv;
    
    
    [~,y] = eig(L);
    m = diag(y);
    lamdaMax = max(m); % 
    
    %%%C
    C1 = zeros(n, c);
    for iv = 1:V
        C1 = C1 + qv(iv)*X{iv}'*Wv{iv}*A;
    end
    C2 = zeros(n, c);
    for iv = 1:V
        C2 = C2 + nv(iv)*(Fv{iv}*Qv{iv}+F);
    end
    
    C = mu*C1+lambda*C2;%

    %%%calculate
    H1 = ones(n,c);
    I=eye(n);

    while 1
        M = 2*(lamdaMax*I-L)*H+2*C;
        [UU,~,VV] = svds(M);
        H = UU*VV';
        
        %disp(norm(H-H1))
        if (norm(H-H1))<1e-3
            break;
        end

        H1 = H;
    end
    Z = H';
    %disp(Z)
    
    % Updating  Qv
    for iv=1:V
        Qv{iv} = inv(Fv{iv}'*Fv{iv})*(Fv{iv}'*Z'-Fv{iv}'*F);
    end
    
    % Updating nv mv pv
    for iv = 1:V
        hv(iv) = norm(X{iv}'* Mv{iv}-F,'fro');
    end
    for iv = 1:V
        pv(iv) = hv(iv)/sum(hv);
        qv(iv) = 0.5/norm(X{iv}-Wv{iv}*A*Z,'fro');
        mv(iv) = 0.5/norm(Gv{iv}+S-U,'fro');
        nv(iv) = 0.5/norm(Fv{iv}*Qv{iv}+F-Z','fro');
        rv(iv) = 0.5/sqrt(trace(Z*LGv{iv}*Z'));
    end
    
    % calculating obj value
    tempObj = 0;
    
    sumXMF = 0;  % Item 1
    for iv = 1:V
        sumXMF = sumXMF + 1.0 / pv(iv) * norm(X{iv}' * Mv{iv} - F, 'fro')^2;
    end
    tempObj = tempObj + sumXMF;
    
    sumMv = 0;  % Item 2
    for iv = 1:V
        sumMv = sumMv + norm(Mv{iv},'fro')^2;
    end
    tempObj = tempObj + beta * sumMv;
    
    
    sumXWAZ = 0;% item 3
    for iv = 1:V
        sumXWAZ = sumXWAZ + qv(iv)*norm(X{iv}-Wv{iv}*A*Z,'fro')^2;
    end
    tempObj = tempObj + mu*sumXWAZ;
    
    sumGSU = 0;% item 4
    for iv = 1:V
        sumGSU = sumGSU + mv(iv)*norm(Gv{iv}+S-U,'fro')^2;
    end
    tempObj = tempObj + alpha*sumGSU;
    
    D = diag(sum(U));% item 5
    Lu = D-U; 
    tempObj = tempObj + theta*trace(F'*Lu*F);
    
    sumFFZ = 0;  % Item 6
    for iv = 1:V
        sumFFZ = sumFFZ + nv(iv) * norm(Fv{iv}*Qv{iv} + F - Z', 'fro')^2;
    end
    tempObj = tempObj + lambda * sumFFZ;
    
    sumZLGvZ = 0;  % Item 7
    for iv = 1:V
        sumZLGvZ = sumZLGvZ + rv(iv) * trace(Z*LGv{iv}* Z');
    end
    tempObj = tempObj + gamma * sumZLGvZ;
    
    obj(iter) = tempObj;
    fprintf(fp, '%.6f,%.6f,%.6f,%.6f,%.6f,%0.6f,%0.6f\n', beta,mu,alpha,lambda,gamma,theta,obj(iter)); 

    % convergence checking
    if iter>1 && abs(obj(iter)-obj(iter-1))/obj(iter-1) < 1e-3 %
        break;
    end
    
    iter = iter+1;
end
fprintf(fp, '\n');
fclose(fp);
disp(iter);
%% ==========================Plot===========================

%plot(obj+1,'color','b');                     % plot object function value according to running times
%title(dataname);
%xlabel('Nummber of iterations');
%ylabel('Object function value');


end