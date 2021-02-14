function S = comul_b_splineSO3(Q,u,N)
% Developed by Abanob SOLIMAN, PhD Student, IBISC Laboratory, France
% Email: abanob.soliman@univ-evry.fr
% Under the supervision of:
% Prof. "Samia Bouchafa Bruneau" <samia.bouchafabruneau@univ-evry.fr>
% Prof. "Dro Désiré Sidibie" <drodesire.sidibie@univ-evry.fr>
% Dr. "fabien bonardi" <fabien.bonardi@univ-evry.fr>
    A = [1,0,0,0];
    S = [];
    U = [];
    M = [];
    m = [];
    add = [];
    
    Q  = compact(Q);
    k  = N + 1;
    
    s = 0:1:k-1;
    n = 0:1:k-1;
    L = 0:1:k-1;

    for i = 1 : k % column
        for j = 1 : k % row
            for l = i : k
                add(l) = (-1)^(L(l)-s(i))*nchoosek(k,L(l)-s(i))*(k-1-L(l))^(k-1-n(j));
            end
            m(i,j) = (nchoosek(k-1,n(j))/(factorial(k-1)))*sum(add(:)); % Basis function matrix
            add(:) = 0;
        end
    end

    for j = 1 : k % column
        for n = 1 : k % row
            M(j,n) = sum(m(j:k,n)); % Comulative Basis function matrix
        end
    end

    for i = 1 : k
        U(i,:) = u.^(i-1); 
    end

    Bt = M*U; % k x pr terms
    for i = 1 : size(Q,1)-k+1
        A  = quat2rotm(Q(i,:));
        for pr = 1 : length(u)
            for j = 1 : k-1
                R1 = quat2rotm(Q(i+j-1,:));
                R2 = quat2rotm(Q(i+j,:));
                d  = logrotm(R1'*R2);
                A  = A*exprotm(Bt(j+1,pr)*d);
            end
            Aq = quatnormalize(rotm2quat(A))';
            S = [S Aq]; % B-spline segment 
            A = quat2rotm(Q(i,:));
        end
    end
    
end