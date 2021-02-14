function S = comul_b_splineSE3(t,u,N)
% Developed by Abanob SOLIMAN, PhD Student, IBISC Laboratory, France
% Email: abanob.soliman@univ-evry.fr
% Under the supervision of:
% Prof. "Samia Bouchafa Bruneau" <samia.bouchafabruneau@univ-evry.fr>
% Prof. "Dro Désiré Sidibie" <drodesire.sidibie@univ-evry.fr>
% Dr. "fabien bonardi" <fabien.bonardi@univ-evry.fr>
    A = [];
    S = [];
    U = [];
    M = [];
    m = [];
    add = [];
    
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
    for i = 1 : size(t,2)-k+1
        A  = [quat2rotm(t(4:end,i)'),t(1:3,i);0 0 0 1];
        for pr = 1 : length(u)
            for j = 1 : k-1
                T1_inv = [quat2rotm(t(4:end,i+j-1)')',-quat2rotm(t(4:end,i+j-1)')'*t(1:3,i+j-1);0 0 0 1];
                T1 = [quat2rotm(t(4:end,i+j-1)'),t(1:3,i+j-1);0 0 0 1];
                T2 = [quat2rotm(t(4:end,i+j)'),t(1:3,i+j);0 0 0 1];
                D  = logrotmT(T1_inv*T2);
                A  = A*exprotmT(Bt(j+1,pr)*D);
            end
            Aq = [A(1:3,4);quatnormalize(rotm2quat(A(1:3,1:3)))'];
            S = [S Aq]; % B-spline segment 
            A = [quat2rotm(t(4:end,i)'),t(1:3,i);0 0 0 1];
        end
    end
    
end