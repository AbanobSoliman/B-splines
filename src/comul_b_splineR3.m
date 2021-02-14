function S = comul_b_splineR3(P,u,N)
% Developed by Abanob SOLIMAN, PhD Student, IBISC Laboratory, France
% Email: abanob.soliman@univ-evry.fr
% Under the supervision of:
% Prof. "Samia Bouchafa Bruneau" <samia.bouchafabruneau@univ-evry.fr>
% Prof. "Dro Désiré Sidibie" <drodesire.sidibie@univ-evry.fr>
% Dr. "fabien bonardi" <fabien.bonardi@univ-evry.fr>
    S = [];
    U = [];
    M = [];
    m = [];
    add = [];
    
    k = N + 1;
    dP = get_inc(P);
    
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

    for i = 1 : size(P,2)-k+1       
        B = [P(:,i),dP(:,i:i+k-2)]*M*U; % B-spline segment            
        S = [S B];          
    end

end