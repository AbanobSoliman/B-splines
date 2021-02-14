function S = b_splineR3(Poses,t,k)
% Developed by Abanob SOLIMAN, PhD Student, IBISC Laboratory, France
% Email: abanob.soliman@univ-evry.fr
% Under the supervision of:
% Prof. "Samia Bouchafa Bruneau" <samia.bouchafabruneau@univ-evry.fr>
% Prof. "Dro Désiré Sidibie" <drodesire.sidibie@univ-evry.fr>
% Dr. "fabien bonardi" <fabien.bonardi@univ-evry.fr>
    S = [];
    n = size(Poses,2) - 1;

    if (k==2)

        disp('Quadratic B-spline!')
        M = 0.5.*[1 -2 1;-2 2 1;1 0 0]; % Basis function matrix
        T = [t.^2;t.^1;t.^0];

        for i = 1:n-k+1
            B = Poses(:,i:i+k)*M*T; % B-spline segment
            S = [S B];
        end

        return

    elseif (k==3)

        disp('Cubic B-spline!')
        M = (1/6).*[-1 3 -3 1;3 -6 3 0;-3 0 3 0;1 4 1 0]; % Basis function matrix
        T = [t.^3;t.^2;t.^1;t.^0]';

        for i = 1:n-k+1
            B = T*M*Poses(:,i:i+k)'; % B-spline segment
            S = [S;B];
        end
        
        S = S';

        return

    else

        disp('This function works with quadratic and cubic B-splines!')
        return 

    end    

end