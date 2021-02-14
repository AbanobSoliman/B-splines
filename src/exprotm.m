function d  = exprotm(R)
% Developed by Abanob SOLIMAN, PhD Student, IBISC Laboratory, France
% Email: abanob.soliman@univ-evry.fr
% Under the supervision of:
% Prof. "Samia Bouchafa Bruneau" <samia.bouchafabruneau@univ-evry.fr>
% Prof. "Dro Désiré Sidibie" <drodesire.sidibie@univ-evry.fr>
% Dr. "fabien bonardi" <fabien.bonardi@univ-evry.fr>
    omega = [R(3,2);R(1,3);R(2,1)];
    
    % Singularity case
    if norm(omega) == 0 
        d = eye(3);
        return
    end
    
    d     = eye(3) + (sin(norm(omega))/norm(omega))*R + ((1-cos(norm(omega)))/(norm(omega)^2))*(R*R);

end