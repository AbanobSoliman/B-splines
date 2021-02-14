function d  = logrotm(R)
% Developed by Abanob SOLIMAN, PhD Student, IBISC Laboratory, France
% Email: abanob.soliman@univ-evry.fr
% Under the supervision of:
% Prof. "Samia Bouchafa Bruneau" <samia.bouchafabruneau@univ-evry.fr>
% Prof. "Dro Désiré Sidibie" <drodesire.sidibie@univ-evry.fr>
% Dr. "fabien bonardi" <fabien.bonardi@univ-evry.fr>
    theta = acos((trace(R)-1)/2);
    
    d     = (theta/(2*sin(theta)))*(R-R');

end