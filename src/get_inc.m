function dX = get_inc(X)
% Developed by Abanob SOLIMAN, PhD Student, IBISC Laboratory, France
% Email: abanob.soliman@univ-evry.fr
% Under the supervision of:
% Prof. "Samia Bouchafa Bruneau" <samia.bouchafabruneau@univ-evry.fr>
% Prof. "Dro Désiré Sidibie" <drodesire.sidibie@univ-evry.fr>
% Dr. "fabien bonardi" <fabien.bonardi@univ-evry.fr>    
    if isa(X,'double')
        
        for i = 1 : size(X,2)-1
            dX(:,i) = X(:,i+1) - X(:,i);
        end
        
        return
        
    elseif isa(X,'quaternion')
        
        for i = 1 : length(X)-1
            dX(i,:) = quatmultiply( compact(X(i)) , quatinv(compact(X(i+1))) );
        end
        
        dX = dX';
        
        return
        
    else
        
        disp('This function gets increments of positions(doubles) and quaternions!');
        return 
        
    end
    
end