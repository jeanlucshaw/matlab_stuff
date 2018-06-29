function [vals,vecs,stds] = jlpca(X)
% Synthax :            [vals,vecs,stds] = jlpca(X)
% 
% Performs principal component analysis on input matrix 'X' . Input matrix
% must be in a form where columns are different features of the data set
% and rows are repetitions of the same measurement.
%
% Also returns the standard deviation of every feature projected in the space
% of the eigenvectors found by the PCA.
%


[m,n] = size(X) ;


%% De-trend
F      = repmat(nan,size(X)) ;
for ii = 1:n
    F(:,ii) = X(:,ii) - mean(X(:,ii),'omitnan') ;
end

%% Build covariance matrix
R = F'*F ;

%% Find eigenvalues
[C,L]     = eig(R) ;

%% sort and normalise
vals      = sort(diag(L),'descend') ;

% eigenvectors in growing importance
vecs = fliplr(C) ;

%% Project in the space of the eigenvectors
Fprime = inv(vecs)*F' ;

%% Get the standard error in each "eigendimension"
stds   = std(Fprime,[],2) ;

end
