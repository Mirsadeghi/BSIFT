function [Y,dim]=comp_pca(X,ratio)

[~,N]=size(X);

% Subtracting the mean
mean_vec=mean(X')';
X_zero=X-mean_vec*ones(1,N);

% Computing the covariance matrix and its eigenvalues/eigenvectors
R=nancov(X_zero');

[V,D]=eig(R);

eigenval=diag(D);
[eigenval,ind]=sort(eigenval,1,'descend');
eigenvec=V(:,ind);

explain=eigenval/sum(eigenval);
tmp = cumsum(explain);
dim = sum(tmp<=ratio);
% Keeping the first m eigenvaules/eigenvectors
eigenvec=eigenvec(:,1:dim);

% Computing the transformation matrix
A=eigenvec(:,1:dim)';

% Computing the transformed data set
Y=A*X;