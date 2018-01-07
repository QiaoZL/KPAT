function [eigvector,Y,eigvalue] = KPCA(X,r,opts)

%

% Kernel Principal Component Analysis

% [eigvector, eigvalue,Y] = KPCA(X,r,opts)

% Input:

% X: d*N data matrix;Each column vector of X is a sample vector.

% r: Dimensionality of reduced space (default: d)

% opts:   Struct value in Matlab. The fields in options that can be set:           

%         KernelType  -  Choices are:

%                  'Gaussian'      - exp{-gamma(|x-y|^2)}

%                  'Polynomial'    - (x'*y)^d

%                  'PolyPlus'      - (x'*y+1)^d

%         gamma       -  parameter for Gaussian kernel

%         d           -  parameter for polynomial kernel

%

% Output:

% eigvector: N*r matrix;Each column is an embedding function, for a new

%            data point (column vector) x,  y = eigvector'*K(x,:)'

%            will be the embedding result of x.

%            K(x,:) = [K(x1,x),K(x2,x),...K(xN,x)]

% eigvalue: The sorted eigvalue of KPCA eigen-problem.

% Y       : Data matrix after the nonlinear transform

X = X';

if nargin<1

  error('Not enough input arguments.')

end

[d,N]=size(X);

if nargin<2

  r=d;
  opts.KernelType = 'Gaussian';
  opts.gamma = 0.5; %1/size(X,1);
  opts.d = 10;

end

%% Ensure r is not bigger than d

if r>d

    r=d;

end;

% Construct the Kernel matrix K
tempY = [];
K = ConstructKernelMatrix(X,tempY,opts);

% Centering kernel matrix

One_N=ones(N)./N;

Kc = K - One_N*K - K*One_N + One_N*K*One_N;

clear One_N;

% Solve the eigenvalue problem N*lamda*alpha = K*alpha

[eigvector, eigvalue] = eig(Kc);

eigvalue = diag(eigvalue);

[junk, index] = sort(-eigvalue);

eigvalue = eigvalue(index);

eigvector = eigvector(:,index);



if r < length(eigvalue)

    eigvalue = eigvalue(1:r);

    eigvector = eigvector(:, 1:r);

end

% maxEigValue = max(abs(eigvalue));

% eigIdx = find(abs(eigvalue)/maxEigValue < 1e-6);

% Normalizing eigenvector

for i=1:length(eigvalue)

    eigvector(:,i)=eigvector(:,i)/sqrt(eigvalue(i));

end;

if nargout >= 3

    % Projecting the data in lower dimensions
    
    Y = eigvector'*K;
    Y = Y';
    % Because the score may be too small to matlab
%     Y = Y * 10000;
    eigvector = eigvector';

end

 


