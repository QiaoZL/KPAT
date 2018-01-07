function K = ConstructKernelMatrix(X,Y,opts)

%

% function K=ConstructKernelMatrix(X,Y,opts)

% Usage:

%   opts.KernelType='Gaussian';

% K = ConstructKernelMatrix(X,[],opts)

%   K = ConstructKernelMatrix(X,Y,opts)

%

% Input:

% X: d*N data matrix;Each column vector of X is a sample vector.

% Y: d*M data matrix;Each column vector of Y is a sample vector.

% opts:   Struct value in Matlab. The fields in options that can be set:                

%         KernelType  -  Choices are:

%                  'Gaussian'      - exp{-gamma(|x-y|^2)}

%                  'Polynomial'    - (x'*y)^d

%                  'PolyPlus'      - (x'*y+1)^d

%         gamma       -  parameter for Gaussian kernel

%         d           -  parameter for polynomial kernel

% Output:

% K N*N or N*M matrix

if nargin<1

  error('Not enough input arguments.')

end

if (~exist('opts','var'))

   opts = [];

else

   if ~isstruct(opts)

       error('parameter error!');

   end

end

N=size(X,2);

if isempty(Y)

    K=zeros(N,N);

else

    M=size(Y,2);

    if size(X,1)~=size(Y,1)

        error('Matrixes X and Y should have the same row dimensionality!');

    end;

    K=zeros(N,M);

end;

%=================================================

if ~isfield(opts,'KernelType')

    opts.KernelType = 'Gaussian';   

end

switch lower(opts.KernelType)

    case {lower('Gaussian')}        %  exp{-gamma(|x-y|^2)}

        if ~isfield(opts,'gamma')

            opts.gamma = 0.5;

        end

    case {lower('Polynomial')}      % (x'*y)^d

        if ~isfield(opts,'d')

            opts.d = 1;

        end

    case {lower('PolyPlus')}      % (x'*y+1)^d

        if ~isfield(opts,'d')

            opts.d = 1;

        end      

    otherwise

       error('KernelType does not exist!');

end

switch lower(opts.KernelType)

    case {lower('Gaussian')}      

        if isempty(Y)

            for i=1:N

               for j=i:N

                   dist = sum(((X(:,i) - X(:,j)).^2));

                    temp=exp(-opts.gamma*dist);

                    K(i,j)=temp;

                    if i~=j

                        K(j,i)=temp;

                    end;

                end

            end

        else

            for i=1:N

               for j=1:M

                    dist = sum(((X(:,i) - Y(:,j)).^2));

                    K(i,j)=exp(-opts.gamma*dist);                  

                end

            end

        end      

    case {lower('Polynomial')}    

        if isempty(Y)

            for i=1:N

                for j=i:N                   

                    temp=(X(:,i)'*X(:,j))^opts.d;

                    K(i,j)=temp;

                    if i~=j

                        K(j,i)=temp;

                    end;

                end

            end

        else

            for i=1:N

                for j=1:M                                      

                    K(i,j)=(X(:,i)'*Y(:,j))^opts.d;

                end

            end

        end      

    case {lower('PolyPlus')}    

        if isempty(Y)

            for i=1:N

                for j=i:N                   

                    temp=(X(:,i)'*X(:,j)+1)^opts.d;

                    K(i,j)=temp;

                    if i~=j

                        K(j,i)=temp;

                    end;

                end

            end

        else

            for i=1:N

                for j=1:M                                      

                    K(i,j)=(X(:,i)'*Y(:,j)+1)^opts.d;

                end

            end

        end      

    otherwise

        error('KernelType does not exist!');

end