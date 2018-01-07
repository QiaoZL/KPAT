
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Just a simple demo to show the functionality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
demo_case = [1,2,3,4];

plot_stuff=0;   % 1 if you want to plot the data
% change to 0 if you donot want to plot anything.

if (plot_stuff) close all; end
clc

size_of_query = 500;
rand('seed',1)
dimen=5;
X=rand(3000,dimen);
rand('seed',2)
Y=rand(size_of_query,dimen);
num_of_points=1;
point=0.6*ones(1,dimen);


sum_time = 0;

disp('##### Build Tree #####');

tree = pat_buildtree(X,plot_stuff);

disp('##### Build Done #####');
for i = 1:size_of_query
    disp('##### N Closest Points #####');
    point=X(i,1:dimen)
%     point=0.6*ones(1,dimen);
    if (plot_stuff); hold on ; end
    if (plot_stuff); plot(point(1),point(2),'g*','MarkerSize',10); end
    
    
    [index_vals,vec_vals, dist_vals, calculationtimes]  = pat_knn(tree,point,num_of_points,plot_stuff);
    index_vals
    calculationtimes
    sum_time = sum_time + calculationtimes;
    if (plot_stuff);
        plot(X(index_vals,1),X(index_vals,2),'g*');
        dist=sqrt(sum((point-X(index_vals(end),1:2)).^2));
        plot(point(1)+dist*cos(0:0.1:2*pi),point(2)+dist*sin(0:0.1:2*pi),'g-','LineWidth',2)
    end
    
end
sum_time



