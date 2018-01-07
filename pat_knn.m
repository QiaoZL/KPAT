function [index_vals,vector_vals,final_nodes,calculationtimes] = pat_knn(tree,point,k,plot_stuff,node_number,b,dlbpower)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% tree        --- the cell array that contains the tree
% point       --- the point of interest
% node_number --- Internal Variable, Donot Define

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
% index_vals  --- the index value in the original matrix that was used to
%                 build the tree
% vector_vals --- the feature vector closest to the given vector
% final_node  --- the node that contains the closest vector, now I return
% the distance between result and query to it.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Donot define the node_number variable -- it is used for internal
% referencin

global best_points_mat
global number_of_points
global tree_cell
global safety_check;

global calculationtime;




dim =size(point,2);
index_vals=0;
vector_vals=0;
final_nodes=0;
debug_val=plot_stuff;

if(debug_val); global h; end

if (nargin==4)
    safety_check=0;
    % for index in tree, index in original mat, dist to point
    best_points_mat=zeros(k,dim+1+1+1);
    number_of_points=0;
    node_number=1;
    b = point;
    dlbpower = 0;
    
    calculationtime = 0;
    
    if(debug_val)
        h = plot(best_points_mat(1:k,1),best_points_mat(1:k,2),'g*');
    end;
    tree_cell=tree;
    clear tree;
end

if (isempty(safety_check))
    error ('Insufficient number of input variables ... please check ');
end



% if the current node is a leaf then output its results
if(strcmp(tree_cell(node_number).type,'leaf'))
    node_check(point,k,node_number,debug_val);
    return;
end

% calculate the projection of boundary point on principal axis
% initialize the sign of search
p = tree_cell(node_number).principalaxis;
projection = b * p;
lowerdone = 'false';
upperdone = 'false';



% if the current node is not a leaf
% check to see if the point is to the left of the split dimension
% if it is to the left then recurse to the left
if (projection <= tree_cell(node_number).hypervector(1,2))
    lowerdone = 'true';
    if (isempty(tree_cell(node_number).left))
        % incase the left node is empty, then output current results
    else
        pat_knn(0,point,k,plot_stuff,tree_cell(node_number).left,b,dlbpower);
    end
elseif (projection > tree_cell(node_number).hypervector(1,2))
    % as the point is to the right of the split dimension
    % recurse to the right
    upperdone = 'true';
    if (isempty(tree_cell(node_number).right))
        % incase the right node is empty, then output current results
    else
        pat_knn(0,point,k,plot_stuff,tree_cell(node_number).right,b,dlbpower);
    end
end




% do the computation to decide if you need to search on the otherside
% initialize variables for main loop
dislbipowerl = 0;
dislbipoweru = 0;
disu=0;
disl=0;
distk=best_points_mat(k,dim+3);
search_otherside= 'true';

if(tree_cell(node_number).index == 39)
    wtf = 1;
end

if (number_of_points<k)
    % if still have not enough result, the otherside must be chec
    
else
    
    if(strcmp(lowerdone,'false'))
        disl = projection - tree_cell(node_number).hypervector(2,1);
    end
    if(strcmp(upperdone,'false'))
        disu = tree_cell(node_number).hypervector(1,2) - projection;
    end
    
    if(strcmp(upperdone,'true')&&strcmp(lowerdone,'false'))
        dislbipowerl = dlbpower + disl^2;
        if (distk < dislbipowerl)
            search_otherside= 'false';
        end
    end
    if(strcmp(lowerdone,'true')&&strcmp(upperdone,'false'))
        dislbipoweru = dlbpower + disu^2;
         if (distk < dislbipoweru)
            search_otherside= 'false';            
         end
    end
    
    
    
    
end



if (strcmp(search_otherside,'true'))
    
    node_check(point,k,node_number,debug_val);
    
    % if the current node is not a leaf
    % check to see if the point is to the left of the split dimension
    % if it is to the left then decide whether to recurse to the right
    if (projection <= tree_cell(node_number).hypervector(1,2))
        if (isempty(tree_cell(node_number).right))
            % incase the right node is empty, then output current results
            
        else
            % as the point is to the right of the split dimension
            % decide whether to recurse to the left
            if (disu ~= 0)
                b = b + (disu*p)';
            end
            pat_knn(0,point,k,plot_stuff,tree_cell(node_number).right,b,dislbipoweru);
        end
    elseif (projection > tree_cell(node_number).hypervector(1,2))
        if (isempty(tree_cell(node_number).left))
            % incase the left node is empty, then output current results
            
        else
            % as the point is to the right of the split dimension
            % decide whether to recurse to the left
            if (disl ~= 0)
                b = b - (disl*p)';
            end
            pat_knn(0,point,k,plot_stuff,tree_cell(node_number).left,b,dislbipowerl);
        end
    end
    
end

if (nargin==4)
    
    vector_vals=best_points_mat(1:number_of_points,1:dim);
    index_vals=best_points_mat(1:number_of_points,dim+1);
    final_nodes=best_points_mat(1:number_of_points,dim+3);
    
    calculationtimes = calculationtime;
    
    clear global best_points_mat;
    clear global number_of_points;
    clear global tree_cell;
    clear global safety_check;
    
end


function []=node_check(point,k,node_number,debug_val)

global best_points_mat
global number_of_points
global tree_cell
global calculationtime

if(debug_val);
    global h;
end

calculationtime = calculationtime + 1;

dim =size(point,2);
distance = sum((point-tree_cell(node_number).nodevector).^2);

if(tree_cell(node_number).index == 39)
    wtf = 1;
end

if (number_of_points==k && best_points_mat(k,dim+3)>distance)
    best_points_mat(k,1:dim)=tree_cell(node_number).nodevector;
    best_points_mat(k,dim+1)=tree_cell(node_number).index;
    best_points_mat(k,dim+2)=node_number;
    best_points_mat(k,dim+3)=distance;
    best_points_mat=sortrows(best_points_mat,dim+3);
    if(debug_val);
        set(h,'XData',best_points_mat(1:k,1))
        set(h,'YData',best_points_mat(1:k,2))
    end
elseif(number_of_points<k)
    number_of_points=number_of_points+1;
    best_points_mat(number_of_points,1:dim)=tree_cell(node_number).nodevector;
    best_points_mat(number_of_points,dim+1)=tree_cell(node_number).index;
    best_points_mat(number_of_points,dim+2)=node_number;
    best_points_mat(number_of_points,dim+3)=distance;
    % once the variable gets filled up then sort the rows
    if(number_of_points==k)
        best_points_mat=sortrows(best_points_mat,dim+3);
    end
    if(debug_val);
        set(h,'XData',best_points_mat(1:k,1))
        set(h,'YData',best_points_mat(1:k,2))
    end
end

return;
