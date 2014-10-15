function [path,pathlength]=dijkstra(G,start_id,end_id,shortorlong)

if(nargin<3)
    shortorlong = 'shortest';
end
    
% graph information
node_ids = [G.nodes.id];
edge_ids = [G.edges.id];
numnodes = length(G.nodes);
numedges = length(G.edges);

% start and end nodes
start_ind = start_id==node_ids;
end_ind = end_id==[G.nodes.id];

% checks
if(~any(start_ind) || ~any(end_ind))
    error('unknown node')
end

% initialize
dist = inf(1,numnodes);         % distance function
dist(start_ind) = 0;
previous = nan(1,numnodes);     % Previous node in optimal path from source                                              %
inQ = true(1,numnodes);         % All nodes in the graph are unoptimized - thus are in Q
dist_between = create_distance_matrix();

if(strcmp(shortorlong,'longest'))
    dist_between = -dist_between;
end

while inQ(end_ind) %any(inQ)
        
    [~,u_ind] = min(dist(inQ));  % node with minimum distance
    
    u = false(1,numnodes);
    inQ_ind = find(inQ);
    u(inQ_ind(u_ind)) = true;    % convert to index array
    clear u_ind inQ_ind
        
    inQ(u) = false;         % remove u from Q
    
    if isinf(dist(u))       % all remaining vertices are inaccessible from source
        break                  
    end
    
    N = get_neighbors(u);
    N(~inQ) = false;        % remove neighbors not in Q
    
    for v_ind=find(N) % each neighbor v of u where v has not yet been removed from Q.

        v = false(1,numnodes);
        v(v_ind) = true;            % convert to index array
        
        alt = dist(u) + dist_between(u,v);
        if alt < dist(v)
            dist(v) = alt;
            previous(v) = find(u);
        end
    end
end

% find path
start_index = find(start_ind);
pathlength = abs(dist(end_ind));
c = find(end_ind);
path = c;
while(c~=start_index)
    c = previous(c);
    path = [c path];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [NN]=get_neighbors(node_ind)
        node = G.nodes(node_ind);
        num_edge_out = length(node.edges_out);
        NN = false(1,numnodes);
        for ii=1:num_edge_out
            ee = G.edges(node.edges_out(ii)==edge_ids);
            nnode_ind = ee.node_id_end==node_ids;
            NN(nnode_ind) = true;
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [DD] = create_distance_matrix()
        DD = inf(numnodes,numnodes);
        for ii=1:numedges
            ee = G.edges(ii);
            nnstart = ee.node_id_start==node_ids;
            nnend = ee.node_id_end==node_ids;
            DD(nnstart,nnend) = ee.val;
        end        
    end

        
end




