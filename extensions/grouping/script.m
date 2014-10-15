clear
%close all

K = 5;
I = 15;

G = graph();

clear all_nodes 
all_edges = [];
for k=0:K
    
    if(k==0)
       all_nodes{1,1} = node('1',[0,1]); % initial node
       clear E
       for i=2:I-K+1
           E(i-1) = edge([0,1],[1,i-1],1);
       end
       all_edges = [all_edges E];
       continue
    end
    
    if(k>0 && k<K)
        for i=1:I-K
            all_nodes{k+1,i} = node(num2str(i+k),[k,i]);
            clear E
            c = 0;
            for ii=1:I-K-i+1
                c = c+1;
                E(c) = edge([k,i],[k+1,i+ii-1],1);  
            end
            if(~isempty(ii))
                all_edges = [all_edges E];
            end
        end        
        continue
    end
    
    if(k==K)
        all_nodes{K+1,I-K} = node(num2str(I),[K,I-K]);
        
        % remove edges not from K-1, not arriving at I
        p_start = reshape([all_edges.p_start],2,length(all_edges));
        p_end = reshape([all_edges.p_end],2,length(all_edges));
        remove_these = p_start(1,:)==K-1 & p_end(2,:)~=I-K;
        all_edges(remove_these)=[];
        clear p_start remove_these
        continue
    end
    
end
G = G.add_nodes(all_nodes');
G = G.add_edges(all_edges);

% put some weights on the edges. Pick euclidean length of the edge
for i=1:length(G.edges)
    e = G.edges(i);
    G.edges(i).val = (sum((e.p_end-e.p_start).^2))^(0.3);
end
clear e

% run shortest path for each one
start_node = G.get_node_id_at([0,1]);
end_node = G.get_node_id_at([K,I-K]);
[s_path,s_L]=dijkstra(G,start_node,end_node,'shortest');
[l_path,l_L]=dijkstra(G,start_node,end_node,'longest');
f=G.plot();
G.highlight_path(f,s_path,'m');
G.highlight_path(f,l_path,'b');
title(['shortest: ' num2str(s_L) ', longest: ' num2str(l_L)])

