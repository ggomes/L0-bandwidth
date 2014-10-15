classdef graph
    
    properties (SetAccess = public, GetAccess = public)
        max_node_id
        max_edge_id
        nodes
        edges
    end
    
    methods ( Access = public )
        
        function obj=graph()
            obj.max_node_id = 0;
            obj.max_edge_id = 0;
            obj.nodes = [];
            obj.edges = [];
        end
        
        function obj=add_nodes(obj,n)
            if(iscell(n))
                n = [n{:}];
            end
            for i=1:length(n)
                n(i).id = obj.max_node_id;
                obj.max_node_id = obj.max_node_id+1;
            end
            
            %%% should check that there are no other nodes at same position
            
            obj.nodes =  [obj.nodes reshape(n,1,length(n))];
        end
        
        function obj=add_edges(obj,E)

            % assign start and end ids to edges
            for i=1:length(E)
                E(i).node_id_start = obj.get_node_id_at(E(i).p_start);
                E(i).node_id_end = obj.get_node_id_at(E(i).p_end);
            end
            
            % assign id
            for i=1:length(E)
                E(i).id = obj.max_edge_id;
                obj.max_edge_id = obj.max_edge_id+1;
            end            
            
            % inform nodes
            node_ids = [obj.nodes.id];
            for i=1:length(E)
                ind = E(i).node_id_start==node_ids;
                obj.nodes(ind).edges_out = [obj.nodes(ind).edges_out E(i).id];
                
                ind = E(i).node_id_end==node_ids;
                obj.nodes(ind).edges_in = [obj.nodes(ind).edges_in E(i).id];
            end
            
            % add to list
            obj.edges =  [obj.edges reshape(E,1,length(E))];            
            
        end
        
        function [id]=get_node_id_at(obj,p)
            ind = p(1)==[obj.nodes.pos_x] & p(2)==[obj.nodes.pos_y];
            if(~any(ind))
                error(['no node found at: (' num2str(p(1)) ',' num2str(p(2)) ')'])
            end
            id = obj.nodes(ind).id;
        end
        
        function [edge_ind]=get_edge_between(obj,n1,n2)
            edge_id = intersect(obj.nodes(n1).edges_out,obj.nodes(n2).edges_in);
            edge_ind = edge_id==[obj.edges.id];
        end
        
        function [f]=plot(obj,f,node_path,highlight,labelon,hlinecolor)
            
            addpath([docroot '/techdoc/creating_plots/examples'])
            
            if(nargin<5)
                labelon = false;
            end
            
            if(nargin<4)
                highlight = false;
            end
            
            if(nargin<3)
                t_nodes = obj.nodes;
                t_edges = obj.edges;
            else
                t_nodes = obj.nodes(node_path);
                for i=1:length(node_path)-1
                    e = obj.get_edge_between(node_path(i),node_path(i+1));
                    t_edges(i) = obj.edges(e);
                end
                clear e
            end
            
            if(nargin<2)
                f=figure;
            end
            
            figure(f)
            
            if(highlight)
                dotmarker = 'o';
                dotsize = 12;
                linewidth = 2;
                if(nargin<6)
                    linecolor = 'm';
                else
                    linecolor = hlinecolor;
                end
            else
                dotmarker = '.';
                dotsize = 14;
                linewidth = 1;
                linecolor = 'k';
            end
                
            X = [t_nodes.pos_x];
            Y = [t_nodes.pos_y];
            
            plot(X,Y,'Marker',dotmarker,'LineStyle','none','MarkerSize',dotsize,'Color',linecolor,'LineWidth',linewidth)
            
            if(labelon)
                for i=1:length(t_nodes)
                    labels{i} = [t_nodes(i).name ' (' num2str(t_nodes(i).id) ')'];
                end
                text(X,Y,labels,'VerticalAlignment','bottom', 'HorizontalAlignment','right')
            end
            
            hold on
            axis(axis)
            
            for i=1:length(t_edges)
                e = t_edges(i);
                bbox = dsxy2figxy(gca, [e.p_start e.p_end-e.p_start]);
                x = bbox(1) + [0 bbox(3)];
                y = bbox(2) + [0 bbox(4)];
                %h=annotation('line',x,y);
                h=annotation('arrow',x,y);
                set(h,'LineWidth',linewidth,'Color',linecolor);
            end
                        
            set(gca,'XTickLabel','')
            set(gca,'YTickLabel','')
            set(gca,'XTick',[])
            set(gca,'YTick',[])
        end
        
        function []=highlight_path(obj,f,path,c)
            hold on
            obj.plot(f,path,true,false,c);
        end
      
    end
    
end