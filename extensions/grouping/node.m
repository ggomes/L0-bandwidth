classdef node
    
    properties (SetAccess = public, GetAccess = public)
        id
        name
        pos_x
        pos_y
        edges_in
        edges_out
    end
    
    methods 
        
        function obj=node(tname,tpos)
            obj.id = nan;
            obj.name = tname;
            obj.pos_x = tpos(1);
            obj.pos_y = tpos(2);
        end

    end
    
    
end