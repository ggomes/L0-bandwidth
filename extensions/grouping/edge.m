classdef edge
    
    properties (SetAccess = public, GetAccess = public)
        id
        p_start
        p_end
        node_id_start
        node_id_end
        val
    end
    
    methods 
        
        function obj = edge(pstart,pend,tval)
            obj.id = nan;
            obj.p_start = pstart;
            obj.p_end = pend;
            obj.node_id_start = nan;
            obj.node_id_end = nan;
            obj.val = tval;           
        end

    end
       
end