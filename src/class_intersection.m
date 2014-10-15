classdef class_intersection < handle
    %CLASS_INTERSECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties( Access = public )
        myArtery    % reference to class_artery
        
        name        % name of the intersection
        
        % parameters for pretimed intersection
        gsplit_o    % [portion of cycle]
        gsplit_i    % [portion of cycle]
        deltasplit  % [portion of cycle]
        
        % pretimed green window
        go          % [sec]
        gi          % [sec]
        
        % gaussian green window
        sigma_o     % [sec]    
        sigma_i     % [sec]
        gamma_o     % [-]
        gamma_i     % [-]
        
    end
 
    properties( Access = public )
        absoffseto      % [sec]
        absoffseti      % [sec]
        reloffseto      % [sec]
        reloffseti      % [sec]
        delta           % [sec] gi_midpoint - go_midpoint
    end
    
    methods( Access = public )
        
        function [obj]=class_intersection(varargin)
            % name = varargin{1}
            % windowtype = varargin{2}
            % for pretimed:
            %     gsplit_o = varargin{3};
            %     gsplit_o = varargin{4};
            %     deltasplit = varargin{5};
            % for actuated gaussian: (not implemented)
            
            if(nargin>=1)
                obj.name = varargin{1};
            end
            
            if(nargin>=2)
                if(~any(strcmp(varargin{2},{'pretimed','gaussian'})))
                    error('unknown intersection type')
                end
            end
            
            if(nargin==5)

                if(any([varargin{3:4}]<0) || any([varargin{3:4}]>1))
                    error('green splits must be in [0,1]')
                end

                if(any([varargin{5}]<-0.5) || any([varargin{5}]>0.5))
                    error('delta split must be in [-0.5,0.5]')
                end

                obj.gsplit_o = varargin{3};
                obj.gsplit_i = varargin{4};
                obj.deltasplit = varargin{5};

                if( 2*abs(obj.deltasplit) >1 )
                    error('delta should be in [-0.5,0.5].')
                end

                if( abs(obj.deltasplit)+(obj.gsplit_o+obj.gsplit_i)/2 >1 )
                    error('splits and delta add up to more than 1.')
                end
                
            end
        end
        
        function [obj]=initialize(obj,windowtype)
            
            obj.delta = obj.deltasplit * obj.myArtery.cycle;
            
            switch windowtype
                case 'pretimed'
                    obj.go    = obj.gsplit_o   * obj.myArtery.cycle;
                    obj.gi    = obj.gsplit_i   * obj.myArtery.cycle;
                case 'gaussian'
                    xgo = obj.gsplit_o * obj.myArtery.cycle;
                    xgi = obj.gsplit_i * obj.myArtery.cycle;
                    obj.sigma_o = xgo/sqrt(2*pi);
                    obj.sigma_i = xgi/sqrt(2*pi);
                    obj.gamma_o = xgo;
                    obj.gamma_i = xgi;
            end
        end
       
        function []=set_reloffset(obj,offo,offi)
            offo = obj.myArtery.modhalf(offo);
            offi = obj.myArtery.modhalf(offi);
%             if(abs(obj.myArtery.modhalf(offi-offo)-obj.delta)<1e-6)
                obj.reloffseto = offo;
                obj.reloffseti = offi;
%             else
%                 warning('set offset failed')
%             end
        end
        
        function []=set_absoffset(obj,offo,offi)
            offo = obj.myArtery.modhalf(offo);
            offi = obj.myArtery.modhalf(offi);
            if(abs(obj.myArtery.modhalf(offi-offo)-obj.delta)<1e-6)
                obj.absoffseto = offo;
                obj.absoffseti = offi;
            else
                obj.absoffseto = nan;
                obj.absoffseti = nan;
                warning('set offset failed')
            end
        end
        
    end
    
end

