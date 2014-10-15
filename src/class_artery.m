classdef class_artery < handle
    %CLASS_ARTERY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties( Access = public )
        n               % [-] number of intersections
        intersection    % 1xn array of class_intersection
        x               % [ft] 1xn intersection positions
        cycle           % [sec] cycle time for coordination
        vi              % [fps] (n-1)x1 inbound speeds
        vo              % [fps] (n-1)x1 outbound speeds
        optbandwidth    % [sec] optimal total bandwidth
        optbo           % [sec] optimal outbound bandwidth
        optbi           % [sec] optimal inbound bandwidth
        goodparameters  % [true/false] parameters validate
        windowtype      % pretimed or gaussian
        pretimed_optimization_method    % 'simplex' or 'ip' (interiot point), or milp
        
        milp_m          % optimization variable for milp
        milp_w          % optimization variable for milp
        milp_wbar       % optimization variable for milp
    end
    
    methods( Access = public )
        
        function [obj]=class_artery(c,windowtype)
            obj.n = 0;
            obj.x = 0;
            obj.intersection = [];
            obj.goodparameters = false;
            if(nargin>=1)
                obj.cycle = c;
            else
                obj.cycle = 100;
            end
            if(nargin>=2)
                obj.windowtype = windowtype;
            else
                obj.windowtype = 'pretimed';
            end
            if(strcmp(obj.windowtype,'pretimed'))
                obj.pretimed_optimization_method = 'ip';
            end
        end
        
        % change the cycle but keep the gren times -> recompute split values
        function [obj]=setCycleKeepGreens(obj,newcycle)
            
            % compute smallest acceptable cycle
            mincycle = -inf;
            for i=1:obj.n
                I = obj.intersection(i);
                mincycle = max([ mincycle (I.go+I.gi)/2 + abs(I.delta) ]);
                mincycle = max([ mincycle 2*abs(I.delta) ]);
            end
            
            if(newcycle<mincycle)
                disp(['Failed: the minimum allowable cycle is ' num2str(mincycle)])
                obj.goodparameters = false;
                return
            end
            
            for i=1:obj.n
                obj.intersection(i).gsplit_o   = obj.intersection(i).go / newcycle;
                obj.intersection(i).gsplit_i   = obj.intersection(i).gi / newcycle;
                obj.intersection(i).deltasplit = obj.intersection(i).delta / newcycle;
            end
            obj.cycle = newcycle;
            obj.goodparameters = true;
        end
        
        % change the cycle but keep the green splits -> recompute absolute values
        function [obj]=setCycleKeepSplits(obj,newcycle)
            obj.cycle = newcycle;
            for i=1:obj.n
                obj.intersection(i).initialize(obj.windowtype);
            end
        end
        
        function []=add_intersection(obj,I)
            if(length(obj.intersection)~=length(obj.vi))
                disp('add a segment first')
            end
            I.myArtery = obj;
            obj.n = obj.n + 1;
            obj.intersection = [obj.intersection I];
        end
        
        % removes intersection and joins segments
        function []=remove_intersection(obj,name)
            ind = strcmp(name,{obj.intersection.name});
            if(~any(ind))
                return
            end
            i = find(ind);
            if(i==1)                % remove first segment
                obj.vi(1) = [];
                obj.vo(1) = [];
                obj.x = obj.x-obj.x(2);
                obj.x(1) = [];
            elseif(i==obj.n)        % remove last segment
                obj.vi(end) = [];
                obj.vo(end) = [];
                obj.x(end) = [];
            else                            % join two segments
                L1 = obj.x(i) - obj.x(i-1);
                L2 = obj.x(i+1) - obj.x(i);
                vi1 = obj.vi(i-1);
                vi2 = obj.vi(i);
                vo1 = obj.vo(i-1);
                vo2 = obj.vo(i);
                obj.vi(i-1) = (L1+L2)/(L1/vi1 + L2/vi2);
                obj.vi(i-1) = (L1+L2)/(L1/vo1 + L2/vo2);
                obj.vi(i) = [];
                obj.vo(i) = [];
                obj.x(i) = [];
            end
            
            obj.n = obj.n-1;
            obj.intersection(ind)=[];
            
        end
        
        function []=add_segment(obj,L,vin,vout)
            if(length(obj.intersection)~=length(obj.vi)+1)
                disp('add an intersection first')
            end
            obj.x = [obj.x obj.x(end)+L];
            obj.vi = [obj.vi vin*5280/3600];
            obj.vo = [obj.vo vout*5280/3600];
        end
        
        function []=setWindowType(obj,str)
            obj.windowtype = str;
            obj.initialize();
        end
        
        function []=initialize(obj)
            for i=1:obj.n
                obj.intersection(i).initialize(obj.windowtype);
            end
            obj.goodparameters = true;
        end
        
        function []=optimize(obj)
            
            if(~obj.goodparameters)
                obj.optbandwidth = nan;
                return;
            end
            
            switch obj.windowtype
                case 'pretimed'
                    switch obj.pretimed_optimization_method
                        case {'ip','simplex'}
                            [obj.optbandwidth,obj.optbo,obj.optbi] = obj.optimize_pretimed_lp();
                        case {'milp_ip','milp_mip'}
                            [obj.optbandwidth,obj.optbo,obj.optbi] = obj.optimize_pretimed_milp();
                    end
                    
                case 'gaussian'
                    [obj.optbandwidth,obj.optbo,obj.optbi] = obj.optimize_gaussian();
            end
            
        end
        
        function []=plot(obj,cf)
            
            if(nargin<2)
                cf = gcf;
            else
                figure(cf)
                clf
            end
            
            offsetO = [obj.intersection.absoffseto];
            offsetI = [obj.intersection.absoffseti];
            reloffsetO = [obj.intersection.reloffseto];
            reloffsetI = [obj.intersection.reloffseti];
            
            switch obj.windowtype
                case 'pretimed'
                    greenO = [obj.intersection.go];
                    greenI = [obj.intersection.gi];
                case 'gaussian'
                    sigmaO = [obj.intersection.sigma_o];
                    sigmaI = [obj.intersection.sigma_i];
                    gammaO = [obj.intersection.gamma_o];
                    gammaI = [obj.intersection.gamma_i];
            end
            
            % compute band sizes
            switch obj.windowtype
                case 'pretimed'
                    [band_o,bandoffset_o] = compute_band_pretimed(obj,reloffsetO,offsetO,greenO);
                    [band_i,bandoffset_i] = compute_band_pretimed(obj,reloffsetI,offsetI,greenI);
                    %title(['bi=' num2str(bi) ' , bo=' num2str(bo) ' , b=' num2str(obj.optbandwidth)])
                case 'gaussian'
                    [band_sigma_o,band_mu_o] = computebands_gaussian(obj,reloffsetO,offsetO,sigmaO);
                    [band_sigma_i,band_mu_i] = computebands_gaussian(obj,reloffsetI,offsetI,sigmaI);
                    %title(['bi=' num2str(bi) ' , bo=' num2str(bo) ' , b=' num2str(obj.optbandwidth)])
            end
            
            % draw bands on each segment
            [to,ti]=obj.segment_travel_times();
            
            for i=1:obj.n-1
                
                switch obj.windowtype
                    case 'pretimed'
                        obj.paint_pretimed_inbound_band( band_i,bandoffset_i(i),ti(i),obj.x(i),obj.x(i+1));
                        obj.paint_pretimed_outbound_band(band_o,bandoffset_o(i),to(i),obj.x(i),obj.x(i+1));
                    case 'gaussian'
                        obj.paint_gaussian_outbound_band(band_sigma_o,band_mu_o(i),to(i),obj.x(i),obj.x(i+1));
                        obj.paint_gaussian_inbound_band( band_sigma_i,band_mu_i(i),ti(i),obj.x(i),obj.x(i+1));
                        
                end
            end
            
            
            % draw green intervals .............
            for i=1:obj.n
                h=line([-1 1]*obj.cycle/2,obj.x(i)*[1 1],'Color',[0 0 0],'LineWidth',1);
                movefront(h);
%                 h=text(-obj.cycle/2,obj.x(i),[' ' obj.intersection(i).name],'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',12,'FontWeight','bold','BackgroundColor',[0.85 0.85 0.85]);
%                 movefront(h);
                switch obj.windowtype
                    case 'pretimed'
                        obj.paint_pretimed_window(obj.x(i),offsetO(i),greenO(i),true);
                        obj.paint_pretimed_window(obj.x(i),offsetI(i),greenI(i),false);
                    case 'gaussian'
                        obj.paint_gaussian_window(obj.x(i),offsetO(i),sigmaO(i),gammaO(i),true);
                        obj.paint_gaussian_window(obj.x(i),offsetI(i),sigmaI(i),gammaI(i),false);
                end
                
            end
            
            minY = min(obj.x)-100;
            maxY = max(obj.x)+100;
            set(gca,'YLim',[0.95*minY 1.05*maxY])
            set(gca,'XLim',[-1 1]*obj.cycle/2)
            grid
%            xlabel('time [sec]')
%            ylabel('position [ft]')
            
        end
        
        function [band,bandoffset]=compute_band_pretimed(obj,reloffset,absoffset,green)
            
            % bandoffset is the absolute offset of the center of the band
            interval = reloffset'*[1 1] + green'/2*[-1 1];
            a = max(interval(:,1));
            b = min(interval(:,2));
            if(b>a)
                band = b-a;
                bandoffset = absoffset' - green'/2 + mean([a b]) - interval(:,1);
            else
                band = 0;
                bandoffset = nan*absoffset;
            end
        end
        
        function [Sigma,Mu]=computebands_gaussian(obj,reloffset,absoffset,sigma)
            % Mu_o applies to the first intersection
            [Mu,Sigma] = obj.gaussian_product(reloffset,sigma);
            Mu = absoffset + ( Mu - reloffset) ;
        end
        
        function [b]=modhalf(obj,a)
            b=mod(a+obj.cycle/2,obj.cycle)-obj.cycle/2;
        end
        
        function [total_band,band_o,band_i]=compute_total_bandwidth(obj)
            
            offsetO = [obj.intersection.absoffseto];
            offsetI = [obj.intersection.absoffseti];
            reloffsetO = [obj.intersection.reloffseto];
            reloffsetI = [obj.intersection.reloffseti];
            
            % compute band sizes
            switch obj.windowtype
                case 'pretimed'
                    greenO = [obj.intersection.go];
                    greenI = [obj.intersection.gi];
                    band_o = compute_band_pretimed(obj,reloffsetO,offsetO,greenO);
                    band_i = compute_band_pretimed(obj,reloffsetI,offsetI,greenI);
                    total_band = band_o + band_i;
                    
                case 'gaussian'
                    sigmaO = [obj.intersection.sigma_o];
                    sigmaI = [obj.intersection.sigma_i];
                    
                    reloffsetO = reloffsetO(2:end);
                    reloffsetI = reloffsetI(2:end);
                    
                    gammaO = [obj.intersection.gamma_o];
                    gammaI = [obj.intersection.gamma_i];
                    
                    [band_sigma_o,band_b_o] = obj.integral_gaussian_product(gammaO,sigmaO);
                    band_o = band_b_o*exp(-reloffsetO*band_sigma_o*reloffsetO'/2);
                    
                    [band_sigma_i,band_b_i] = obj.integral_gaussian_product(gammaI,sigmaI);
                    band_i = band_b_i*exp(-reloffsetI*band_sigma_i*reloffsetI'/2);
                    
                    total_band = band_o + band_i;
            end
            
        end
        
    end
    
    methods ( Access = public )
        
        function [totalbandwidth,bo,bi]=optimize_pretimed_lp(obj)
            
            Z=[];   % can be used to hold additional outputs
            
            % obj.initialize(); done in "load"
            
            % compute travel times .......................
            [to,ti]=segment_travel_times(obj);
            
            % compute translated internal offsets .........
            delta0 = translated_internal_offsets(obj,to,ti);
            
            % mean green intervals ....................
            go = [obj.intersection.go];
            gi = [obj.intersection.gi];
            Go = (repmat(go,obj.n,1)+repmat(go',1,obj.n))/2;
            Gi = (repmat(gi,obj.n,1)+repmat(gi',1,obj.n))/2;
            gstaro = min(go);
            gstari = min(gi);
            
            % set up linear program ....................
            % x = [w2 ... wn b bbar]'
            % max ( b + bbar )
            % s.t. b \pm (wi - wj) < Go(i,j) for all pairs i,j in 2...n i~=j
            %      b \pm wi < Go(1,i)     i=2...n
            %      bbar \pm (wi - wj) < Gi(i,j) \mp delta0(i)-delta0(j) for all pairs i,j in 2...n i~=j
            %      bbar \pm wi < Gi(1,i) \mp (delta0(i)-delta0(1))   i=2...n
            %      b < min(go)
            %      bbar < min(gi)
            numvar = obj.n+1;
            A = [];
            B = [];
            b = obj.n;                  % index to "b" in row
            bbar = obj.n+1;             % index to "bbar" in row
            for i=2:obj.n
                
                wi = i-1;               %index to "wi" in row
                
                for j=2:obj.n
                    if(i==j)
                        continue
                    end
                    
                    wj = j-1;          %index to "wj" in row
                    
                    % b + wi - wj < Go(i,j) for all pairs i,j in 2...n i~=j
                    row = zeros(1,numvar);
                    row(b) = 1;
                    row(wi) = 1;
                    row(wj) = -1;
                    A = [A;row];
                    B = [B;Go(i,j)];
                    
                    % bbar + wi - wj < Gi(i,j) - delta0(i) + delta0(j) for all pairs i,j in 2...n i~=j
                    row = zeros(1,numvar);
                    row(bbar) = 1;
                    row(wi) = 1;
                    row(wj) = -1;
                    A = [A;row];
                    B = [B; Gi(i,j) - delta0(i) + delta0(j) ];
                end
                
                % b + wi < Go(1,i)     i=2...n
                row = zeros(1,numvar);
                row(b) = 1;
                row(wi) = 1;
                A = [A;row];
                B = [B; Go(1,i)];
                
                % b - wi < Go(1,i)     i=2...n
                row = zeros(1,numvar);
                row(b) = 1;
                row(wi) = -1;
                A = [A;row];
                B = [B; Go(1,i)];
                
                % bbar + wi < Gi(1,i) - (delta0(i)-delta0(1))   i=2...n
                row = zeros(1,numvar);
                row(bbar) = 1;
                row(wi) = 1;
                A = [A;row];
                B = [B;Gi(1,i)-(delta0(i)-delta0(1))];
                
                % bbar - wi < Gi(i) + (delta0(i)-delta0(1))   i=2...n
                row = zeros(1,numvar);
                row(bbar) = 1;
                row(wi) = -1;
                A = [A;row];
                B = [B;Gi(1,i)+(delta0(i)-delta0(1))];
                
            end
            
            if(strcmp(obj.pretimed_optimization_method,'simplex'))
                [OmegaLstar,fLstar] = linprog([zeros(1,obj.n-1) -1 -1], ... % objective function
                    A,B, ...                               % inequality constraints
                    [],[], ...                             % equality constraints
                    [-inf(obj.n-1,1);0;0],...              % lower bound
                    [ inf(obj.n-1,1);gstaro;gstari],...    % upper bound
                    [], ...
                    optimset('LargeScale','off','Simplex','on') );
            else
                [OmegaLstar,fLstar] = linprog([zeros(1,obj.n-1) -1 -1], ... % objective function
                    A,B, ...                               % inequality constraints
                    [],[], ...                             % equality constraints
                    [-inf(obj.n-1,1);0;0],...              % lower bound
                    [ inf(obj.n-1,1);gstaro;gstari]);      % upper bound
            end
            
            gstar = max([min(go) min(gi)]);
            
            % apply the theorem ..............
            if( -fLstar >= gstar )
                % disp('double band')
                OmegaNstar = OmegaLstar;
            else
                % disp('single band')
                OmegaNstar = nan(numvar,1);
                if( gstaro >= gstari )   % outbound is dominant
                    OmegaNstar(1:obj.n-1) = 0;
                    OmegaNstar(obj.n) = gstar;
                    OmegaNstar(obj.n+1) = 0;
                else                    % inbound is dominant
                    OmegaNstar(1:obj.n-1) = delta0(1)-delta0(2:end);
                    OmegaNstar(obj.n) = 0;
                    OmegaNstar(obj.n+1) = gstar;
                end
                
            end
            
            bo = OmegaNstar(end-1);
            bi = OmegaNstar(end);
            totalbandwidth = bo+bi;
            
            % extract relative offsets ............
            omegaO = [0;OmegaNstar(1:obj.n-1)];
            omegaI = omegaO + delta0' - delta0(1);
            omegaO = obj.modhalf(omegaO);
            omegaI = obj.modhalf(omegaI);
            
            % translate back into absolute offsets ...............
            [thetaO,thetaI] = obj.relative2absolute(omegaO,omegaI,to,ti);
            
            %  copy to intersections ................
            for i=1:obj.n
                obj.intersection(i).set_reloffset(omegaO(i),omegaI(i));
                obj.intersection(i).set_absoffset(thetaO(i),thetaI(i));
            end
            
        end
        
        function [totalbandwidth,bo,bi]=optimize_pretimed_milp(obj)

            addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'milp_solvers')))
            
            % compute travel times .......................
            [to,ti]=segment_travel_times(obj);
            
            to = abs(to/obj.cycle);
            ti = abs(ti/obj.cycle);

            go = [obj.intersection.go]/obj.cycle;
            gi = [obj.intersection.gi]/obj.cycle;
            
            ro = 1 - go;
            ri = 1 - gi;
            
            delta = -[obj.intersection.delta]/obj.cycle;
            
            % set up mixed integer linear program ....................
            % x = [m1...mn w1...wn wbar1...wbarn b bbar]'
            % max ( b + bbar )
            % s.t. 
            % w_i+wbar_i-w_{i+1}-wbar_{i+1}-m_i =
            % -t_i-tbar_i-delta_i+delta_{i+1}-0.5(ri+rbari)+0.5(r_{i+1}+rbar_{i+1})
            % i=1...n-1
            % wi+b leq c-ri    i=1...n
            % wbari+bbar leq c-rbari   i=1...n
            % mi is integer
            % m,w,wbar,b,bbar geq 0    i=1...n
            %

            mindex = 0;
            windex = obj.n;
            wbarindex = 2*obj.n;
            bindex = 3*obj.n+1;
            bbarindex = 3*obj.n + 2;
            
            numvar = 3*obj.n + 2;
            
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            
            for i=1:obj.n
                
                % wi+b leq c-ri    i=1...n
                row = zeros(1,numvar);
                row(windex+i) = 1;
                row(bindex) = 1;
                A = [A;row];
                b = [b; go(i)];
                
                % wbari+bbar leq c-rbari   i=1...n
                row = zeros(1,numvar);
                row(wbarindex+i) = 1;
                row(bbarindex) = 1;
                A = [A;row];
                b = [b; gi(i)];

                % w_i+wbar_i-w_{i+1}-wbar_{i+1}-m_i =
                % -t_i-tbar_i-delta_i+delta_{i+1}-0.5(ri+rbari)+0.5(r_{i+1}+rbar_{i+1})
                if(i<obj.n)
                    row = zeros(1,numvar);
                    row(windex+i) = 1;
                    row(wbarindex+i) = 1;
                    
                    row(windex+i+1) = -1;
                    row(wbarindex+i+1) = -1;
                    row(mindex+i) = -1;
                    
                    rhs = -to(i)-ti(i)-delta(i)+delta(i+1)-0.5*(ro(i)+ri(i))+0.5*(ro(i+1)+ri(i+1));
                    
                    Aeq = [Aeq;row];
                    beq = [beq; rhs];
                    
                end
                
            end
            
            all_mindex = false(1,numvar);
            all_mindex(1:obj.n) = true;
            
            f = zeros(1,numvar);
            f(bindex)=-1;
            f(bbarindex)=-1;
            
            lb = zeros(numvar,1);
            ub = inf(numvar,1);
            
            switch(obj.pretimed_optimization_method)
                case 'milp_ip'
                    X = IP(f,A,b,Aeq,beq,lb,ub,find(all_mindex),2^-24);
                case 'milp_mip'
                    X = miprog(f,A,b,Aeq,beq,lb,ub,all_mindex);
            end
            
            obj.milp_m = X(mindex+(1:obj.n))';
            obj.milp_w = X(windex+(1:obj.n))';
            obj.milp_wbar = X(wbarindex+(1:obj.n))';
            
            bo = X(bindex);
            bi = X(bbarindex);
            totalbandwidth = bo+bi;
            
        end
        
        function [totalbandwidth,bo,bi]=optimize_gaussian(obj)
            
            % compute travel times .......................
            [to,ti]=segment_travel_times(obj);
            
            % compute translated internal offsets .........
            delta0 = translated_internal_offsets(obj,to,ti);
            delta_bold = delta0(1)-delta0(2:end);
            
            % mean green intervals ....................
            sigma_o = [obj.intersection.sigma_o];
            sigma_i = [obj.intersection.sigma_i];
            
            gamma_o = [obj.intersection.gamma_o];
            gamma_i = [obj.intersection.gamma_i];
            
            [Sigma_o,bo_o] = obj.integral_gaussian_product(gamma_o,sigma_o);
            [Sigma_i,bo_i] = obj.integral_gaussian_product(gamma_i,sigma_i);
            
            alpha_o = delta_bold*Sigma_o*delta_bold' / 2;
            alpha_i = delta_bold*Sigma_i*delta_bold' / 2;
            
            e_o = obj.find_max_mixture(bo_o,alpha_o,bo_i,alpha_i,0);
            e_i = obj.find_max_mixture(bo_o,alpha_o,bo_i,alpha_i,1);
            
            Y_eo = bo_o*exp(-alpha_o*e_o.^2) + bo_i*exp(-alpha_i*(1-e_o).^2);
            Y_ei = bo_o*exp(-alpha_o*e_i.^2) + bo_i*exp(-alpha_i*(1-e_i).^2);
            
            if(Y_eo(end)>=Y_ei(end))
                e_opt = e_o(end);
            else
                e_opt = e_i(end);
            end
            
            omegaO_star = e_opt * delta_bold;
            omegaI_star = omegaO_star - delta_bold;
            
            bo = bo_o*exp(-omegaO_star*Sigma_o*omegaO_star'/2);
            bi = bo_i*exp(-omegaI_star*Sigma_i*omegaI_star'/2);
            totalbandwidth = bo+bi;
            
            %plotmixture(bo_o,alpha_o,bo_i,alpha_i,e_o,e_i)
            
            % extract relative offsets ............
            omegaO = [0;omegaO_star'];
            omegaI = omegaO + delta0' - delta0(1);
            omegaO = obj.modhalf(omegaO);
            omegaI = obj.modhalf(omegaI);
            
            % translate back into absolute offsets ...............
            [thetaO,thetaI] = obj.relative2absolute(omegaO,omegaI,to,ti);
            
            %  copy to intersections ................
            for i=1:obj.n
                obj.intersection(i).set_reloffset(omegaO(i),omegaI(i));
                obj.intersection(i).set_absoffset(thetaO(i),thetaI(i));
            end
            
        end
        
        function [to,ti] = segment_travel_times(obj)
            % Equations (1) and (2).
            ti = nan(1,obj.n-1);            % [sec]
            to = nan(1,obj.n-1);            % [sec]
            for i=1:obj.n-1
                ti(i) = (obj.x(i)-obj.x(i+1)) / obj.vi(i);
                to(i) = (obj.x(i+1)-obj.x(i)) / obj.vo(i);
            end
        end
        
        function [delta0] = translated_internal_offsets(obj,to,ti)
            % Equations (33) and (34)
            delta = obj.modhalf([obj.intersection.delta]);      % [sec]
            delta0 = nan(1,obj.n);                              % [sec]
            sumtito = 0;
            for i=1:obj.n-1
                delta0(i) = delta(i) + sumtito;
                sumtito = sumtito + to(i) - ti(i);
            end
            delta0(end) = delta(end) + sumtito;
            %delta0 = obj.modhalf(delta0);
            delta0 = obj.modhalf(delta0-delta0(1))+delta0(1);
                % this is an alternative to straight modhalf which attempts
                % to put all inbound intervals as close to delta0(1) as
                % possible.
            clear sumtito
        end
        
        function [Sigma,bo] = integral_gaussian_product(obj,gamma,sigma)
            
            A = zeros(obj.n-1);
            for i=2:obj.n
                for j=2:obj.n
                    A(i-1,j-1) = (sigma(i)*sigma(j))^-2;
                end
            end
            Sigma = diag(sigma(2:end).^-2) - A/sum(sigma.^-2);
            
            bo = prod(gamma)/sqrt( ((2*pi)^(obj.n-1)) * sum(sigma.^-2) * prod(sigma.^2) );
            
        end
        
        function [Mu,Sigma] = gaussian_product(obj,mu,sigma)
            SigmaSquare = 1/sum(sigma.^-2);               % Eq (73)
            Mu = SigmaSquare * sum( mu .* sigma.^-2 );    % Eq (72)
            Sigma = sqrt(SigmaSquare);
        end
        
        function [e]=find_max_mixture(obj,bo_o,alpha_o,bo_i,alpha_i,start)

            lambda_o = 1/sqrt(2*alpha_o);
            lambda_i = 1/sqrt(2*alpha_i);
            
            k=1;
            e(k) = start;
            maxnumsteps = 10000;
            while(k<maxnumsteps)
                if(e(k)>lambda_o && e(k)<1-lambda_i)
                    if(start==0)
                        e(k+1) = lambda_o - 1e-3;
                    else
                        e(k+1) = 1-lambda_i + 1e-3;
                    end
                else
                    xi_d_o  = -2*alpha_o * bo_o * e(k) * exp(-alpha_o*e(k)^2);                                  % Eq. (64)
                    xi_dd_o =  2*alpha_o * bo_o * ( 2*alpha_o*e(k)^2-1 ) * exp(-2*alpha_o*e(k)^2);              % Eq. (65)
                    xi_d_i  =  2*alpha_i * bo_i * (1-e(k)) * exp(-alpha_i*(1-e(k))^2);                          % Eq. (66)
                    xi_dd_i =  2*alpha_i * bo_i * ( 2*alpha_i*(1-e(k))^2-1 ) * exp(-2*alpha_i*(1-e(k))^2);      % Eq. (67)
                    stepsize = -(xi_d_o+xi_d_i)/(min([0 xi_dd_o]) + min([0 xi_dd_i]));                          % Eq. (68)
                    if(abs(stepsize)>50)
                        disp('large step')
                    end
                    if( abs(stepsize)<1e-10 )
                        break;
                    end
                    e(k+1) = e(k) + stepsize;
                end
                k=k+1;
                
            end
            
        end
        
        function [thetaO,thetaI] = relative2absolute(obj,omegaO,omegaI,to,ti)
            thetaI = nan(1,obj.n);
            thetaO = nan(1,obj.n);
            thetaO(1) = 0;
            thetaI(1) = obj.modhalf(obj.intersection(1).delta);
            for i=2:obj.n
                thetaO(i) = omegaO(i) + thetaO(1) - omegaO(1) + sum(to(1:i-1));
                thetaI(i) = omegaI(i) + thetaI(1) - omegaI(1) + sum(ti(1:i-1));
            end
            thetaO = obj.modhalf(thetaO);
            thetaI = obj.modhalf(thetaI);
        end
        
        function [omegaO,omegaI] = absolute2reltaive(obj,thetaO,thetaI,to,ti)
            omegaI = nan(1,obj.n);
            omegaO = nan(1,obj.n);
            omegaO(1) = 0;
            omegaI(1) = 0;
            for i=2:obj.n
                omegaO(i) = thetaO(i) - thetaO(1) - sum(to(1:i-1));
                omegaI(i) = thetaI(i) - thetaI(1) - sum(ti(1:i-1));
            end
            omegaO = obj.modhalf(omegaO);
            omegaI = obj.modhalf(omegaI);
        end
        
        function []=paint_pretimed_window(obj,x,offset,green,isout)
            minY = min(obj.x)-100;
            maxY = max(obj.x)+100;
            w = (maxY-minY)/50;
            xx = offset-green/2;
            if(isout)
                z = 0.7;
                x = x-w;
            else
                z = 1;
            end
            if(xx<-obj.cycle/2 || xx+green>obj.cycle/2)  % split
                xx = obj.modhalf(xx);
                h=rectangle('Position',[xx x green w],'FaceColor',z*[1 1 1]);
                movefront(h);
                h=rectangle('Position',[xx-obj.cycle x green w],'FaceColor',z*[1 1 1]);
                movefront(h);
            else
                h=rectangle('Position',[xx x green w],'FaceColor',z*[1 1 1]);
                movefront(h);
            end
        end
        
        function []=paint_gaussian_window(obj,x,offset,sigma,gamma,isout)
            xx = linspace(-obj.cycle/2,obj.cycle/2);
            y = obj.eval_gauss(xx,offset,sigma,gamma) + ...
                obj.eval_gauss(xx,offset+obj.cycle,sigma,gamma) + ...
                obj.eval_gauss(xx,offset-obj.cycle,sigma,gamma);
            y = y - min(y);
            y = y/max(y);
            minY = min(obj.x)-100;
            maxY = max(obj.x)+100;
            w = (maxY-minY)/50;
            y = w*y;
            hold on
            if(isout)
                z = 0.7;
                y = x - y;
                h=jbfill(gcf,xx,x+0*y,y,z*[1 1 1],[0 0 0],1,1);
            else
                z = 1;
                y = x + y;
                h=jbfill(gcf,xx,y,x+0*y,z*[1 1 1],[0 0 0],1,1);
            end
            movefront(h);
            
        end
        
        function []=paint_pretimed_inbound_band(obj,band,bandoffset,t,x_i,x_ip)
            if(band>0)
                done = false;
                while( bandoffset + band/2 >= obj.cycle/2 )      % start from leftmost
                    bandoffset=bandoffset-obj.cycle;
                end
                while(~done)
                    p1 = bandoffset - band/2;
                    p2 = bandoffset + band/2;
                    p3 = p1 + t;
                    p4 = p2 + t;
                    patch([p1 p3 p4 p2],[x_i x_ip x_ip x_i],'k','FaceColor',0.5*ones(1,3),'EdgeAlpha',0,'FaceAlpha',0.2);
                    if( min([p3 p4]) < obj.cycle/2 )
                        bandoffset = bandoffset + obj.cycle;        % keep on increasing bandoffset until out of the picture
                    else
                        done = true;
                    end
                end
            end
        end
        
        function []=paint_pretimed_outbound_band(obj,band,bandoffset,t,x_i,x_ip)
            if(band>0)
                done = false;
                while( bandoffset - band/2 <= -obj.cycle/2 )      % start from rightmost
                    bandoffset = bandoffset + obj.cycle;
                end
                while(~done)
                    p1 = bandoffset - band/2;
                    p2 = bandoffset + band/2;
                    p3 = p1 + t;
                    p4 = p2 + t;
                    patch([p1 p3 p4 p2],[x_i x_ip x_ip x_i],'k','FaceColor',0.5*ones(1,3),'EdgeAlpha',0,'FaceAlpha',0.2);
                    if( max([p3 p4]) > -obj.cycle/2)
                        bandoffset = bandoffset - obj.cycle;
                    else
                        done = true;
                    end
                    
                end
            end
        end
        
        function []=paint_gaussian_inbound_band(obj,sigma,bandoffset,t,x_i,x_ip)
            
            band = 3*sigma;
            nn = round(band*10);
            zero = zeros(nn,1);
            done = false;
            
            while( bandoffset + band/2 >= obj.cycle/2 )      % start from leftmost
                bandoffset = bandoffset - obj.cycle;
            end
            while(~done)
                p1 = bandoffset - band/2;
                p2 = bandoffset + band/2;
                p3 = p1 + t;
                p4 = p2 + t;
                
                verts = [ [ zero+x_i  linspace(p1,p2,nn)' ] ; ...
                    [ zero+x_ip linspace(p1,p2,nn)'+t ] ];
                verts(:,[2 1]) = verts(:,[1 2]);
                ind = verts(:,2)==1;
                verts(ind,1)=verts(ind,1)+3;
                faces = [(1:nn-1)' (2:nn)' (nn+2:2*nn)' (nn+1:2*nn-1)'];
                p = patch('Faces',faces,'Vertices',verts,'EdgeColor','none');
                
                cdata = linspace(0,band,nn-1)';
                cdata = 1-exp( -(cdata-band/2).^2/sigma );
                cdata = cdata / max(cdata);
                set(p,'FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','scaled','FaceAlpha',0.8)
                if( min([p3 p4]) < obj.cycle/2)
                    bandoffset = bandoffset + obj.cycle;
                else
                    done = true;
                end
                
            end
            set(gca,'CLim',[0 1])
            colormap('gray')
            
        end
        
        function []=paint_gaussian_outbound_band(obj,sigma,bandoffset,t,x_i,x_ip)
            
            band = 3*sigma;
            nn = round(band*10);
            zero = zeros(nn,1);
            done = false;
            
            while( bandoffset - band/2 <= -obj.cycle/2 )      % start from rightmost
                bandoffset = bandoffset + obj.cycle;
            end
            while(~done)
                p1 = bandoffset - band/2;
                p2 = bandoffset + band/2;
                p3 = p1 + t;
                p4 = p2 + t;
                
                verts = [ [ zero+x_i  linspace(p1,p2,nn)' ] ; ...
                    [ zero+x_ip linspace(p1,p2,nn)'+t ] ];
                verts(:,[2 1]) = verts(:,[1 2]);
                ind = verts(:,2)==1;
                verts(ind,1)=verts(ind,1)+3;
                faces = [(1:nn-1)' (2:nn)' (nn+2:2*nn)' (nn+1:2*nn-1)'];
                p = patch('Faces',faces,'Vertices',verts,'EdgeColor','none');
                
                cdata = linspace(0,band,nn-1)';
                cdata = 1-exp( -(cdata-band/2).^2/sigma );
                cdata = cdata / max(cdata);
                set(p,'FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','scaled','FaceAlpha',0.8)
                if( max([p3 p4]) > -obj.cycle/2)
                    bandoffset = bandoffset - obj.cycle;
                else
                    done = true;
                end
                
            end
            set(gca,'CLim',[0 1])
            colormap('gray')
            
        end
        
        function [y]=eval_gauss(obj,x,mu,sigma,gamma)
            y = gamma * exp( -((x-mu).^2)/2/sigma^2 )/sqrt(2*pi)/sigma;
        end
        
        function [str]=toString(obj)
            str = '';
            
            str = sprintf('%swindowtype\t\t%s\n',str,obj.windowtype);
            str = sprintf('%smethod\t\t\t%s\n',str,obj.pretimed_optimization_method);
            str = sprintf('%scycle\t\t\t%f\n',str,obj.cycle);
            str = sprintf('%sband(o,i,t)\t\t(%f,%f,%f)\n',str,obj.optbo,obj.optbi,obj.optbandwidth);
            str = sprintf('%sSegments\n',str);
            str = sprintf('%s\tx\t\tvo\t\t\t\tvi\n',str);
            for i=1:length(obj.intersection)-1
                str = sprintf('%s\t%d\t\t%f\t\t%f\n',str,i,obj.vo(i),obj.vi(i));
            end
        
            for i=1:length(obj.intersection)
                I = obj.intersection(i);
                str = sprintf('%s\n%s\n',str,I.name);
                
                                
                str = sprintf('%s\tdelta\t\t\t%f\n',str,I.delta);

                if(strcmp(obj.windowtype,'pretimed'))
                    str = sprintf('%s\tgreen(o,i)\t\t(%f,%f)\n',str,I.go,I.gi);
                else
                    str = sprintf('%s\tsigma(o,i)\t\t(%f,%f)\n',str,I.sigma_o,I.sigma_i);
                    str = sprintf('%s\tgamma(o,i)\t\t(%f,%f)\n',str,I.gamma_o,I.gamma_i);
                end
                str = sprintf('%s\tabsoffset(o,i)\t(%f,%f)\n',str,I.absoffseto,I.absoffseti);
                str = sprintf('%s\treloffseto(o,i)\t(%f,%f)\n',str,I.reloffseto,I.reloffseti);
            end
        end
    end
    
end

