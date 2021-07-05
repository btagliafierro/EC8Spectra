classdef designSpectrum < handle   

    
    properties(SetAccess = private)        

        S
        Tb
        Tc
        Td
        
    end
    properties
        Sda  
        time
        pseudoAcc       
        lambda
        period
        multiplier=1.00
        ag = 0.25*9.81   %%% ag * g * IF
        importanceFactor=1.00   
        behaviourFactor=1.0
        type    = 'Type 1'
        soil    = 'B'
    end
    
    methods
        function obj = designSpectrum(soil,type)
            % Calculate target spectrum EC8
                           
            if nargin<1
            else
                obj.soil=soil;
                obj.type=type;
            end
            soil=obj.soil;
            if strcmpi(obj.type,'Type 1')
                if strcmp(soil,'A')
                    S=1.00; Tb=0.15; Tc=0.4; Td=2.;
                elseif strcmp(soil,'B')
                    S=1.2; Tb=0.15; Tc=0.5; Td=2.;
                elseif strcmp(soil,'C')
                    S=1.15; Tb=0.20; Tc=0.6; Td=2.;
                elseif strcmp(soil,'D')
                    S=1.35; Tb=0.20; Tc=0.8; Td=2.;
                elseif strcmp(soil,'E')
                    S=1.4; Tb=0.15; Tc=0.5; Td=2.;
                end
            elseif strcmpi(obj.type,'Type 2')
                if strcmp(soil,'A')
                    S=1.00; Tb=0.05; Tc=0.25; Td=1.2;
                elseif strcmp(soil,'B')
                    S=1.35; Tb=0.05; Tc=0.25; Td=1.2;
                elseif strcmp(soil,'C')
                    S=1.50; Tb=0.10; Tc=0.25; Td=1.2;
                elseif strcmp(soil,'D')
                    S=1.80; Tb=0.10; Tc=0.30; Td=1.2;
                elseif strcmp(soil,'E')
                    S=1.60; Tb=0.05; Tc=0.25; Td=1.2;
                end
            end
            obj.S=S;
            obj.Tb=Tb;
            obj.Tc=Tc;
            obj.Td=Td;
        end
        
        function set.type(obj,data)
            obj.type=data;
        end
        function set.soil(obj,data)
            obj.soil=data;
        end
        function set.period(obj,data)
            obj.period=data;
        end
        function set.multiplier(obj,data)
            obj.multiplier=data;
        end  
        function set.ag(obj,data)
            obj.ag=data;
        end
        function set.importanceFactor(obj,data)
            obj.importanceFactor=data;
        end 
        
        function m=get.Sda(obj) 
            index=find(obj.time>=obj.period,1);
            m=obj.pseudoAcc(index);            
        end
        
        function m=get.pseudoAcc(obj) 
            dt=0.001;
            beta=0.2;            
            acceleration=obj.ag*obj.importanceFactor*obj.multiplier;
            q=obj.behaviourFactor;
            T1=@(T) acceleration*obj.S*(2/3+T/obj.Tb*(2.5/q-2/3));
            T2=@(T) acceleration*obj.S*2.5/q+T*0;
            T3=@(T) max(acceleration*obj.S*2.5*(obj.Tc./T)/q , beta*acceleration);
            T4=@(T) max(acceleration*obj.S*2.5*(obj.Tc*obj.Td./T.^2)/q , beta*acceleration);
            n=10;n1=1;
            EC8=[T1(dt/n1:dt/n1:obj.Tb) T2(obj.Tb+dt/n1:dt/n1:obj.Tc) T3(obj.Tc+dt*n/2:dt*n/2:obj.Td) T4(obj.Td+dt*n:dt*n:10)];
            m=obj.multiplier*EC8;
        end
        function m=get.time(obj)   
            dt=0.001;
            n=10;n1=1;
            T=[(dt:dt/n1:obj.Tb) (obj.Tb+dt/n1:dt/n1:obj.Tc) (obj.Tc+dt*n/2:dt*n/2:obj.Td) (obj.Td+dt*n:dt*n:10)];
            m=T;
        end

        
    end
end

