classdef elasticSpectrum < handle
        
    properties(SetAccess = private)     
        S
        Tb
        Tc
        Td
    end
    properties
        type
        soil    
        ag      
        damping
        Sda  
        time
        pseudoAcc       
        lambda
        period
    end
    
    methods
        function obj = elasticSpectrum(type,soil)
            % Calculate target spectrum EC8
            if nargin<1
            else
                obj.soil=soil;
                obj.type=type;
            end
            soil=obj.soil;
            if strcmpi(obj.type,'Type_1')
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
            elseif strcmpi(obj.type,'Type_2')
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
        
        function set.damping(obj,temp)             
            obj.damping=temp;            
        end
        function set.type(obj,temp)             
            obj.type=temp;            
        end
        function set.soil(obj,temp)             
            obj.soil=temp;            
        end
        function set.period(obj,data)
            obj.period=data;
        end        
        function set.ag(obj,data)
            obj.ag=data;
        end      
        function m=get.Sda(obj) 
            index=find(obj.time>=obj.period,1);
            m=obj.pseudoAcc(index);            
        end
        
        function m=get.pseudoAcc(obj)            
            dt=0.001;   
            eta=(10/(5+obj.damping*100))^0.5;
            T1=@(T) obj.ag*obj.S*(1+T/obj.Tb*(2.5*eta-1));
            T2=@(T) obj.ag*obj.S*eta*2.5*T./T;
            T3=@(T) obj.ag*obj.S.*eta*2.5*(obj.Tc./T);
            T4=@(T) obj.ag*obj.S*eta*2.5*(obj.Tc*obj.Td./T.^2);
            n=10;n1=1;
            EC8=[T1(dt/n1:dt/n1:obj.Tb) T2(obj.Tb+dt/n1:dt/n1:obj.Tc) T3(obj.Tc+dt*n/2:dt*n/2:obj.Td) T4(obj.Td+dt*n:dt*n:10)];

            m=EC8;
        end
        
        function m=get.time(obj)   
            dt=0.001;
            n=10;n1=1;
            T=[(dt/n1:dt/n1:obj.Tb) (obj.Tb+dt/n1:dt/n1:obj.Tc) (obj.Tc+dt*n/2:dt*n/2:obj.Td) (obj.Td+dt*n:dt*n:10)];
            m=T;
        end
        
    end
end

