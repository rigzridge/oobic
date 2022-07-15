classdef LoadingBar < handle
    
    properties 
        
        Figure
        Axes
        Index
        Limit
        Color
        
    end % properties
    
    methods
        
        function lb = LoadingBar(lim,cmap)
        % Constructor    
            lb.Limit = lim;
            lb.Color  = eval(sprintf('%s(%d)',cmap,lim));
            lb.Figure = figure('Position',[500 400 60 40],...
                                'Resize','off',...
                                'MenuBar','none',...
                                'WindowStyle','modal');
            lb.Axes = axes('XLimMode','manual','YLimMode','manual',...
                            'DrawMode','fast',...
                            'Parent',lb.Figure);
        end % LoadingBar
        
        function lb = TickBar(lb,index)
        % Change bar
            barh(lb.Figure,index,'FaceColor',lb.Color(index,:));
            axis(lb.Axes,[0 lb.Limit 0.75 1]);
            axis(lb.Axes,'off');
            pause(eps);
        end % TickBar
        
        
        function lb = KillBar(lb)
        % End loading bar
            delete(lb.Figure);
            %eval(sprintf('clear(''%s'')','lb'))
        end % KillBar
        
        
        function lb = Demo(lb)
        %
            N = lb.Limit;
            %u = lb.LoadingBar(N,'hsv');
            for k=1:N
                lb = lb.TickBar(k);
                %pause(0.05);
                pause(eps)
            end
            lb.KillBar;
        end
        
    end
    
end

