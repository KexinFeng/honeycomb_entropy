classdef ColoredObject
    properties
        color
    end
    
    methods
        function obj = ColoredObject(color)
            obj.color = color;
        end
        
        function displayColor(obj)
            disp(['Color: ' obj.color]);
        end
    end
end


