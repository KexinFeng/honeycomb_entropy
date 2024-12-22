classdef Rectangle < ColoredObject & NamedObject & IShape
    properties
        width
        height
    end
    
    methods
        function obj = Rectangle(name, color, width, height)
            obj@ColoredObject(color);
            obj@NamedObject(name);
            obj.width = width;
            obj.height = height;
        end
        
        function area = calculateArea(obj)
            area = obj.width * obj.height;
        end
        
        function perimeter = calculatePerimeter(obj)
            perimeter = 2 * (obj.width + obj.height);
        end
    end
end
