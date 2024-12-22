classdef (Abstract) IShape
    methods (Abstract)
        area = calculateArea(obj)
        perimeter = calculatePerimeter(obj)
    end
end
