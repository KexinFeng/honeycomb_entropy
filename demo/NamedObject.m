classdef NamedObject
    properties
        name
    end
    
    methods
        function obj = NamedObject(name)
            obj.name = name;
        end
        
        function displayName(obj)
            disp(['Name: ' obj.name]);
        end
    end
end