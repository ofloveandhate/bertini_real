classdef scene_manipulator
	
	properties
      sceneprop1;
   end
   methods
       function obj = scene_manipulator(new_prop_val)
           if  nargin > 0
              obj.sceneprop1 = new_prop_val;
           end
       end
   end
end