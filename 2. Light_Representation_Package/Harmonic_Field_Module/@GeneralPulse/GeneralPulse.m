classdef GeneralPulse
    %GeneralPulse defines general field as set of harmonic fields
    
    properties
        PulseHarmonicFieldSet % Array of harmonic fields with spectral weight included
        % Central wavelength is in the middle = floor(nField/2)
        Direction % default Direction of the pulsed beam
    end
    
    methods
        % constructor
         function NewGeneralPulse = GeneralPulse(harmonicFieldSet,direction)
              if nargin == 0
                  harmonicFieldSet =  HarmonicFieldSet;
                  direction = [0,0,1]';
              elseif nargin == 1
                  direction = [0,0,1]';
              else
                  
              end
              NewGeneralPulse.PulseHarmonicFieldSet = harmonicFieldSet;
              NewGeneralPulse.Direction = direction;              
         end
         
    end
     methods(Static)
        function newObj = InputGUI()
            newObj = ObjectInputDialog(MyHandle(GeneralPulse()));
        end
    end   
end

