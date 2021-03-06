function FNO = getObjectFNo(OS)
    % getObjectFN: returns object space FN of the optical system
    switch OS.SystemApertureType
        case 1 % given 'Enterance Pupil Diameter'

        case 2 % given 'Object Space NA'

        case 3 % given 'Object Space F#'
           FNO = OS.SystemApertureValue;
        case 4 % given 'Image Space NA'

        case 5 % given 'Image Space F#'                 

        otherwise

    end            
end