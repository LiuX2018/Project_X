function [fieldNames,fieldFormat,uniqueParamStruct] = getApertureUniqueParameters( variableInputArgument )
    %getApertureUniqueParameters Returns the field names, formats, and current
    %struct of all unique parameters which are specific to this Aperture type
    % Inputs:
    % variableInputArgument:
    %       1. Struct/object --> currentAperture
    %       2. Char --> ApertureType
    % Outputs:
    %   [fieldNames,fieldFormat,uniqueParamStruct]
    
    % <<<<<<<<<<<<<<<<<<<<<<<<< Author Section >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %   Written By: Worku, Norman Girma
    %   Advisor: Prof. Herbert Gross
    %	Optical System Design and Simulation Research Group
    %   Institute of Applied Physics
    %   Friedrich-Schiller-University of Jena
    
    % <<<<<<<<<<<<<<<<<<< Change History Section >>>>>>>>>>>>>>>>>>>>>>>>>>
    % Date----------Modified By ---------Modification Detail--------Remark
    % Jun 17,2015   Worku, Norman G.     Original Version
    
    if nargin == 0
        returnDefault = 1;
        apertureType = 'FloatingCircularAperture';
    elseif isAperture(variableInputArgument)
        currentAperture = variableInputArgument;
        returnDefault = 0;
        apertureType = currentAperture.Type;
    elseif ischar(variableInputArgument)
        returnDefault = 1;
        apertureType = variableInputArgument;
    else
        disp('Error: Invalid input to getApertureUniqueParameters. So it is just ignored.');
        returnDefault = 1;
        apertureType = 'FloatingCircularAperture';
    end
    
    % Connect the surface definition function
    apertureDefinitionHandle = str2func(apertureType);
    returnFlag = 1;
    [fieldNames,fieldFormat,defaultUniqueParamStruct] = apertureDefinitionHandle(returnFlag);
    if returnDefault
        uniqueParamStruct = defaultUniqueParamStruct;
    else
        uniqueParamStruct = currentAperture.UniqueParameters;
    end  
end

