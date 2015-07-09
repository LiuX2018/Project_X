function [ dispName, imgFileName, surfDescription ] = getSurfaceAboutInfo( currentComponent )
    %getSurfaceAboutInfo: Returns the information related to the surface
    %type of current surface
    % Inputs:
    %   (currentSurface)
    % Outputs:
    %   [dispName, imgFileName, surfDescription]
    
    % <<<<<<<<<<<<<<<<<<<<<<<<< Author Section >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %   Written By: Worku, Norman Girma
    %   Advisor: Prof. Herbert Gross
    %	Optical System Design and Simulation Research Group
    %   Institute of Applied Physics
    %   Friedrich-Schiller-University of Jena
    
    % <<<<<<<<<<<<<<<<<<< Change History Section >>>>>>>>>>>>>>>>>>>>>>>>>>
    % Date----------Modified By ---------Modification Detail--------Remark
    % Jun 17,2015   Worku, Norman G.     Original Version
    
    surfaceDefinitionFileName = currentComponent.Type;
    % Connect the surface definition function
    surfaceDefinitionHandle = str2func(surfaceDefinitionFileName);
    returnFlag = 'ABTS'; % About the surface
    [ dispName, imgFileName, surfDescription] = surfaceDefinitionHandle(returnFlag);
end
