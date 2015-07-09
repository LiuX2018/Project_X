function [ gridType ] = getGridType( currentSurface )
    %GETGRIDTYPE Returns the Grid Type of the surface
    % Inputs:
    %   (currentSurface)
    % Outputs:
    %   [gridType]
    
    % <<<<<<<<<<<<<<<<<<<<<<<<< Author Section >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %   Written By: Worku, Norman Girma
    %   Advisor: Prof. Herbert Gross
    %	Optical System Design and Simulation Research Group
    %   Institute of Applied Physics
    %   Friedrich-Schiller-University of Jena
    
    % <<<<<<<<<<<<<<<<<<< Change History Section >>>>>>>>>>>>>>>>>>>>>>>>>>
    % Date----------Modified By ---------Modification Detail--------Remark
    % Jun 17,2015   Worku, Norman G.     Original Version
    
    % It depends on the outer aperture shape
    outerApertShape = currentSurface.Aperture.OuterShape;
    switch lower(outerApertShape)
        case {'circular','elliptical'}
            gridType = 'Polar';
        case 'rectangular'
            gridType = 'Rectangular';
        otherwise
            disp('Error: Unsapported aperture outer shape. So the default rectangular grid is used.');
            gridType = 'Rectangular';
    end
    
%     surfaceType = currentSurface.Type;
%     % Connect the surface definition function
%     surfaceDefinitionHandle = str2func(surfaceType);
%     returnFlag = 'GRTY';
%     gridType = surfaceDefinitionHandle(returnFlag);
end

