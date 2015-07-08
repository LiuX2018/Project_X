function [ isInsideTheAdditionalEdge ] = IsInsideTheAdditionalEdge( surfAperture, xyVector )
    % IsInsideTheOuterAperture Returns 1 if the xyVector is insied the outer aperture and
    % 0 otherwise.
    % Inputs:
    %   ( surfAperture, xyVector )
    % Outputs:
    %    [isInsideTheAdditionalEdge]
    
    % <<<<<<<<<<<<<<<<<<<<<<<<< Author Section >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %   Written By: Worku, Norman Girma
    %   Advisor: Prof. Herbert Gross
    %	Optical System Design and Simulation Research Group
    %   Institute of Applied Physics
    %   Friedrich-Schiller-University of Jena
    
    % <<<<<<<<<<<<<<<<<<< Change History Section >>>>>>>>>>>>>>>>>>>>>>>>>>
    % Date----------Modified By ---------Modification Detail--------Remark
    % Jun 19,2015   Worku, Norman G.     Original Version
    
    apertureType = surfAperture.Type;
    apertDecX = surfAperture.Decenter(1);
    apertDecY = surfAperture.Decenter(2);
    apertRotInDeg = surfAperture.Rotation;
    
    % First decenter and then rotate the given points
    xyVector_Decenter = [xyVector(:,1) - apertDecX , xyVector(:,2) - apertDecY];
    
    initial_r = computeNormOfMatrix( xyVector_Decenter, 2 );
    initial_angleInRad = atan(xyVector_Decenter(:,2)./xyVector_Decenter(:,1));
    new_angleInRad = initial_angleInRad - apertRotInDeg*pi/180;
    
    xyVector_final = [initial_r.*cos(new_angleInRad), initial_r.*sin(new_angleInRad)];
    

    my_eps = 0;
    
%     xyVector_final = xyVector;
    % Now connect to the aperture defintion function and compute the
    % maximum Radius in x and y
    apertureDefinitionHandle = str2func(apertureType);
    returnFlag = 2; % maximumRadiusXY
    apertureParameters = surfAperture.UniqueParameters;
    [ maximumRadiusXY] = ...
        apertureDefinitionHandle(returnFlag,apertureParameters,xyVector_final);
    maximumRadiusXYWithEdge = maximumRadiusXY*(1+surfAperture.AdditionalEdge);
    switch lower(surfAperture.OuterShape)
        case {'elliptical','circular'}
            semiDiamX1 = maximumRadiusXY(1);
            semiDiamY1 = maximumRadiusXY(2);
            
            semiDiamX2 = maximumRadiusXYWithEdge(1);
            semiDiamY2 = maximumRadiusXYWithEdge(2);
            
            pointX = xyVector(:,1);
            pointY = xyVector(:,2);
            umInsideTheInnerEllipse = ((((pointX).^2)/semiDiamX1^2) + (((pointY).^2)/semiDiamY1^2) < 1 + my_eps);
            umInsideTheOuterEllipse  = ((((pointX).^2)/semiDiamX2^2) + (((pointY).^2)/semiDiamY2^2) < 1 + my_eps);
            isInsideTheAdditionalEdge = umInsideTheOuterEllipse & ~umInsideTheInnerEllipse ;
            
        case 'rectangular'
            semiDiamX1 = maximumRadiusXY(1);
            semiDiamY1 = maximumRadiusXY(2);
            
            semiDiamX2 = maximumRadiusXYWithEdge(1);
            semiDiamY2 = maximumRadiusXYWithEdge(2);
            
            pointX = xyVector(:,1);
            pointY = xyVector(:,2);
            umInsideTheInnerRectangle = (abs(pointX) < semiDiamX1 + my_eps & abs(pointY) < semiDiamY1 + my_eps);
            umInsideTheOuterRectangle = (abs(pointX) < semiDiamX2 + my_eps & abs(pointY) < semiDiamY2 + my_eps);
            isInsideTheAdditionalEdge = umInsideTheOuterRectangle & ~umInsideTheInnerRectangle;
        otherwise
    end
    
end

