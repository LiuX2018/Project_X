function [ampTs,ampTp,powTs,powTp,jonesTransMatrix] = getTransmissionCoefficients...
        (coating,wavLenInUm,incAngleInDeg,indexBefore,indexAfter,referenceWavLenInUm)
    % getTransmissionCoefficients: Computes the amplitude and power coefficients of
    % transmission using general Fresnel's equations. The function is vectorized so
    % it can work on multiple sets of inputs once at the same time. i.e incAngle 
    % or wavLen becomes array)
    % Inputs:
    %   (coating,wavLenInUm,incAngleInDeg,substrateGlass,claddingGlass,referenceWavLenInUm)
    % Outputs:
    %   [ampTs,ampTp,powTs,powTp,jonesTransMatrix]
    
    % <<<<<<<<<<<<<<<<<<<<<<<<< Author Section >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %   Written By: Worku, Norman Girma
    %   Advisor: Prof. Herbert Gross
    %	Optical System Design and Simulation Research Group
    %   Institute of Applied Physics
    %   Friedrich-Schiller-University of Jena
    
    % <<<<<<<<<<<<<<<<<<< Change History Section >>>>>>>>>>>>>>>>>>>>>>>>>>
    % Date----------Modified By ---------Modification Detail--------Remark
    % Oct 14,2013   Worku, Norman G.     Original Version       Version 3.0
    % Jan 18,2014   Worku, Norman G.     Vectorized inputs and outputs
    % Jun 17,2015   Worku, Norman G.     Support the user defined coating
    %                                    definitions
    
    % <<<<<<<<<<<<<<<<<<<<< Main Code Section >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nRayAngle = size(incAngleInDeg,2);
    nRayWav = size(wavLenInUm,2);
    if nRayAngle == 1
        nRay = nRayWav;
        incAngleInDeg = repmat(incAngleInDeg,[1,nRay]);
    elseif nRayWav == 1
        nRay = nRayAngle;
        wavLenInUm = repmat(wavLenInUm,[1,nRay]);
    elseif nRayAngle == nRayWav % Both wavLen and incAngle for all rays given
        nRay = nRayAngle;
    else
        disp(['Error: The size of Incident Angle and Wavelength should '...
            'be equal or one of them should be 1.']);
        return;
    end
    
    coatingType = coating.Type;
    coatingParameters = coating.Parameters;
    
    % Connect to Coating Defintion function
    if strcmpi(coatingType,'None')
        coatingDefinitionHandle = @NullCoating;
    else
        coatingDefinitionHandle = str2func(coatingType);
    end
    returnFlag = 2; % Jones matrix
    wavLen = wavLenInUm*10^-6;
    referenceWavLen = referenceWavLenInUm*10^-6;
    [ ampTransJonesMatrix,ampRefJonesMatrix,powTransJonesMatrix,powRefJonesMatrix] = ...
        coatingDefinitionHandle(returnFlag,coatingParameters,wavLen,referenceWavLen,...
        incAngleInDeg,indexBefore,indexAfter);
    
    ampTs = ampTransJonesMatrix(1,1,:);
    ampTp = ampTransJonesMatrix(2,2,:);
    powTs = powTransJonesMatrix(1,1,:);
    powTp = powTransJonesMatrix(2,2,:);
    jonesTransMatrix = ampTransJonesMatrix;
end
