function [multipleRayTracerResult,pupilMeshGrid,outsidePupilIndices ] = ...
        multipleRayTracer(optSystem,wavLenInWavUnit,fieldPointXYInLensUnit,...
        nRay1,nRay2,pupSamplingType,considerPolarization,...
        considerSurfAperture,recordIntermediateResults,endSurface) %
    % Trace bundle of rays through an optical system based on the pupil 
    % sampling specified. Multiple rays can be defined with wavLenInWavUnit (1XnWav),
    % fieldPointXYInLensUnit (2XnField) and the total number of ray will be nRay*nWav*nField
    % That is all rays from each field point with each of wavelegths will be
    % traced. And the result will be 4 dimensional matrix (nNonDummySurface X nRay X nField X nWav). 

    
	% <<<<<<<<<<<<<<<<<<<<<<<<< Author Section >>>>>>>>>>>>>>>>>>>>>>>>>>>>
	%   Written By: Worku, Norman Girma  
	%   Advisor: Prof. Herbert Gross
	%   Part of the RAYTRACE_TOOLBOX V3.0 (OOP Version)
	%	Optical System Design and Simulation Research Group
	%   Institute of Applied Physics
	%   Friedrich-Schiller-University of Jena   
							 
	% <<<<<<<<<<<<<<<<<<< Change History Section >>>>>>>>>>>>>>>>>>>>>>>>>>
	% Date----------Modified By ---------Modification Detail--------Remark
	% Oct 14,2013   Worku, Norman G.     Original Version       Version 3.0
    % Jan 18,2014   Worku, Norman G.     Vectorized version
    
	% <<<<<<<<<<<<<<<<<<<<< Main Code Section >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    % Deafault inputs
    if nargin < 1
        disp('Error: The function multipleRayTracer needs atleast optical system');
        return;
    elseif nargin == 1
        % Take all field points and primary wavelength   
        nField = optSystem.getNumberOfFieldPoints;
        fieldPointMatrix = optSystem.FieldPointMatrix;
        fieldPointXYInLensUnit = (fieldPointMatrix(:,1:2))';
        wavLenInWavUnit = repmat(getPrimaryWavelength(optSystem)/getWavelengthUnitFactor(optSystem),[1,nField]);

        nRay1 = 3;
        nRay2 = 3;
        pupSamplingType = 'Cartesian';
        considerPolarization = 0;
        considerSurfAperture = 1;        
        recordIntermediateResults = 0;        
        endSurface = optSystem.getNumberOfSurfaces;    
    elseif nargin == 2 
        % Take all field points and given wavelength   
        nField = optSystem.getNumberOfFieldPoints;
        fieldPointMatrix = optSystem.FieldPointMatrix;
        fieldPointXYInLensUnit = (fieldPointMatrix(:,1:2))';
        wavLenInWavUnit = repmat(wavLenInWavUnit(1),[1,nField]);
        
        nRay1 = 3;
        nRay2 = 3;
        pupSamplingType ='Cartesian';
        considerPolarization = 0;
        considerSurfAperture = 1;  
        recordIntermediateResults = 0;        
        endSurface = optSystem.getNumberOfSurfaces;
    elseif nargin == 3
        nRay1 = 3;
        nRay2 = 3;
        pupSamplingType ='Cartesian';
        considerPolarization = 0;
        considerSurfAperture = 1;  
        recordIntermediateResults = 0;       
        endSurface = optSystem.getNumberOfSurfaces;
     elseif nargin == 4
        nRay2 = 3;
        pupSamplingType ='Cartesian';
        considerPolarization = 0;
        considerSurfAperture = 1;
        recordIntermediateResults = 0;        
        endSurface = optSystem.getNumberOfSurfaces;
    elseif nargin == 5
        pupSamplingType ='Cartesian';
        considerPolarization = 0;
        considerSurfAperture = 1;  
        recordIntermediateResults = 0;        
        endSurface = optSystem.getNumberOfSurfaces;
    elseif nargin == 6
        considerPolarization = 0;
        considerSurfAperture = 1; 
        recordIntermediateResults = 0;        
        endSurface = optSystem.getNumberOfSurfaces;
    elseif nargin == 7
        considerSurfAperture = 1; 
        recordIntermediateResults = 0;        
        endSurface = optSystem.getNumberOfSurfaces;
    elseif nargin == 8
        recordIntermediateResults = 0;       
        endSurface = optSystem.getNumberOfSurfaces;
    elseif nargin == 9        
        endSurface = optSystem.getNumberOfSurfaces;  
    elseif nargin == 11
     
    end
    
    tic
    
    % Determine the number of non dummy surfaces used for ray tracing 
    % That can be used for final reshaping of the ray trace result matrix.
    NonDummySurfaceIndices = optSystem.getNonDummySurfaceIndices;
    startNonDummyIndex = 1; % Ray trace start from object surface
    indicesBeforeEndSurf = find(NonDummySurfaceIndices<=endSurface);
    endNonDummyIndex = indicesBeforeEndSurf(end);
    
    nNonDummySurfaceConsidered = endNonDummyIndex - startNonDummyIndex + 1;
    
    %% Compute initial ray bundle parameters
    pupilShape = 'Circular';
    nField = size(fieldPointXYInLensUnit,2);
    nWav  = size(wavLenInWavUnit,2);
    
    [ initialRayBundle, pupilSamplingPoints,pupilMeshGrid,...
     outsidePupilIndices   ] = ...
    computeInitialRayArray( optSystem, wavLenInWavUnit,...
    fieldPointXYInLensUnit, nRay1,nRay2,pupSamplingType);
    nRayTotal = size(pupilSamplingPoints,2);
  
    %===============RAYTRACE For Bundle of Ray========================
    rayTraceResult = rayTracer(optSystem,initialRayBundle,considerPolarization,considerSurfAperture,recordIntermediateResults,endSurface); 
    if recordIntermediateResults
        multipleRayTracerResult = reshape(rayTraceResult,[nNonDummySurfaceConsidered,nRayTotal,nField,nWav]); %(nNonDummySurface X nRay X nField X nWav)
    else
        multipleRayTracerResult = reshape(rayTraceResult,[2,nRayTotal,nField,nWav]); %(2 X nRay X nField X nWav)
    end
    pupilCoordinates = pupilSamplingPoints;
    timeElapsed =  toc;
    disp(['Ray Bundle Trace Completed. Polarized  = ',num2str(considerPolarization), ...
        ', Total Number  = ', num2str(nRayTotal*nField*nWav), ', Time Elapsed = ', ...
        num2str(timeElapsed)]); 
end  