function [ incidenceAngle ] = getAllSurfaceIncidenceAngle( allSurfaceRayTraceResult,...
        rayPupilIndices,rayFieldIndices,rayWavelengthIndices)
    %getAllSurfaceIncidenceAngle: Returns the exit ray direction of a specific
    % ray specified by (rayPupilIndex,rayFieldIndex,rayWavIndex) for all surfaces
    % Input:
    %   allSurfaceRayTraceResult: vector of raytrace result (size = nSurf)
    %   (rayPupilIndex,rayFieldIndex,rayWavIndex) : Indices specifying a given
    %   ray
    % Output:
    %   incidenceAngle: is (1 X nSurface X nPupilPointsRequested X nFieldRequested X nWavRequested)
    
    if nargin == 0
        disp(['Error: The function  getAllSurfaceIncidenceAngle requires ',...
            'atleast the surface trace result struct as argument.']);
        incidenceAngle = NaN;
        return;
    elseif nargin == 1
        rayPupilIndices = 0; % All
        rayFieldIndices = 0; % All
        rayWavelengthIndices = 0; % All
    elseif nargin == 2
        rayFieldIndices = 0; % All
        rayWavelengthIndices = 0; % All
    elseif nargin == 3
        rayWavelengthIndices = 0; % All
    else
        
    end
    
    requestedResultFieldName = 'IncidenceAngle';
    requestedFieldFirstDim = 1;
    incidenceAngle = getRayTraceResultFieldForAllSurfaces( ...
        allSurfaceRayTraceResult,requestedResultFieldName,requestedFieldFirstDim,...
        rayPupilIndices,rayFieldIndices,rayWavelengthIndices);
end

