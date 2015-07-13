function rayPathMatrix = computeRayPathMatrix...
        (optSystem,wavLen,fieldPointXY,PupSamplingType,nRay1,nRay2)
    % computeRayPathMatrix: computes the ray path coordinates for all field
    % points and wavelengths by performing multiple ray trace. The
    % wavInd,fldInd can be vectors. Then the output rayPathMatrix will be
    % of 5D dimension. (3 X nSurf X nPupilPoints X nFieldPoints X nWav)
    
    considerPolarization = 0;
    considerSurfAperture = 1;
    recordIntermediateResults = 1;
    % The polarizedRayTraceResult will be a a vector of RayTraceResult
    % object of size = nSurf.
    [polarizedRayTracerResult] =  ...
        multipleRayTracer(optSystem,wavLen,fieldPointXY,nRay1,nRay2,PupSamplingType,considerPolarization,...
        considerSurfAperture,recordIntermediateResults);
    
    nSurface = length(polarizedRayTracerResult);
    [ exitRayPositions ] = getAllSurfaceExitRayPosition( polarizedRayTracerResult);
    [ rayIntersectionPoints ] = getAllSurfaceRayIntersectionPoint( polarizedRayTracerResult);
    
    rayPathMatrix(:,[1:2:2*nSurface],:,:,:) = rayIntersectionPoints;
    rayPathMatrix(:,[2:2:2*nSurface],:,:,:)  = exitRayPositions;
end
