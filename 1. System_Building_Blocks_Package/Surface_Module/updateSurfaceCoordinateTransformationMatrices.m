function updatedSurfaceArray = updateSurfaceCoordinateTransformationMatrices(surfaceArray)
    %updateSurfaceCoordinateTransformationMatrices Computes the surface
    %coordinate transformation matrices of each surface from the given
    %parameters of the surface thickness, tilt and decenter parameters.
    % If the first surface is an object surface then the reference will be
    % the second surface but otherwise the first surface is used as
    % global reference.
    
    if surfaceArray(1).ObjectSurface
        % Object surface
        objThickness = abs(surfaceArray(1).Thickness);
        if objThickness > 10^10 % Replace Inf with INF_OBJ_Z = 0 for graphing
            objThickness = 0;
        end
        % since global coord but shifted by objThickness
        refCoordinateTM =   [1,0,0,0;
            0,1,0,0;
            0,0,1,-objThickness;
            0,0,0,1];
    else
        refCoordinateTM =   [1,0,0,0;
            0,1,0,0;
            0,0,1,0;
            0,0,0,1];
    end
    
    %     updatedSystem = savedOpticalSystem;
    %     if ~IsSurfaceBased(savedOpticalSystem)
    %         return;
    %     end
    
    updatedSurfaceArray = surfaceArray;
    nSurface = length(surfaceArray);
    
    
    surfaceCoordinateTM = refCoordinateTM;
    referenceCoordinateTM = refCoordinateTM;
    % set surface property
    updatedSurfaceArray(1).SurfaceCoordinateTM = surfaceCoordinateTM;
    updatedSurfaceArray(1).ReferenceCoordinateTM = referenceCoordinateTM;
    for ss = 2:nSurface
        currentSurface = updatedSurfaceArray(ss);
        previousSurface = updatedSurfaceArray(ss-1);
        % Update the surface coordinates and positions
        % compute position from decenter and thickness
        
        prevRefCoordinateTM = previousSurface.ReferenceCoordinateTM;
        prevSurfCoordinateTM = previousSurface.SurfaceCoordinateTM;
        prevThickness = previousSurface.Thickness;
        if prevThickness > 10^10 % Replace Inf with INF_OBJ_Z = 0 for object distance
            prevThickness = 0;
        end
        [surfaceCoordinateTM,referenceCoordinateTM] = currentSurface.TiltAndDecenter(prevRefCoordinateTM,...
            prevSurfCoordinateTM,prevThickness);
        % set surface property
        currentSurface.SurfaceCoordinateTM = surfaceCoordinateTM;
        currentSurface.ReferenceCoordinateTM = referenceCoordinateTM;
        updatedSurfaceArray(ss) = currentSurface;
    end
    
end

