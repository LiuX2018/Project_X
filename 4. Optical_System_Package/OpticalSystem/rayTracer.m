function rayTracerResultReshaped = rayTracer(optSystem, objectRayMatrix,rayTraceOptionStruct,endSurface,nRayPupil,nField,nWav)
    % rayTracer: main function of polarized ray tracer from object surface to the
    % endSurface (inclusive)
    % The function is vectorized so it can work on multiple sets of
    % inputs once at the same time. That is array of ray objects.
    % Inputs
    %   optSystem: data type "OpticalSystem"
    %   ObjectRay: data type "Ray" or array of Ray object
    %   startSurf,endSurf: Indices of start and end surface. ObjectRay is
    %   assumed to be given just after the start surface and it will be traced
    %   till end surface (inclusive)
    %   rayTraceOptionStruct: struct with options indicating what to
    %   compute and consider during ray trace. (See RayTraceOptionStruct()
    %   function for more details)
    % Output:
    %   rayTracerResult: (array if all surface results are recorded) of "RayTraceResult" or can be
    %   matrix of RayTraceResult objects if the input is array of Ray
    %   object. Size : nSurface X nTotalRay
    
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
    % Oct 22,2014   Worku, Norman G.     Avoid recording of intermediate surface
    %                                    results if not neccessary or stated
    % <<<<<<<<<<<<<<<<<<<<< Main Code Section >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
%     profile on
    % Deafualt arguments
    if nargin < 2
        disp(['Error: Missing Input. The function rayTracer needs '...
            ' atleast the optical system and an object ray as input argument.']);
        rayTracerResult = NaN;
        return;
    elseif nargin == 2
        rayTraceOptionStruct = RayTraceOptionStruct();
        endSurface = getNumberOfSurfaces(optSystem);
        nRayPupil = 1;
        nField = 1;
        nWav = 1;
    elseif nargin == 3
        endSurface = getNumberOfSurfaces(optSystem);
        nRayPupil = 1;
        nField = 1;
        nWav = 1;
    elseif nargin == 4
        nRayPupil = 1;
        nField = 1;
        nWav = 1;
    elseif nargin == 5
        nField = 1;
        nWav = 1;
    elseif nargin == 6
        nWav = 1;
    else
        
    end
    
    considerPolarization = rayTraceOptionStruct.ConsiderPolarization;
    recordIntermediateResults = rayTraceOptionStruct.RecordIntermediateResults;
    considerSurfAperture = rayTraceOptionStruct.ConsiderSurfAperture;
    
    computeGeometricalPathLength = rayTraceOptionStruct.ComputeGeometricalPathLength;
    computeOpticalPathLength = rayTraceOptionStruct.ComputeOpticalPathLength;
    computeGroupPathLength = rayTraceOptionStruct.ComputeGroupPathLength;
    
    computeRefractiveIndex = rayTraceOptionStruct.ComputeRefractiveIndex;
    computeRefractiveIndexFirstDerivative = rayTraceOptionStruct.ComputeRefractiveIndexFirstDerivative ;
    computeRefractiveIndexSecondDerivative  = rayTraceOptionStruct.ComputeRefractiveIndexSecondDerivative ;
    computeGroupIndex = rayTraceOptionStruct.ComputeGroupIndex;
    
    
    % Determine the size of ObjectRay matrix so as to return rayTraceResult
    % array of the same matrix
    sizeOfInputObjectRay = size(objectRayMatrix);
    
    objectRayBundle = objectRayMatrix(:);
    % Determine the start and end indeces in non dummy surface array
    %     [surfaceArray,nSurface, nonDummySurfaceArray, nonDummySurfaceArray,...
    %         nonDummySurfaceIndices, nNonDummySurface] = getSurfaceArrayAllInformation(optSystem);
    
    [ NonDummySurfaceArray,nNonDummySurface,NonDummySurfaceIndices,...
        surfaceArray,nSurface ] = getNonDummySurfaceArray( optSystem );
    
    %     nSurface = getNumberOfSurfaces(optSystem);
    %     nNonDummySurface = getNumberOfNonDummySurfaces(optSystem);
    %     NonDummySurfaceArray = getNonDummySurfaceArray(optSystem);
    %     NonDummySurfaceIndices = getNonDummySurfaceIndices(optSystem);
    
    startNonDummyIndex = 1;
    indicesBeforeEndSurf = find(NonDummySurfaceIndices<=endSurface);
    endNonDummyIndex = indicesBeforeEndSurf(end);
    
    % Set info to be displayed on the command window
    dispTIRStatus = 0;
    dispNoIntersectionStatus = 0;
    dispOutOfApertureStatus = 0;
    
    lensUnitFactor = getLensUnitFactor(optSystem);
    wavUnitFactor = getWavelengthUnitFactor(optSystem);
    primaryWavlenInM = getPrimaryWavelength(optSystem);
    
%     nRay = length(objectRayBundle);
    nRay = size(objectRayBundle.Position,2);
    
    wavlenInM  = objectRayBundle.Wavelength; % wavlen is in m for Ray object
    wavlenInWavlenUnit = wavlenInM/wavUnitFactor;
    CurrentRayDirection = objectRayBundle.Direction;
    currentRayPositionInM = objectRayBundle.Position;
%     currentRayPupilCoordinate = objectRayBundle.PupilCoordinate;
    
    
    % Since all ray tracing is done in lens units convert the object position
    % from meter to lens unit.
    CurrentRayPosition = currentRayPositionInM/lensUnitFactor;
    % Initialize the total phase to 0,
    totalPhase = zeros([1,nRay]);
    
    if considerPolarization
        % initialize toatl P, Q and coatingJonesMatrix matrices to unity;
        CoatingJonesMatrix = repmat(eye(2),[1,1,nRay]);
        TotalPMatrix = repmat(eye(3),[1,1,nRay]);
        TotalQMatrix = repmat(eye(3),[1,1,nRay]);
        CoatingPMatrix = repmat(eye(3),[1,1,nRay]);
        CoatingQMatrix = repmat(eye(3),[1,1,nRay]);
    end
    
    % Add initial ray at the first surface to ray trace result structure
    % (object)
    %  All 1 D array are recorded as 3x1 vertical vector
    infThickness = abs(CurrentRayPosition(3,:)) > 10^10;
    CurrentRayPosition(3,infThickness) = 0;
    
    currentSurfaceNormal = repmat([0;0;1],[1,nRay]); % assume plane object surface.
    
    RayIntersectionPoint = CurrentRayPosition;
    ExitRayPosition = CurrentRayPosition;
    IncidentRayDirection = [NaN;NaN;NaN]*ones(1,nRay);%repmat([NaN;NaN;NaN],[1,nRay]); % no ray is incident to object surface
    ExitRayDirection = CurrentRayDirection;
    
    SurfaceNormal = currentSurfaceNormal; % [0,0,1]; % assume plane object surface.
    
    %     IncidenceAngle = NaN([1,nRay]);
    %     ExitAngle = computeAngleBetweenVectors...
    %         (currentSurfaceNormal,CurrentRayDirection);
    
    GeometricalPathLength = zeros([1,nRay]);
    OpticalPathLength = zeros([1,nRay]);
    AdditionalPathLength = zeros([1,nRay]);
    % Absorption can be added here in the future
    
    GroupPathLength = zeros([1,nRay]);
    TotalPathLength = zeros([1,nRay]);
    TotalOpticalPathLength = zeros([1,nRay]);
    TotalGroupPathLength = zeros([1,nRay]);
    
    NoIntersectionPoint = zeros([1,nRay]);
    OutOfAperture = zeros([1,nRay]);
    TotalInternalReflection = zeros([1,nRay]);
    
    RefractiveIndex = ones([1,nRay]);
    RefractiveIndexFirstDerivative  = zeros([1,nRay]);
    RefractiveIndexSecondDerivative  = zeros([1,nRay]);
    GroupRefractiveIndex  = zeros([1,nRay]);
    
    if considerPolarization
        if recordIntermediateResults
            multipleRayTracerResultAll(startNonDummyIndex,:) = RayTraceResult(nRayPupil,nField,nWav,...
                RayIntersectionPoint,ExitRayPosition,SurfaceNormal,...
                IncidentRayDirection,ExitRayDirection,NoIntersectionPoint,OutOfAperture,TotalInternalReflection,GeometricalPathLength,AdditionalPathLength,OpticalPathLength,...
                GroupPathLength,TotalPathLength,TotalOpticalPathLength,TotalGroupPathLength,...
                RefractiveIndex,RefractiveIndexFirstDerivative,RefractiveIndexSecondDerivative,...
                GroupRefractiveIndex,CoatingJonesMatrix,CoatingPMatrix,CoatingQMatrix,TotalPMatrix,TotalQMatrix);
%             multipleRayTracerResultAll(startNonDummyIndex,:).RayPupilCoordinates = currentRayPupilCoordinate;
        else
            multipleRayTracerResultFinal(startNonDummyIndex,:) = RayTraceResult(nRayPupil,nField,nWav,...
                RayIntersectionPoint,ExitRayPosition,SurfaceNormal,...
                IncidentRayDirection,ExitRayDirection,NoIntersectionPoint,OutOfAperture,TotalInternalReflection,GeometricalPathLength,AdditionalPathLength,OpticalPathLength,...
                GroupPathLength,TotalPathLength,TotalOpticalPathLength,TotalGroupPathLength,...
                RefractiveIndex,RefractiveIndexFirstDerivative,RefractiveIndexSecondDerivative,...
                GroupRefractiveIndex,CoatingJonesMatrix,CoatingPMatrix,CoatingQMatrix,TotalPMatrix,TotalQMatrix);
%             multipleRayTracerResultFinal(startNonDummyIndex,:).RayPupilCoordinates = currentRayPupilCoordinate;
        end
    else
        if recordIntermediateResults
            multipleRayTracerResultAll(startNonDummyIndex,:) = RayTraceResult(nRayPupil,nField,nWav,...
                RayIntersectionPoint,ExitRayPosition,SurfaceNormal,...
                IncidentRayDirection,ExitRayDirection,NoIntersectionPoint,OutOfAperture,TotalInternalReflection,GeometricalPathLength,AdditionalPathLength,OpticalPathLength,...
                GroupPathLength,TotalPathLength,TotalOpticalPathLength,TotalGroupPathLength,...
                RefractiveIndex,RefractiveIndexFirstDerivative,RefractiveIndexSecondDerivative,...
                GroupRefractiveIndex);
%              multipleRayTracerResultAll(startNonDummyIndex,:).RayPupilCoordinates = currentRayPupilCoordinate;
        else
            multipleRayTracerResultFinal(startNonDummyIndex,:) = RayTraceResult(nRayPupil,nField,nWav,...
                RayIntersectionPoint,ExitRayPosition,SurfaceNormal,...
                IncidentRayDirection,ExitRayDirection,NoIntersectionPoint,OutOfAperture,TotalInternalReflection,GeometricalPathLength,AdditionalPathLength,OpticalPathLength,...
                GroupPathLength,TotalPathLength,TotalOpticalPathLength,TotalGroupPathLength,...
                RefractiveIndex,RefractiveIndexFirstDerivative,RefractiveIndexSecondDerivative,...
                GroupRefractiveIndex);
%              multipleRayTracerResultFinal(startNonDummyIndex,:).RayPupilCoordinates = currentRayPupilCoordinate;
        end
    end
    
    
    
    firstGlass = NonDummySurfaceArray(1).Glass;
    
    % compute refractive indices
    if computeGroupPathLength || computeGroupIndex
        [groupIndexAfter,indexAfter, firstDerivative_IndexAfter] = getGroupRefractiveIndex(firstGlass,wavlenInM) ;
    elseif computeRefractiveIndexFirstDerivative
        groupIndexAfter = zeros([1,nRay]);
        indexAfter = getRefractiveIndex(firstGlass,wavlenInM);
        firstDerivative_IndexAfter = getRefractiveIndex(firstGlass,wavlenInM,1);
    else
        groupIndexAfter = zeros([1,nRay]);
        indexAfter = getRefractiveIndex(firstGlass,wavlenInM);
        firstDerivative_IndexAfter = zeros([1,nRay]);
    end
    if computeRefractiveIndexSecondDerivative
        secondDerivative_IndexAfter = getRefractiveIndex(firstGlass,wavlenInM,2);
    else
        secondDerivative_IndexAfter = zeros([1,nRay]);
    end
    
    
    %     if computeGroupPathLength
    %         [groupIndexAfter,indexAfter, firstDerivative_IndexAfter] = getGroupRefractiveIndex(firstGlass,wavlenInM) ;
    %     else
    %         groupIndexAfter = NaN;
    %         indexAfter = getRefractiveIndex(firstGlass,wavlenInM);
    %         firstDerivative_IndexAfter = NaN;
    %     end
    
    % To keep track of mirrored coordinates (after mirror occuring odd times)
    mirroredCoordinate = 0;
    
    for surfaceIndex = startNonDummyIndex+1:1:endNonDummyIndex
        surfaceType = NonDummySurfaceArray(surfaceIndex).Type;
        % Transfer current ray location, direction and polarization vector
        % to local coordinates corresponding to kth surface
        globalRayExitPosition = CurrentRayPosition;
        globalRayDirection = CurrentRayDirection;
        
        %% New Code for Global to Local Coordinate transformation using Coordinate Transfer Matrix
        prevThickness = NonDummySurfaceArray(surfaceIndex-1).Thickness;
        if abs(prevThickness) > 10^10
            prevThickness = 0;
        end
        surfaceCoordinateTM = NonDummySurfaceArray(surfaceIndex).SurfaceCoordinateTM;
        [localRayExitPosition,localRayDirection] = globalToLocalCoordinate(...
            globalRayExitPosition,globalRayDirection,surfaceCoordinateTM);
        
        %% End of New code
        
        % from this point onwards  currentRay become localRay
        CurrentRayPosition = localRayExitPosition;
        CurrentRayDirection = localRayDirection;
        
        % Real ray trace
        glassAfter = NonDummySurfaceArray(surfaceIndex).Glass;
        
        indexBefore = indexAfter;
%         firstDerivative_IndexBefore = firstDerivative_IndexAfter;
%         secondDerivative_IndexBefore = secondDerivative_IndexAfter;
        groupIndexBefore = groupIndexAfter ;
        %         indexAfter = getRefractiveIndex(glassAfter,wavlenInM);
        %
        %         firstDerivative_IndexAfter = getRefractiveIndex(glassAfter,wavlenInM,1);
        %         secondDerivative_IndexAfter = getRefractiveIndex(glassAfter,wavlenInM,2);
        
        %         [groupIndexAfter,indexAfter, firstDerivative_IndexAfter] = getGroupRefractiveIndex(glassAfter,wavlenInM);
        
        % compute refractive indices
        if computeGroupPathLength || computeGroupIndex
            [groupIndexAfter,indexAfter, firstDerivative_IndexAfter] = getGroupRefractiveIndex(glassAfter,wavlenInM) ;
        elseif computeRefractiveIndexFirstDerivative
            groupIndexAfter = zeros([1,nRay]);
            indexAfter = getRefractiveIndex(glassAfter,wavlenInM);
            firstDerivative_IndexAfter = getRefractiveIndex(glassAfter,wavlenInM,1);
        else
            groupIndexAfter = zeros([1,nRay]);
            indexAfter = getRefractiveIndex(glassAfter,wavlenInM);
            firstDerivative_IndexAfter = zeros([1,nRay]);
        end
        if computeRefractiveIndexSecondDerivative
            secondDerivative_IndexAfter = getRefractiveIndex(glassAfter,wavlenInM,2);
        else
            secondDerivative_IndexAfter = zeros([1,nRay]);
        end
        
        
        %
        %     if computeGroupPathLength
        %         [groupIndexAfter,indexAfter, firstDerivative_IndexAfter] = getGroupRefractiveIndex(glassAfter,wavlenInM) ;
        %     else
        %         groupIndexAfter = NaN;
        %         indexAfter = getRefractiveIndex(glassAfter,wavlenInM);
        %         firstDerivative_IndexAfter = NaN;
        %     end
        %
        
        rayInitialPosition = CurrentRayPosition;
        rayDirection = CurrentRayDirection;
        
        if  strcmpi(glassAfter.Name,'Mirror') % reflection
            reflection = 1;
        else
            reflection = 0;
        end
        
        [ GeometricalPathLength,additionalPath,localRayIntersectionPoint, ...
            localSurfaceNormal,localExitRayPosition,localExitRayDirection,...
            TIR,NoIntersectionPoint]   = traceRealRaysToThisSurface(...
            NonDummySurfaceArray(surfaceIndex),rayInitialPosition,rayDirection,...
            indexBefore,indexAfter,wavlenInM,primaryWavlenInM,reflection);
        
        AdditionalPathLength = additionalPath;
        
        TotalPathLength = TotalPathLength + GeometricalPathLength;
        
        if computeOpticalPathLength
            OpticalPathLength = indexBefore.*(GeometricalPathLength + additionalPath);
            TotalOpticalPathLength = TotalOpticalPathLength + OpticalPathLength;
        else
            OpticalPathLength = zeros([1,nRay]);
            TotalOpticalPathLength = zeros([1,nRay]);
        end
        
        if computeGroupPathLength
            GroupPathLength = groupIndexBefore.*(GeometricalPathLength + additionalPath);
            TotalGroupPathLength = TotalGroupPathLength + GroupPathLength;
        else
            GroupPathLength = zeros([1,nRay]);
            TotalGroupPathLength = zeros([1,nRay]);
        end

        if dispNoIntersectionStatus
            totalNoIntersection = sum(NoIntersectionPoint);
            if totalNoIntersection > 0
                disp([num2str(totalNoIntersection) , ' rays do not intersect surface ',...
                    num2str(surfaceIndex),'. They have no intersection point. They are',...
                    ' ignored in other computations.']);
            else
                disp(['All of ',num2str(nRay),' rays intersect surface ',num2str(surfaceIndex),'.']);
            end
        end
        
        OutOfAperture = zeros([1,nRay]);
        if considerSurfAperture
            % check whether the intersection point is inside the aperture
            surfAperture = NonDummySurfaceArray(surfaceIndex).Aperture;
            pointX = localRayIntersectionPoint(1,:);
            pointY = localRayIntersectionPoint(2,:);
            xyVector = [pointX',pointY'];
            [ isInsideTheMainAperture ] = IsInsideTheMainAperture( surfAperture, xyVector );
            
            OutOfAperture = ~isInsideTheMainAperture';
        end
        
        totalOutOfAperture = sum(double(OutOfAperture));
        if dispOutOfApertureStatus
            if totalOutOfAperture > 0
                disp([num2str(totalOutOfAperture) , ' rays are vignated at surface ',...
                    num2str(surfaceIndex),'. They intersect the surface outside the ' ,...
                    'aperture area. They are ignored in other computations.']);
            else
                disp(['All of ',num2str(nRay),' rays pass through surface ',num2str(surfaceIndex),'.']);
            end
        end
        
        CurrentRayPosition = localExitRayPosition;
        
        % If the current coordinate is mirrored coordinate
        % => Mirrored Z axis => the normal vector shall be in -ve direction
        if mirroredCoordinate
            localSurfaceNormal = -localSurfaceNormal;
        end
        localIncidentRayDirection = CurrentRayDirection;
        
        % localIncidenceAngle = computeAngleBetweenVectors (localSurfaceNormal,...
        %     localIncidentRayDirection); % new function replacing yi's 5th function
        % localIncidenceAngle = zeros(1,nRay);
        
        if surfaceIndex < nNonDummySurface
            if considerPolarization
                coatingType = NonDummySurfaceArray(surfaceIndex).Coating.Type;
            end
            Kqm1 = localIncidentRayDirection;
            
            CurrentRayDirection = localExitRayDirection;
            
            % localExitAngle = computeAngleBetweenVectors(localSurfaceNormal,...
            %     localExitRayDirection);
            % localExitAngle = zeros(1,nRay);
            
            Kq = CurrentRayDirection;
            nTotIR = sum(TIR);
            if dispTIRStatus
                if nTotIR > 0
                    disp([num2str(nTotIR) , ' rays encountered Total Internal Reflections at ' ,...
                        'surface ', num2str(surfaceIndex),'.']);
                else
                    disp(['All of ',num2str(nRay),' rays refracted through surface ',...
                        num2str(surfaceIndex),' with out Total Internal Reflections.']);
                end
            end
            wavLenInUm = wavlenInM*10^6;
            if strcmpi(NonDummySurfaceArray(surfaceIndex).Glass.Name,'Mirror')
                if considerPolarization
                    localIncidenceAngle = computeAngleBetweenVectors (localSurfaceNormal,...
                        localIncidentRayDirection);
                    % jones matrix for multiple rays with localIncidenceAngle
                    % Wavelength shall be converted to um
                    
                    primaryWaveLenInUm = primaryWavlenInM*10^6;
                    [ampRs,ampRp,powRs,powRp,jonesMatrix]=...
                        getReflectionCoefficients(NonDummySurfaceArray(surfaceIndex).Coating,...
                        wavLenInUm,localIncidenceAngle,...
                        indexBefore,indexAfter,primaryWaveLenInUm);
                    
                end
                mirroredCoordinate = ~mirroredCoordinate;
            elseif ~strcmpi(NonDummySurfaceArray(surfaceIndex).Glass.Name,'Mirror')
                if considerPolarization
                    localIncidenceAngle = computeAngleBetweenVectors (localSurfaceNormal,...
                        localIncidentRayDirection);
                    % new code
                    % Wavelength shall be converted to um
                    
                    primaryWaveLenInUm = primaryWavlenInM*10^6;
                    [ampTs,ampTp,powTs,powTp,jonesMatrix] = ...
                        getTransmissionCoefficients(...
                        NonDummySurfaceArray(surfaceIndex).Coating, ...
                        wavLenInUm,localIncidenceAngle,...
                        indexBefore,indexAfter,primaryWaveLenInUm);
                end
            end
            
            RefractiveIndex = indexAfter;
            RefractiveIndexFirstDerivative = firstDerivative_IndexAfter;
            RefractiveIndexSecondDerivative = secondDerivative_IndexAfter;
            %  RefractiveIndexSecondDerivative = zeros(1,nRay);
            GroupRefractiveIndex = groupIndexAfter;
            
            % compute the new P matrix
            if considerPolarization
                CoatingJonesMatrix = jonesMatrix;
                CoatingPMatrix = convertJonesMatrixToPolarizationMatrix(jonesMatrix,Kqm1,Kq);
                % compute the total P matrix
                TotalPMatrix = multiplySliced3DMatrices(CoatingPMatrix,TotalPMatrix);
                % compute the parallel ray transport matrix Q
                qJonesMatrix = repmat([1,0;0,1],[1,1,nRay]);   % unity jones matrix for non polarizing element
                CoatingQMatrix = convertJonesMatrixToPolarizationMatrix(qJonesMatrix,Kqm1,Kq);
                % compute the total Q matrix
                TotalQMatrix = multiplySliced3DMatrices(CoatingQMatrix,TotalQMatrix);
            end
        elseif surfaceIndex==nNonDummySurface
            % for image surface no refraction
            %             localExitRayDirection  = localIncidentRayDirection;
            %             localExitAngle = localIncidenceAngle;
            if considerPolarization
                CoatingJonesMatrix = repmat(eye(2),[1,1,nRay]);
                CoatingPMatrix = repmat(eye(3),[1,1,nRay]);
                CoatingQMatrix = repmat(eye(3),[1,1,nRay]);
            end
        end
        % Finally Transfer back to Global coordinate system
        % Local to Global Coordinate transformation
        [GlobalRayIntersectionPoint,GlobalExitRayPosition,GlobalSurfaceNormal,...
            GlobalIncidentRayDirection,GlobalExitRayDirection] = ...
            localToGlobalCoordinate(localRayIntersectionPoint,localExitRayPosition,...
            localSurfaceNormal,localIncidentRayDirection,localExitRayDirection,surfaceCoordinateTM);
        
        %         globalIncidenceAngle = localIncidenceAngle;
        %         globalExitAngle = localExitAngle;
        
        % Add signs to the angles
        % +Ve angles: CCW when observed from (SurfaceNormal X RayDirection)
        % from direction with -ve x component
        % -Ve angles: CW when observed from (SurfaceNormal X RayDirection)
        % from direction with -ve x component
        % This makes the angles valid for rays in yz planes with conventional sign convention.
        % normalToPlaneOfPropagationIncident =  ...
        %     compute3dCross(GlobalSurfaceNormal,GlobalIncidentRayDirection);
        %
        % normalToPlaneOfPropagationExit =  ...
        %     compute3dCross(GlobalSurfaceNormal,GlobalExitRayDirection);
        
        %     if strcmpi(NonDummySurfaceArray(surfaceIndex).Glass.Name,'Mirror')
        %         % Invert the incidence direction in case of mirrror
        %         %if mirroredCoordinate
        %         signOfIncidence = -1;
        %     else
        %         signOfIncidence = 1;
        %     end
        
        % GlobalIncidenceAngleSigned = globalIncidenceAngle.* ...
        %     -sign(normalToPlaneOfPropagationIncident(1));
        % GlobalIncidenceAngleSigned = globalIncidenceAngle;
        
        % GlobalExitAngleSigned = globalExitAngle.* ...
        %     -sign(normalToPlaneOfPropagationExit(1));
        % GlobalExitAngleSigned = globalExitAngle;
        
        % Record all neccessary outputs to trace the ray
        if recordIntermediateResults
            if considerPolarization
                multipleRayTracerResultAll(surfaceIndex,:) = RayTraceResult(nRayPupil,nField,nWav,...
                    GlobalRayIntersectionPoint,GlobalExitRayPosition,GlobalSurfaceNormal,...
                    GlobalIncidentRayDirection,GlobalExitRayDirection,...
                    NoIntersectionPoint,OutOfAperture,TIR,GeometricalPathLength,AdditionalPathLength,OpticalPathLength,...
                    GroupPathLength,TotalPathLength,TotalOpticalPathLength,TotalGroupPathLength,...
                    RefractiveIndex,RefractiveIndexFirstDerivative,RefractiveIndexSecondDerivative,...
                    GroupRefractiveIndex,CoatingJonesMatrix,CoatingPMatrix,CoatingQMatrix,TotalPMatrix,TotalQMatrix);
            else
                multipleRayTracerResultAll(surfaceIndex,:) = RayTraceResult(nRayPupil,nField,nWav,...
                    GlobalRayIntersectionPoint,GlobalExitRayPosition,GlobalSurfaceNormal,...
                    GlobalIncidentRayDirection,GlobalExitRayDirection,NoIntersectionPoint,...
                    OutOfAperture,TIR,GeometricalPathLength,AdditionalPathLength,OpticalPathLength,...
                    GroupPathLength,TotalPathLength,TotalOpticalPathLength,TotalGroupPathLength,...
                    RefractiveIndex,RefractiveIndexFirstDerivative,RefractiveIndexSecondDerivative,...
                    GroupRefractiveIndex);
            end
        else
            if surfaceIndex == endNonDummyIndex
                if considerPolarization
                    multipleRayTracerResultFinal(2,:) = RayTraceResult(nRayPupil,nField,nWav,...
                        GlobalRayIntersectionPoint,GlobalExitRayPosition,GlobalSurfaceNormal,GlobalIncidentRayDirection,...
                        GlobalExitRayDirection,NoIntersectionPoint,OutOfAperture,TIR,GeometricalPathLength,AdditionalPathLength,OpticalPathLength,...
                        GroupPathLength,TotalPathLength,TotalOpticalPathLength,TotalGroupPathLength,...
                        RefractiveIndex,RefractiveIndexFirstDerivative,RefractiveIndexSecondDerivative,...
                        GroupRefractiveIndex,CoatingJonesMatrix,CoatingPMatrix,CoatingQMatrix,TotalPMatrix,TotalQMatrix);
                else
                    multipleRayTracerResultFinal(2,:) = RayTraceResult(nRayPupil,nField,nWav,...
                        GlobalRayIntersectionPoint,GlobalExitRayPosition,GlobalSurfaceNormal,...
                        GlobalIncidentRayDirection,GlobalExitRayDirection,NoIntersectionPoint,...
                        OutOfAperture,TIR,GeometricalPathLength,AdditionalPathLength,OpticalPathLength,...
                        GroupPathLength,TotalPathLength,TotalOpticalPathLength,TotalGroupPathLength,...
                        RefractiveIndex,RefractiveIndexFirstDerivative,RefractiveIndexSecondDerivative,...
                        GroupRefractiveIndex);
                end
            end
        end
        
        % now the currentnt ray becomes global ray agin
        CurrentRayPosition = GlobalExitRayPosition;
        CurrentRayDirection = GlobalExitRayDirection;
        if considerSurfAperture
            % Make the localRayIntersectionPoint = NaN for those rays outside the
            % aperture
            CurrentRayPosition(:,OutOfAperture) = NaN;
            CurrentRayDirection(:,OutOfAperture) = NaN;
        end
    end
    % The ray trace result of the non dummy surfaces only.
    if recordIntermediateResults
        rayTracerResult = multipleRayTracerResultAll; % All results
    else
        rayTracerResult = multipleRayTracerResultFinal; % only 1st and last surface
    end
    
    rayTracerResultReshaped = rayTracerResult;
%     profile viewer
end
