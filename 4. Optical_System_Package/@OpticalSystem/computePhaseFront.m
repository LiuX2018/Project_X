function [ meshGridX,meshGridY,phaseFront] =...
    computePhaseFront( optSystem,wavLenInWavUnit,fieldPointXYInLensUnit,refSurfaceIndex,...
         offsetInTime, nSamplingPoints1,nSamplingPoints2,gridType)
% locationOfPhaseFrontAlongChiefRay: In lens unit
% gridType: Polar or Cartesian used as pupil sampling type
%  perform scalar/polarization ray trace depending on  jonesVec

c = 299792458;
PupSamplingType = gridType;
considerSurfAperture = 1;
recordIntermediateResults = 1;
considerPolarization = 0;

[rayTracerResult,pupilMeshGrid,outsidePupilIndices] = ...
    multipleRayTracer(optSystem,wavLenInWavUnit,fieldPointXYInLensUnit,nSamplingPoints1,...
    nSamplingPoints2,PupSamplingType,considerPolarization,considerSurfAperture,...
    recordIntermediateResults);
X = pupilMeshGrid(:,:,1);
Y = pupilMeshGrid(:,:,2);
if isempty(outsidePupilIndices)
    % Find cheif ray index
    cheifRayIndex = find((X).^2+(Y).^2==0);   
else
    % Find cheif ray index
    cheifRayIndex = find((X(~outsidePupilIndices)).^2+(Y(~outsidePupilIndices)).^2==0);
end

% take first value in case of multiple ocuurances
cheifRayIndex = cheifRayIndex(1);

% Determine the surface just before the requested distance

totalPathLengthVector = [rayTracerResult(:,cheifRayIndex).TotalPathLength];
% surfacesBeforeCurrentPoint = find(totalPathLengthVector <= locationOfPhaseFrontAlongChiefRay);

%prevSurfaceIndex = surfacesBeforeCurrentPoint(end);
%refractiveIndexOfLastMedium = optSystem.getSurfaceArray(prevSurfaceIndex).Glass.getRefractiveIndex(wavLenInSI);

% 
if offsetInTime == 0
   prevSurfaceIndex = refSurfaceIndex; 
elseif offsetInTime > 0 % 
   prevSurfaceIndex = refSurfaceIndex;
else
    prevSurfaceIndex = refSurfaceIndex-1; 
end

wavLenInSI = wavLenInWavUnit * getWavelengthUnitFactor(optSystem);
currentRefractiveIndex = optSystem.getSurfaceArray(prevSurfaceIndex).Glass.getRefractiveIndex(wavLenInSI);
locationOfPhaseFrontAlongChiefRay = totalPathLengthVector(refSurfaceIndex)+...
    (c*offsetInTime)/(getLensUnitFactor(optSystem));

% Compute the points on each ray with equal optical pathlength
totalOpticalPathLengthToPhaseFront = (rayTracerResult(prevSurfaceIndex,cheifRayIndex).TotalOpticalPathLength)+...
    (locationOfPhaseFrontAlongChiefRay-rayTracerResult(prevSurfaceIndex,cheifRayIndex).TotalPathLength)*currentRefractiveIndex;
distanceFromLastSurfaceToPhaseFront = (totalOpticalPathLengthToPhaseFront-...
    [rayTracerResult(prevSurfaceIndex,:).TotalOpticalPathLength])/currentRefractiveIndex;
pointsOfThePhaseFront =  [rayTracerResult(prevSurfaceIndex,:).RayIntersectionPoint]+...
    repmat(distanceFromLastSurfaceToPhaseFront,[3,1]).*[rayTracerResult(prevSurfaceIndex,:).ExitRayDirection];


%% Plot the phase front
% Interpolate points to the mesh grid for surf plot
x = pointsOfThePhaseFront(1,:);
y = pointsOfThePhaseFront(2,:);
z = pointsOfThePhaseFront(3,:);

    switch lower(gridType)
        case lower({'Tangential','Sagittal'})
            Xq = x;
            Yq = y;
            Zq = z;   
        case lower({'Cartesian','Rectangular'})
            nPointsX = nSamplingPoints1;
            nPointsY = nSamplingPoints2;        
            if nPointsX == 1
                xgv = 0;
            else
                xgv = linspace(min(x),max(x),nPointsX);
            end
            if nPointsY == 1
                ygv = 0;
            else
                ygv = linspace(min(y),max(y),nPointsY);
            end
            [Xq,Yq] = meshgrid(xgv,ygv);
          
        case lower('Polar')
            nPointsRadial = nSamplingPoints1;
            nPointsAngular = nSamplingPoints2;
            % Radius of the largest circle circumscribing the aperture
            maxR = max(sqrt(((x).^2+(y).^2)));
            % Draw a circle with radiaus maxR and then cut out the part required
            % using the given X and Y ranges
            if nPointsRadial == 1
                r = 0;
            else
                r = (linspace(-maxR,maxR,nPointsRadial*2+1))';
            end
            if nPointsAngular == 1
                phi = 0;
            else
                phi = linspace(0,2*pi,nPointsAngular);
            end
            Xq = r*cos(phi);
            Yq = r*sin(phi);
    end

Zq = griddata(x,y,z,Xq,Yq);

if isempty(Zq)
    % The underlying triangulation is empty - the points may be collinear. 
    meshGridX = x;
    meshGridY = y;
    phaseFront = z;   
else
    % Return the phase front surface in an interpolated grid mesh
    meshGridX = Xq;
    meshGridY = Yq;
    phaseFront = Zq;
end
end

