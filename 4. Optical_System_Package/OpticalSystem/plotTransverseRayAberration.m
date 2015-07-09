function plotTransverseRayAberration(optSystem,surfIndex,wavLen,...
    fieldPointXY,numberOfRays1,numberOfRays2,sagittalAberrComp,tangentialAberrComp,...
    plotPanelHandle)
% Displays the transverse ray aberration of sagittal and tangetial ray fans
% on any surface with respect to the cheif ray. 

% <<<<<<<<<<<<<<<<<<<<<<<<< Author Section >>>>>>>>>>>>>>>>>>>>>>>>>>>>
%   Written By: Worku, Norman Girma
%   Advisor: Prof. Herbert Gross
%   Part of the RAYTRACE_TOOLBOX V3.0 (OOP Version)
%	Optical System Design and Simulation Research Group
%   Institute of Applied Physics
%   Friedrich-Schiller-University of Jena

% <<<<<<<<<<<<<<<<<<< Change History Section >>>>>>>>>>>>>>>>>>>>>>>>>>
% Date----------Modified By ---------Modification Detail--------Remark
% May 21,2014   Worku, Norman G.     Original Version      As example of extension

% <<<<<<<<<<<<<<<<<<<<< Main Code Section >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Default Inputs
if nargin < 8
    disp('Error: The function requires atleast 6 arguments, optSystem,',...
        ' surfIndex, wavLen, fieldPointXY, numberOfRays1, numberOfRays2',...
        'sagittalAberrComp and tangentialAberrComp.');
    return;
elseif nargin == 8
    axesHandle = axes('Parent',figure,'Units','normalized',...
        'Position',[0.1,0.1,0.8,0.8]);
else
    
end
NonDummySurfaceIndices = optSystem.getNonDummySurfaceIndices;
if ~find(ismember(NonDummySurfaceIndices,surfIndex))
    disp('Error:No ray trace data for Dummy surfaces');
    return;    
end
nonDummySurfacesBeforeCurrentSurface = sum(NonDummySurfaceIndices<=surfIndex);
dummySurfacesBeforeCurrentSurface = surfIndex - nonDummySurfacesBeforeCurrentSurface;        
surfIndexWithOutDummy = surfIndex-dummySurfacesBeforeCurrentSurface;
% Assign different symbals and colors for lines of d/t wavelengths
% and feild points respectively
availablelineColor = repmat(['b','k','r','g','c','m','y'],1,20); % 7*20 = 140
lineColorList = availablelineColor(1:size(wavLen,2)*size(fieldPointXY,2));

%cla(axesHandle,'reset')
delete(allchild(plotPanelHandle));

% polarizedRayTracerResult =  nSurfWithoutDummy X nRay X nField X nWav
[sagittalRayTracerResult] = multipleRayTracer(optSystem,wavLen,...
    fieldPointXY,numberOfRays1,numberOfRays2,'Sagittal');
[tangentialRayTracerResult] = multipleRayTracer(optSystem,wavLen,...
    fieldPointXY,numberOfRays1,numberOfRays2,'Tangential');

nSurfaceWithOutDummy = size(sagittalRayTracerResult,1);
nRayTotalSagittal = size(sagittalRayTracerResult,2);
nRayTotalTangential = size(tangentialRayTracerResult,2);
nField = size(sagittalRayTracerResult,3);
nWav = size(sagittalRayTracerResult,4);
% trace the cheif ray with either primary wavelength (for multiple 
% wavelength analysis)or the specified wavelength (for single wavelegnth).
if nWav > 1
    % Use the primary wavelength for the cheif ray 
    cheifRayWavLenInMet = getPrimaryWavelength(optSystem);
%     % Change wavlegth back to the system wavelength unit
%     wavUnitFactor = optSystem.getWavelengthUnitFactor;
%     cheifRayWavLenInM = cheifRayWavLenInMet/wavUnitFactor;
elseif nWav == 1 
    % Use the specified wavelength for the cheif ray 
    cheifRayWavLenInMet = wavLen*getPrimaryWavelength(optSystem);
else
end
considerSurfAperture = 0;
% cheifRayTraceResult: nSurf X nField as each field point has different
% cheif rays
chiefRayTraceResult = traceChiefRay(optSystem,fieldPointXY,cheifRayWavLenInMet);

% All raytracing are done with out recording intermediate results so the
% index of last surface will be 2
rayTraceResultIndexOfLastSurface = 2;

% Use different color for diffrent wavelengths and different field points.         
for wavIndex = 1:nWav
    for fieldIndex = 1:nField
        lineIndex = fieldIndex + (wavIndex-1)*nField;
        
        cheifRayIntersection = [chiefRayTraceResult(rayTraceResultIndexOfLastSurface,fieldIndex).RayIntersectionPoint];
        
        cheifRayIntersectionsSagittal = repmat(cheifRayIntersection,[1,nRayTotalSagittal]);
        cheifRayIntersectionsTangential = repmat(cheifRayIntersection,[1,nRayTotalTangential]);
        sagittalFanIntersectionPoints = ...
            [sagittalRayTracerResult(rayTraceResultIndexOfLastSurface,:,fieldIndex,wavIndex).RayIntersectionPoint];
        tangentialFanIntersectionPoints = ...
            [tangentialRayTracerResult(rayTraceResultIndexOfLastSurface,:,fieldIndex,wavIndex).RayIntersectionPoint];

        if strcmpi(sagittalAberrComp,'X Aberration')
            sagY(lineIndex,:) = sagittalFanIntersectionPoints(1,:) - cheifRayIntersectionsSagittal(1,:);
            yLabelSag = 'EX';
        elseif strcmpi(sagittalAberrComp,'Y Aberration')
            sagY(lineIndex,:) = sagittalFanIntersectionPoints(2,:) - cheifRayIntersectionsSagittal(2,:);
            yLabelSag = 'EY';
        else
        end
            
        if strcmpi(tangentialAberrComp,'X Aberration')
            tanY(lineIndex,:) = tangentialFanIntersectionPoints(1,:) - cheifRayIntersectionsTangential(1,:);
            yLabelTan = 'EX';
        elseif strcmpi(tangentialAberrComp,'Y Aberration')
            tanY(lineIndex,:) = tangentialFanIntersectionPoints(2,:) - cheifRayIntersectionsTangential(2,:);
            yLabelTan = 'EY';
        else
        end           
        
        xLabelTan = 'PY';
        xLabelSag = 'PX';
        sagX(lineIndex,:) = linspace(-1,1,size(sagY,2));
        tanX(lineIndex,:) = linspace(-1,1,size(tanY,2));

        legendText{lineIndex} = ['Field: [',num2str(fieldPointXY(1,fieldIndex)),',',...
            num2str(fieldPointXY(2,fieldIndex)),']',...
            ' Wav: ',num2str(wavLen(wavIndex))];
    end
end

% Generate two new panel for sagittal and tangential fans.
tangentialPlotPanel = uipanel('Parent',plotPanelHandle,...
    'Units','Normalized',...
    'Position',[0.52,0.1,0.45,0.8],...
    'Title',[char(tangentialAberrComp),' for Tangential Fan']);
tangentialPlotAxes = axes('Parent',tangentialPlotPanel,...
    'Units','Normalized',...
    'Position',[0.1,0.2,0.88,0.6]);

for tanKK = 1:lineIndex
    currentLineColor = lineColorList(tanKK);
    plot(tangentialPlotAxes,tanX(tanKK,:),tanY(tanKK,:),currentLineColor);        
    hold on;
end
grid on;
xlabel(tangentialPlotAxes,xLabelTan,'FontSize',12);
ylabel(tangentialPlotAxes,yLabelTan,'FontSize',12);
legend(tangentialPlotAxes,legendText)

sagittalPlotPanel = uipanel('Parent',plotPanelHandle,...
    'Units','Normalized',...
    'Position',[0.03,0.1,0.45,0.8],...
    'Title',[char(sagittalAberrComp),' for Sagittal Fan']);

sagittalPlotAxes = axes('Parent',sagittalPlotPanel,...
    'Units','Normalized',...
    'Position',[0.1,0.2,0.88,0.6]);
        
for sagKK = 1:lineIndex
    currentLineColor = lineColorList(sagKK);
    plot(sagittalPlotAxes,sagX(sagKK,:),sagY(sagKK,:),currentLineColor);
    hold on;
end
grid on;
xlabel(sagittalPlotAxes,xLabelSag,'FontSize',12);
ylabel(sagittalPlotAxes,yLabelSag,'FontSize',12);
legend(sagittalPlotAxes,legendText)

end