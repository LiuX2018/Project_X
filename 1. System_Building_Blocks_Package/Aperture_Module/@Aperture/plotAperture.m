function [  h ] = plotAperture( surfAperture,nPoints1,nPoints2,xyCenterPoint,axesHandle )
    %PLOTAPERTURE Plots the aperture in the given axis
    % Three regions are colored differently
    % The actual aperture region: White, The outer aperture region: gray
    % and the edge region : dark
    
    [outerApertShape,maximumRadiusXY] =  getOuterAperture( surfAperture );
    % Radius of the largest circle circumscribing the aperture
        maxR = sqrt((maximumRadiusXY(1))^2+(maximumRadiusXY(2))^2);
        maxRWithEdge = (1+surfAperture.AdditionalEdge)*maxR;
    if strcmpi(outerApertShape,'Circular')||strcmpi(outerApertShape,'Elliptical')  
        % Draw a circle with radiaus maxR and then cut out the part required
        % using the given X and Y ranges
        %     r = (linspace(-maxR,maxR,nPoints1))';
        r = (linspace(-maxRWithEdge,maxRWithEdge,nPoints1))';
        phi = (linspace(-pi,pi,nPoints2));
        xMesh = r*cos(phi);
        yMesh = r*sin(phi);
    else % rectangular
        xgv = linspace(-maxRWithEdge,maxRWithEdge,nPoints1);
        ygv = linspace(-maxRWithEdge,maxRWithEdge,nPoints2);
        [xMesh,yMesh] = meshgrid(xgv,ygv);
     end
        nRow = size(xMesh,1);
    nCol = size(xMesh,2);
    % The x and y points will be completely in unrotated and undecentered
    % aperture coordinate which is assumed by IsInsideTheOuterAperture function
    % below.
    xyVector = [xMesh(:),yMesh(:)];

    [ isInsideTheOuterAperture ] = IsInsideTheOuterAperture( surfAperture, xyVector );
    [ isInsideTheMainAperture ] = IsInsideTheMainAperture( surfAperture, xyVector );
    [ isInsideTheAdditionalEdge ] = IsInsideTheAdditionalEdge( surfAperture, xyVector );   
    colorData = yMesh*NaN;
    % Remove those points outside aperture + edge
    isInsideTheOuterApertureWithEdge = isInsideTheOuterAperture | isInsideTheAdditionalEdge;
        colorData(isInsideTheOuterAperture) = 1;
        colorData(isInsideTheMainAperture) = 0;
        colorData(isInsideTheAdditionalEdge) = 0.5;
        
    xMesh(~isInsideTheOuterApertureWithEdge) = NaN;
    yMesh(~isInsideTheOuterApertureWithEdge) = NaN;
    colorData(~isInsideTheOuterApertureWithEdge) = NaN;
    
    % Remove the all NaN rows and columns
    insideDrawn = reshape(isInsideTheOuterApertureWithEdge,[nRow,nCol]);
    xMesh(:,~any(insideDrawn,1)) = [];
    yMesh(:,~any(insideDrawn,1)) = [];
    colorData(:,~any(insideDrawn,1)) = [];
    
    xMesh(~any(insideDrawn,2),:) = [];
    yMesh(~any(insideDrawn,2),:) = [];
    colorData(~any(insideDrawn,2),:) = [];
    
    % Replace the zeros in colorData wihth NaN so that they will be drwn
    % white
    colorData(colorData==0) = NaN;
    h = pcolor(xMesh,yMesh,colorData);
    set(h, 'EdgeColor', 'none');
end

