function NewGaussianBeam = ScalarGaussianBeam(centralRay,waistRadiusInX,...
        waistRadiusInY,distanceFromWaistInX,distanceFromWaistInY,...
        peakAmplitude,localXDirection,localYDirection)
    % ScalarGaussianBeam Struct:
    %
    %   To represent all scalar (with no polarization) gaussian beam objects
    %   The class supports constructors to construct an array of gaussian beam
    %   objects from array of its properties.
    %
    % Properties: 7 and methods: 0
    %
    % Example Calls:
    %
    % newScalarGaussianBeam = ScalarGaussianBeam()
    %   Returns a default scalar gaussian beam.
    %
    % newScalarGaussianBeam = ScalarGaussianBeam(centralRay,waistRadiusInX,...
    %            waistRadiusInY,distanceFromWaistInX,distanceFromWaistInY,...
    %            peakAmplitude,localXDirection,localYDirection)
    %   Returns a scalar gaussian beam with the given properties.
    
    % <<<<<<<<<<<<<<<<<<<<<<<<< Author Section >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %   Written By: Worku, Norman Girma
    %   Advisor: Prof. Herbert Gross
    %	Optical System Design and Simulation Research Group
    %   Institute of Applied Physics
    %   Friedrich-Schiller-University of Jena
    
    % <<<<<<<<<<<<<<<<<<< Change History Section >>>>>>>>>>>>>>>>>>>>>>>>>>
    % Date----------Modified By ---------Modification Detail--------Remark
    % Nov 07,2014   Worku, Norman G.     Original Version
    
    
    if nargin == 0
        % Empty constructor
        NewGaussianBeam.CentralRay = ScalarRayBundle();
        NewGaussianBeam.WaistRadiusInX = 0.05; % Default value in lens units
        NewGaussianBeam.WaistRadiusInY = 0.05;
        NewGaussianBeam.DistanceFromWaistInX = 0;
        NewGaussianBeam.DistanceFromWaistInY = 0;
        NewGaussianBeam.PeakAmplitude = 1;
        NewGaussianBeam.LocalXDirection = [1,0,0]';
        NewGaussianBeam.LocalYDirection = [0,1,0]';
    else
        if nargin == 1
            waistRadiusInX = 0.05; % Default value in lens units
            waistRadiusInY = 0.05;
            distanceFromWaistInX = 0;
            distanceFromWaistInY = 0;
            peakAmplitude = 1;
            localXDirection = [1,0,0]';
            localYDirection = [0,1,0]';
        elseif nargin == 2
            waistRadiusInY = 0.05;
            distanceFromWaistInX = 0;
            distanceFromWaistInY = 0;
            peakAmplitude = 1;
            localXDirection = [1,0,0]';
            localYDirection = [0,1,0]';
        elseif nargin == 3
            distanceFromWaistInX = 0;
            distanceFromWaistInY = 0;
            peakAmplitude = 1;
            localXDirection = [1,0,0]';
            localYDirection = [0,1,0]';
        elseif nargin == 4
            distanceFromWaistInY = 0;
            peakAmplitude = 1;
            localXDirection = [1,0,0]';
            localYDirection = [0,1,0]';
        elseif nargin == 5
            peakAmplitude = 1;
            localXDirection = [1,0,0]';
            localYDirection = [0,1,0]';
        elseif nargin == 6
            localXDirection = [1,0,0]';
            localYDirection = [0,1,0]';
        elseif nargin == 7
            localYDirection = [0,1,0]';
        end
        
        % If the inputs are arrays the NewGaussianBeam becomes object array
        % Determine the size of each inputs
        nCentralRay = size(centralRay,2);
        nWaistRadiusInX = size(waistRadiusInX,2);
        nWaistRadiusInY = size(waistRadiusInY,2);
        nDistanceFromWaistInX = size(distanceFromWaistInX,2);
        nDistanceFromWaistInY = size(distanceFromWaistInY,2);
        nPeakAmplitude = size(peakAmplitude,2);
        nLocalXDirection = size(localXDirection,2);
        nLocalYDirection = size(localYDirection,2);
        
        % The number of array to be initialized is maximum of all inputs
        nMax = max([nCentralRay,nWaistRadiusInX,nWaistRadiusInY,...
            nDistanceFromWaistInX,nPeakAmplitude,nLocalXDirection,nLocalYDirection]);
        
        % Fill the smaller inputs to have nMax size by repeating their
        % last element
        if nCentralRay < nMax
            centralRay = cat(2,centralRay,repmat(centralRay(end),...
                [1,nMax-nCentralRay]));
        end
        if nWaistRadiusInX < nMax
            waistRadiusInX = cat(2,waistRadiusInX,...
                repmat(waistRadiusInX(end), [1,nMax-nWaistRadiusInX]));
        end
        if nWaistRadiusInY < nMax
            waistRadiusInY = cat(2,waistRadiusInY,...
                repmat(waistRadiusInY(end),[1,nMax-nWaistRadiusInY]));
        end
        
        if nDistanceFromWaistInX < nMax
            distanceFromWaistInX = cat(2,distanceFromWaistInX,...
                repmat(distanceFromWaistInX(end),[1,nMax-nDistanceFromWaistInX]));
        end
        if nDistanceFromWaistInY < nMax
            distanceFromWaistInY = cat(2,distanceFromWaistInY,...
                repmat(distanceFromWaistInY(end),[1,nMax-nDistanceFromWaistInY]));
        end
        if nPeakAmplitude < nMax
            peakAmplitude = cat(2,peakAmplitude,...
                repmat(peakAmplitude(end),[1,nMax-nPeakAmplitude]));
        end
        if nLocalXDirection < nMax
            localXDirection = cat(2,localXDirection,...
                repmat(localXDirection(:,end),[1,nMax-nLocalXDirection]));
        end
        if nLocalYDirection < nMax
            localYDirection = cat(2,localYDirection,...
                repmat(localYDirection(:,end),[1,nMax-nLocalYDirection]));
        end
        
        % Preallocate object array
        NewGaussianBeam(nMax) = ScalarGaussianBeam;
        for kk = 1:nMax
            NewGaussianBeam(kk).CentralRay = centralRay(kk);
            NewGaussianBeam(kk).WaistRadiusInX = waistRadiusInX(kk);
            NewGaussianBeam(kk).WaistRadiusInY = waistRadiusInY(kk);
            NewGaussianBeam(kk).DistanceFromWaistInX = distanceFromWaistInX(kk);
            NewGaussianBeam(kk).DistanceFromWaistInY = distanceFromWaistInY(kk);
            NewGaussianBeam(kk).PeakAmplitude = peakAmplitude(kk);
            NewGaussianBeam(kk).LocalXDirection = localXDirection(:,kk);
            NewGaussianBeam(kk).LocalYDirection = localYDirection(:,kk);
        end
    end
end

