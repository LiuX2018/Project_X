function [ returnDataStruct] = Dummy(returnFlag,surfaceParameters,inputDataStruct)
    %DUMMY Dummy surface definition
    % surfaceParameters = empty
    % rayPosition and rayDirections: 3xN matrix of ray positions [x,y,z] and directions [dx,dy,dz]
    % paraxialRayParameter: 2xN matrix of paraxial heigh and angle [y,u]
    % reverseTracing: boolean indicating direction of ray trace (currently
    % supported only for paraxial tracing)
    
    % returnFlag : An integer indicating what is requested. Depending on it the
    % returnDataStruct will have different fields
    % 1: About the surface
    %   inputDataStruct:
    %       empty
    %   Output Struct:
    %       returnDataStruct.Name
    %       returnDataStruct.ImageFullFileName
    %       returnDataStruct.Description
    % 2: Surface specific 'UniqueSurfaceParameters' table field names and initial values in Surface Editor GUI
    %   inputDataStruct:
    %       empty
    %   Output Struct:
    %       returnDataStruct.UniqueParametersStructFieldNames
    %       returnDataStruct.UniqueParametersStructFieldTypes
    %       returnDataStruct.DefaultUniqueParametersStruct
    % 3: Surface specific 'Extra Data' table field names and initial values in Surface Editor GUI
    %   inputDataStruct:
    %       empty
    %   Output Struct:
    %       returnDataStruct.UniqueExtraDataFieldNames
    %       returnDataStruct.DefaultUniqueExtraData
    % 4: Return the surface sag at given xyGridPoints computed from rayPosition % Used for plotting the surface
    %   inputDataStruct:
    %       inputDataStruct.xyMeshGrid
    %   Output Struct:
    %       returnDataStruct.MainSag
    %       returnDataStruct.AlternativeSag
    % 5: Paraxial ray trace results (Ray height and angle)
    %   inputDataStruct:
    %       inputDataStruct.InputParaxialRayParameters
    %       inputDataStruct.IndexBefore
    %       inputDataStruct.IndexAfter
    %       inputDataStruct.Wavelength
    %       inputDataStruct.ReflectionFlag
    %       inputDataStruct.ReverseTracingFlag
    %   Output Struct:
    %       returnDataStruct.OutputParaxialRayParameters
    % 6: Real Ray trace results (Surface intersection points, Surface normal, Path length to the surface intersection points, Exit ray direction, Exit ray position for given
    %   inputDataStruct:
    %       inputDataStruct.InitialRayPosition
    %       inputDataStruct.InitialRayDirection
    %       inputDataStruct.IndexBefore
    %       inputDataStruct.IndexAfter
    %       inputDataStruct.Wavelength
    %       inputDataStruct.ReferenceWavelength
    %       inputDataStruct.ReflectionFlag
    %       inputDataStruct.ReverseTracingFlag
    %   Output Struct:
    %       returnDataStruct.GeometricalPathLength
    %       returnDataStruct.AdditionalPathLength
    %       returnDataStruct.LocalRayIntersectionPoint
    %       returnDataStruct.LocalSurfaceNormal
    %       returnDataStruct.LocalExitRayPosition
    %       returnDataStruct.LocalExitRayDirection
    %       returnDataStruct.TotalInternalReflectionFlag
    %       returnDataStruct.NoIntersectionPointFlag
    
    % <<<<<<<<<<<<<<<<<<<<<<<<< Author Section >>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %   Written By: Worku, Norman Girma
    %   Advisor: Prof. Herbert Gross
    %	Optical System Design and Simulation Research Group
    %   Institute of Applied Physics
    %   Friedrich-Schiller-University of Jena
    
    % <<<<<<<<<<<<<<<<<<< Change History Section >>>>>>>>>>>>>>>>>>>>>>>>>>
    % Date----------Modified By ---------Modification Detail--------Remark
    % Jun 17,2015   Worku, Norman G.     Original Version
    % Jul 10,2015   Worku, Norman G.     input and output are made struct    
    %% Default input vaalues
    if nargin == 0
        disp('Error: The function Dummy() needs atleat the return type.');
        returnDataStruct = struct;
        return;
    elseif nargin == 1 || nargin == 2
        if returnFlag == 1 || returnFlag == 2 || returnFlag == 3 || returnFlag == 4
            inputDataStruct = struct();
        else
            disp('Error: Missing input argument for Dummy().');
            returnDataStruct = struct();
            return;
        end
    elseif nargin == 3
        % This is fine
    else
        
    end
    switch returnFlag
        case 1 % About the surface
            surfName = {'Dummy','DUMY'}; % display name
            % look for image description in the current folder and return
            % full address
            [pathstr,name,ext] = fileparts(mfilename('fullpath'));
            imageFullFileName = {[pathstr,'\Surface.jpg']};  % Image file name
            description = {['Dummy Surface: Is is just surface with no effects ',...
             ' on the optical system.']};  % Text description
            
            returnDataStruct = struct();
            returnDataStruct.Name = surfName;
            returnDataStruct.ImageFullFileName = imageFullFileName;
            returnDataStruct.Description =  description;
        case 2 % Surface specific 'UniqueSurfaceParameters'
            uniqueParametersStructFieldNames = {'Unused'};
            uniqueParametersStructFieldTypes = {{'numeric'}};
            defaultUniqueParametersStruct = struct();
            defaultUniqueParametersStruct.Unused = 0;
                        
            returnDataStruct = struct();
            returnDataStruct.UniqueParametersStructFieldNames = uniqueParametersStructFieldNames;
            returnDataStruct.UniqueParametersStructFieldTypes = uniqueParametersStructFieldTypes;
            returnDataStruct.DefaultUniqueParametersStruct = defaultUniqueParametersStruct;
        case 3 % Surface specific 'Extra Data' table
            uniqueExtraDataFieldNames = {'Unused'};
            defaultUniqueExtraData = {[0]};
            
            returnDataStruct = struct();
            returnDataStruct.UniqueExtraDataFieldNames = uniqueExtraDataFieldNames;
            returnDataStruct.DefaultUniqueExtraData = defaultUniqueExtraData;
        case 4 % Surface sag at given xyGridPoints
            disp('Error: Sag of Dummy surface can not be computed.');
            returnDataStruct = struct();
        case 5 % Paraxial ray trace results
            disp('Error: Paraxial ray trace results of Dummy surface can not be computed.');
            returnDataStruct = struct();
        case 6 % Real Ray trace results
            disp('Error: Real ray trace results of Dummy surface can not be computed.');
            returnDataStruct = struct();
        otherwise
    end
end

