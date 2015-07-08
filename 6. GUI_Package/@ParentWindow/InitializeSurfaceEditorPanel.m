function InitializeSurfaceEditorPanel( parentWindow )
    %INITIALIZESURFACEEDITORPANEL Define and initialized the uicontrols of the
    % Surface Editor Panel
    % Member of ParentWindow class
    aodHandles = parentWindow.ParentHandles;
    fontSize = aodHandles.FontSize;
    fontName = aodHandles.FontName;
    
    %% Divide the area in to surface list panel, and surf detail
    % ( surf figure, surf description and surface parameters) panel
    aodHandles.panelSurfaceList = uipanel(...
        'Parent',aodHandles.panelSurfaceEditorMain,...
        'FontSize',fontSize,'FontName', fontName,...
        'Title','Surface List',...
        'units','normalized',...
        'Position',[0.0,0.05,0.5,0.95]);
    
    aodHandles.panelSurfaceDetail = uipanel(...
        'Parent',aodHandles.panelSurfaceEditorMain,...
        'FontSize',fontSize,'FontName', fontName,...
        'units','normalized',...
        'Position',[0.5,0.05,0.5,0.95]);
    
    aodHandles.surfaceParametersTabGroup = uitabgroup(...
        'Parent', aodHandles.panelSurfaceDetail, ...
        'Units', 'Normalized', ...
        'Position', [0,0.0,1.0,1.0]);
    aodHandles.surfBasicDataTab = ...
        uitab(aodHandles.surfaceParametersTabGroup, 'title','Standard Data');
    aodHandles.surfApertureDataTab = ...
        uitab(aodHandles.surfaceParametersTabGroup, 'title','Aperture Data');
    aodHandles.surfTiltDecenterDataTab = ...
        uitab(aodHandles.surfaceParametersTabGroup, 'title','Tilt Decenter Data');
    
    
    % Initialize the ui table for surfacelist
    aodHandles.tblSurfaceList = uitable('Parent',aodHandles.panelSurfaceList,...
        'FontSize',fontSize,'FontName', fontName,...
        'units','normalized','Position',[0 0 1 0.93]);
    
    % Command buttons for adding and removing surfaces
    aodHandles.btnInsertSurface = uicontrol( ...
        'Parent',aodHandles.panelSurfaceList,...
        'Style', 'pushbutton', ...
        'FontSize',fontSize,'FontName', fontName,...
        'units','normalized',...
        'Position',[0.01,0.94,0.2,0.05],...
        'String','Insert',...
        'Callback',{@btnInsertSurface_Callback,parentWindow});
    aodHandles.btnRemoveSurface = uicontrol( ...
        'Parent',aodHandles.panelSurfaceList,...
        'Style', 'pushbutton', ...
        'FontSize',fontSize,'FontName', fontName,...
        'units','normalized',...
        'Position',[0.21,0.94,0.2,0.05],...
        'String','Remove',...
        'Callback',{@btnRemoveSurface_Callback,parentWindow});

    updatedSystem = aodHandles.OpticalSystem;
    updatedSystem.SurfaceArray = updateSurfaceCoordinateTransformationMatrices(aodHandles.OpticalSystem.SurfaceArray);
    
    if IsSurfaceBased(updatedSystem)
        nSurface = getNumberOfSurfaces(updatedSystem);
        stopSurfaceIndex = getStopSurfaceIndex(updatedSystem);
    else
        nSurface = 3;
        stopSurfaceIndex = 2;
    end
    
    
    aodHandles.lblStopSurfaceIndex = uicontrol( ...
        'Parent',aodHandles.panelSurfaceList,...
        'Tag', 'lblStopSurfaceIndex', ...
        'Style', 'text', ...
        'HorizontalAlignment','left',...
        'FontSize',fontSize,'FontName', 'FixedWidth',...
        'String', 'Stop Surface ',...
        'units','normalized',...
        'Position',[0.46,0.935,0.30,0.04]);
    
    aodHandles.popStopSurfaceIndex = uicontrol( ...
        'Parent',aodHandles.panelSurfaceList,...
        'Tag', 'popStopSurfaceIndex', ...
        'FontSize',fontSize,'FontName', 'FixedWidth',...
        'Style', 'popupmenu', ...
        'BackgroundColor', [1 1 1], ...
        'String', [num2cell(1:nSurface)],...
        'Value',stopSurfaceIndex,...
        'units','normalized',...
        'Position',[0.70,0.93,0.15,0.05]);
    
    % Set  celledit and cellselection callbacks
    set(aodHandles.tblSurfaceList,...
        'CellSelectionCallback',{@tblSurfaceList_CellSelectionCallback,parentWindow},...
        'CellEditCallback',{@tblSurfaceList_CellEditCallback,parentWindow});
    
    set(aodHandles.popStopSurfaceIndex,...
        'Callback',{@popStopSurfaceIndex_Callback,parentWindow});
    
    % Initialize the ui table for surface parameters
    aodHandles.tblSurfaceBasicParameters = uitable('Parent',aodHandles.surfBasicDataTab,...
        'FontSize',fontSize,'FontName', fontName,...
        'units','normalized','Position',[0 0 1 1]);
    aodHandles.tblSurfaceApertureParameters = uitable('Parent',aodHandles.surfApertureDataTab,...
        'FontSize',fontSize,'FontName', fontName,...
        'units','normalized','Position',[0 0 1 1]);
    aodHandles.tblSurfaceTiltDecenterParameters = uitable('Parent',aodHandles.surfTiltDecenterDataTab,...
        'FontSize',fontSize,'FontName', fontName,...
        'units','normalized','Position',[0 0 1 1]);
    
    %     Set  celledit and cellselection callbacks
    set(aodHandles.tblSurfaceBasicParameters,...
        'CellSelectionCallback',{@tblSurfaceBasicParameters_CellSelectionCallback,parentWindow},...
        'CellEditCallback',{@tblSurfaceBasicParameters_CellEditCallback,parentWindow});
    set(aodHandles.tblSurfaceApertureParameters,...
        'CellSelectionCallback',{@tblSurfaceApertureParameters_CellSelectionCallback,parentWindow},...
        'CellEditCallback',{@tblSurfaceApertureParameters_CellEditCallback,parentWindow});
    set(aodHandles.tblSurfaceTiltDecenterParameters,...
        'CellSelectionCallback',{@tblSurfaceTiltDecenterParameters_CellSelectionCallback,parentWindow},...
        'CellEditCallback',{@tblSurfaceTiltDecenterParameters_CellEditCallback,parentWindow});
    
    supportedSurfaces = GetSupportedSurfaces();
    columnName1 =   {'Surface','Type', 'Name/Note'};
    columnWidth1 = {70,120, 100};
    columnEditable1 =  [false,true ,true];
    columnFormat1 =  {'char',{supportedSurfaces{:}}, 'char'};
    initialTable1 = {'OBJ','OBJECT','Object';...
        'STOP','Standard',  '';...
        'IMG','IMAGE', 'Image'};
    set(aodHandles.tblSurfaceList, ...
        'ColumnFormat',columnFormat1,...
        'Data', initialTable1,'ColumnEditable', columnEditable1,...
        'ColumnName', columnName1,'ColumnWidth',columnWidth1);
    
    % reset surface parameter tables
    columnName2 =   {'Parameters', 'Value', 'Solve'};
    columnWidth2 = {150, 100 70};
    columnEditable2 =  [false true false];
    columnFormat2 =  {'char', 'char','char'};
    initialTable2 = {'Param 1','0',' '};
    set(aodHandles.tblSurfaceBasicParameters, ...
        'ColumnFormat',columnFormat2,...
        'Data', initialTable2,'ColumnEditable', columnEditable2,...
        'ColumnName', columnName2,'ColumnWidth',columnWidth2);
    set(aodHandles.tblSurfaceApertureParameters, ...
        'ColumnFormat',columnFormat2,...
        'Data', initialTable2,'ColumnEditable', columnEditable2,...
        'ColumnName', columnName2,'ColumnWidth',columnWidth2);
    set(aodHandles.tblSurfaceTiltDecenterParameters, ...
        'ColumnFormat',columnFormat2,...
        'Data', initialTable2,'ColumnEditable', columnEditable2,...
        'ColumnName', columnName2,'ColumnWidth',columnWidth2);
    
    % Give all tables initial data
    parentWindow.ParentHandles = aodHandles;
    updateSystemConfigurationWindow( parentWindow );
    updateQuickLayoutPanel(parentWindow,1);
end

% Local functions
function btnInsertSurface_Callback(~,~,parentWindow)
    if ~checkTheCurrentSystemDefinitionType(parentWindow.ParentHandles)
        return
    end
    global CURRENT_SELECTED_SURFACE
    global CAN_ADD_SURFACE
    if CAN_ADD_SURFACE
        insertPosition = CURRENT_SELECTED_SURFACE;
        InsertNewSurface(parentWindow,'Standard','Standard',insertPosition);
        aodHandles = parentWindow.ParentHandles;
    end
    updateQuickLayoutPanel(parentWindow,CURRENT_SELECTED_SURFACE);
end



function btnRemoveSurface_Callback(~,~,parentWindow)
    if ~checkTheCurrentSystemDefinitionType(parentWindow.ParentHandles)
        return
    end
    global CAN_REMOVE_SURFACE
    global CURRENT_SELECTED_SURFACE
    
    aodHandles = parentWindow.ParentHandles;
    if CAN_REMOVE_SURFACE
        % Confirm action
        % Construct a questdlg with three options
        choice = questdlg('Are you sure to delete the surface?', ...
            'Confirm Deletion', ...
            'Yes','No','Yes');
        % Handle response
        switch choice
            case 'Yes'
                % Delete the surface
                removePosition = CURRENT_SELECTED_SURFACE;
                RemoveSurface(parentWindow,removePosition);
                aodHandles = parentWindow.ParentHandles;
            case 'No'
                % Mark the delete box again
        end
    else
        % Mark the delete box again
    end
    parentWindow.ParentHandles = aodHandles;
    updateQuickLayoutPanel(parentWindow,CURRENT_SELECTED_SURFACE);
end

function popStopSurfaceIndex_Callback(hObject, eventdata,parentWindow)
    aodHandles = parentWindow.ParentHandles;
    currentOpticalSystem = aodHandles.OpticalSystem;
    prevStopIndex = getStopSurfaceIndex(currentOpticalSystem );
    newStopIndex = get(hObject,'Value');
    if newStopIndex == 1 || newStopIndex == getNumberOfSurfaces(currentOpticalSystem) % Object and image can not be Stop
        newStopIndex = prevStopIndex;
    end
    if newStopIndex ~= prevStopIndex
        aodHandles.OpticalSystem.SurfaceArray(prevStopIndex).Stop = 0;
        aodHandles.OpticalSystem.SurfaceArray(newStopIndex).Stop = 1;
    end
    parentWindow.ParentHandles = aodHandles;
    updateSurfaceOrComponentEditorPanel( parentWindow , newStopIndex);
end

% Cell select and % CellEdit Callback
% --- Executes when selected cell(s) is changed in aodHandles.tblSurfaceLis.
function tblSurfaceList_CellSelectionCallback(~, eventdata,parentWindow)
    global CURRENT_SELECTED_SURFACE
    global CAN_ADD_SURFACE
    global CAN_REMOVE_SURFACE
    aodHandles = parentWindow.ParentHandles;
    
    selectedCell = eventdata.Indices;
    if isempty(selectedCell)
        return
    end
    CURRENT_SELECTED_SURFACE = selectedCell(1);
    
    tblData = get(aodHandles.tblSurfaceList,'data');
    sizeTblData = size(tblData);
    
    if CURRENT_SELECTED_SURFACE == 1 % object surf
        CAN_ADD_SURFACE = 0;
        CAN_REMOVE_SURFACE = 0;
        % make the 2nd column uneditable
        columnEditable1 = [false,false,true];
        set(aodHandles.tblSurfaceList,'ColumnEditable', columnEditable1);
    elseif CURRENT_SELECTED_SURFACE == sizeTblData(1)% image surf
        CAN_ADD_SURFACE = 1;
        CAN_REMOVE_SURFACE = 0;
        % make the 2nd column uneditable
        columnEditable1 = [false,false,true];
        set(aodHandles.tblSurfaceList,'ColumnEditable', columnEditable1);
    elseif sizeTblData(1) == 3 % only 3 surfaces left
        CAN_ADD_SURFACE = 1;
        CAN_REMOVE_SURFACE = 0;
        % make the 2nd column editable
        columnEditable1 = [false,true,true];
        set(aodHandles.tblSurfaceList,'ColumnEditable', columnEditable1);
    else
        CAN_ADD_SURFACE = 1;
        CAN_REMOVE_SURFACE = 1;
        % make the 2nd column editable
        columnEditable1 = [false,true,true];
        set(aodHandles.tblSurfaceList,'ColumnEditable', columnEditable1);
    end
    
    % Show the surface parameters in the surface detail window
    parentWindow.ParentHandles = aodHandles;
    updateQuickLayoutPanel(parentWindow,CURRENT_SELECTED_SURFACE)
    updateSurfaceOrComponentEditorPanel( parentWindow , CURRENT_SELECTED_SURFACE);
end
% --- Executes when entered data in editable cell(s) in aodHandles.tblSurfaceLis.
function tblSurfaceList_CellEditCallback(~, eventdata,parentWindow)
    % hObject    handle to aodHandles.tblSurfaceList (see GCBO)
    % eventdata  structure with the following fields (see UITABLE)
    %	Indices: row and column indices of the cell(s) edited
    %	PreviousData: previous data for the cell(s) edited
    %	EditData: string(s) entered by the user
    %	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
    %	Error: error string when failed to convert EditData to appropriate value for Data
    % parentWindow: object with structure with aodHandles and user data (see GUIDATA)
    
    global CURRENT_SELECTED_SURFACE
    aodHandles = parentWindow.ParentHandles;
    editedCellIndex = eventdata.Indices;
    if ~isempty(editedCellIndex)
        editedRow = editedCellIndex(1,1);
        editedCol = editedCellIndex(1,2);
    else
        return;
    end
    CURRENT_SELECTED_SURFACE = editedRow;
    if editedCol== 2 % surface type changed
        if strcmpi(eventdata.PreviousData,'OBJECT')||strcmpi(eventdata.PreviousData,'IMAGE')
            % for object or image surf surf change the surf type back to OBJECT or
            % IMAGE
            tblData1 = get(aodHandles.tblSurfaceList,'data');
            sizeTblData1 = size(tblData1);
            nSurface = sizeTblData1(1);
            tblData1{1,2} = 'OBJECT';
            tblData1{nSurface,2} = 'IMAGE';
            set(aodHandles.tblSurfaceList, 'Data', tblData1);
        else
            if ~checkTheCurrentSystemDefinitionType(aodHandles)
                return
            end
            % reset the surface type in the surface detail window
            selectedSurfaceType = eventdata.NewData;
            newSurface = Surface(selectedSurfaceType);
            % Add the new surfonet to the temporary surfaceArray
            aodHandles.OpticalSystem.SurfaceArray(editedRow) = newSurface;
        end
        
    elseif editedCol== 3 % surface comment changed
        aodHandles.OpticalSystem.SurfaceArray(editedRow).Comment = eventdata.NewData;
    end
    parentWindow.ParentHandles = aodHandles;
    updateQuickLayoutPanel(parentWindow,CURRENT_SELECTED_SURFACE);
    updateSurfaceOrComponentEditorPanel( parentWindow , CURRENT_SELECTED_SURFACE)
    
end



function tblSurfaceBasicParameters_CellSelectionCallback(~, eventdata,parentWindow)
    % hObject    handle to aodHandles.tblSurfaceList (see GCBO)
    % eventdata  structure with the following fields (see UITABLE)
    %	Indices: row and column indices of the cell(s) edited
    %	PreviousData: previous data for the cell(s) edited
    %	EditData: string(s) entered by the user
    %	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
    %	Error: error string when failed to convert EditData to appropriate value for Data
    % parentWindow: object with structure with aodHandles and user data (see GUIDATA)
    
    global CURRENT_SELECTED_SURFACE
    aodHandles = parentWindow.ParentHandles;
    selectedSurface =  aodHandles.OpticalSystem.SurfaceArray(CURRENT_SELECTED_SURFACE);
    selectedCellIndex = eventdata.Indices;
    if ~isempty(selectedCellIndex)
        selectedRow = selectedCellIndex(1,1);
        selectedCol = selectedCellIndex(1,2);
    else
        return;
    end
    
    if selectedCol ~= 2
        return;
    end
    
    [ paramNames,paramTypes,paramValues,paramValuesDisp] = ...
        getSurfaceParameters(selectedSurface,'Basic',selectedRow);
    
    myType = paramTypes{1};
    myName = paramNames{1};
    if  iscell(myType) && length(myType)>1
        % type is choice of popmenu
        nChoice = length(myType);
        choice = menu(myName,myType(1:nChoice));
        if choice == 0
            choice = 1;
        end
        newParam = myType{choice};
        newParamDisp = newParam;
        % Update the surface parameter and surface editor
        if selectedRow <=3
            selectedSurface.(myName) = newParam;
        else
            selectedSurface.UniqueParameters.(myName) = newParam;
        end
    else
        if strcmpi('logical',myType)
            % type is choice of popmenu true or false
            trueFalse = {'1','0'};
            choice = menu(myName,'False','True');
            if choice == 0
                choice = 1;
            end
            newParam = trueFalse{choice};
            newParamDisp = newParam;
            
            % Update the surface parameter and surface editor
            if selectedRow <=3
                selectedSurface.(myName) = newParam;
            else
                selectedSurface.UniqueParameters.(myName) = newParam;
            end
        elseif strcmpi('numeric',myType) || strcmpi('char',myType)||...
                strcmpi('Glass',myType) || strcmpi('Coating',myType)
            
        else
            disp(['Error: Unknown parameter type: ',myType]);
            return;
        end
    end
    aodHandles.OpticalSystem.SurfaceArray(CURRENT_SELECTED_SURFACE) = selectedSurface;
    parentWindow.ParentHandles = aodHandles;
    updateSurfaceOrComponentEditorPanel( parentWindow , CURRENT_SELECTED_SURFACE)
    
end

function tblSurfaceApertureParameters_CellSelectionCallback(~, eventdata,parentWindow)
    % hObject    handle to aodHandles.tblSurfaceList (see GCBO)
    % eventdata  structure with the following fields (see UITABLE)
    %	Indices: row and column indices of the cell(s) edited
    %	PreviousData: previous data for the cell(s) edited
    %	EditData: string(s) entered by the user
    %	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
    %	Error: error string when failed to convert EditData to appropriate value for Data
    % parentWindow: object with structure with aodHandles and user data (see GUIDATA)
    
    global CURRENT_SELECTED_SURFACE
    aodHandles = parentWindow.ParentHandles;
    selectedSurface =  aodHandles.OpticalSystem.SurfaceArray(CURRENT_SELECTED_SURFACE);
    selectedCellIndex = eventdata.Indices;
    if ~isempty(selectedCellIndex)
        selectedRow = selectedCellIndex(1,1);
        selectedCol = selectedCellIndex(1,2);
    else
        return;
    end
    if selectedCol ~= 2
        return;
    end
    [ paramNames,paramTypes,paramValues,paramValuesDisp] = ...
        getSurfaceParameters(selectedSurface,'Aperture',selectedRow);
    
    myType = paramTypes{1};
    myName = paramNames{1};
    if  iscell(myType)&& length(myType)>1
        % type is choice of popmenu
        nChoice = length(myType);
        choice = menu(myName,myType(1:nChoice));
        if choice == 0
            choice = 1;
        end
        newParam = myType{choice};
        newParamDisp = newParam;
        % Update the surface parameter and surface editor
        if selectedRow == 1
            % A new type with defualt parameter
            selectedSurface.Aperture = Aperture(newParam);
        elseif selectedRow <=7
            selectedSurface.Aperture.(myName) = newParam;
        else
            selectedSurface.Aperture.UniqueParameters.(myName) = newParam;
        end
    else
        if strcmpi('logical',myType)
            % type is choice of popmenu true or false
            trueFalse = {'0','1'};
            choice = menu(myName,'False','True');
            if choice == 0
                choice = 1;
            end
            newParam = trueFalse{choice};
            
            newParamDisp = newParam;
            
            % Update the surface parameter and surface editor
            if selectedRow == 1
                % A new type with defualt parameter
                selectedSurface.Aperture = Aperture(newParam);
            elseif selectedRow <=7
                selectedSurface.Aperture.(myName) = newParam;
            else
                selectedSurface.Aperture.UniqueParameters.(myName) = newParam;
            end
        elseif strcmpi('numeric',myType) || strcmpi('char',myType)
            
        else
            disp(['Error: Unknown parameter type: ',myType]);
            return;
        end
    end
    aodHandles.OpticalSystem.SurfaceArray(CURRENT_SELECTED_SURFACE) = selectedSurface;
    parentWindow.ParentHandles = aodHandles;
    updateSurfaceOrComponentEditorPanel( parentWindow , CURRENT_SELECTED_SURFACE)
end

function tblSurfaceTiltDecenterParameters_CellSelectionCallback(~, eventdata,parentWindow)
    % hObject    handle to aodHandles.tblSurfaceTiltDecenterParameters (see GCBO)
    % eventdata  structure with the following fields (see UITABLE)
    %	Indices: row and column indices of the cell(s) edited
    %	PreviousData: previous data for the cell(s) edited
    %	EditData: string(s) entered by the user
    %	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
    %	Error: error string when failed to convert EditData to appropriate value for Data
    % parentWindow: object with structure with aodHandles and user data (see GUIDATA)
    
    global CURRENT_SELECTED_SURFACE
    aodHandles = parentWindow.ParentHandles;
    selectedSurface =  aodHandles.OpticalSystem.SurfaceArray(CURRENT_SELECTED_SURFACE);
    selectedCellIndex = eventdata.Indices;
    if ~isempty(selectedCellIndex)
        selectedRow = selectedCellIndex(1,1);
        selectedCol = selectedCellIndex(1,2);
    else
        return;
    end
    if selectedCol ~= 2
        return;
    end
    [ paramNames,paramTypes,paramValues,paramValuesDisp] = ...
        getSurfaceParameters(selectedSurface,'TiltDecenter',selectedRow);
    
    myType = paramTypes{1};
    myName = paramNames{1};
    if  iscell(myType)&& length(myType)>1
        % type is choice of popmenu
        nChoice = length(myType);
        choice = menu(myName,myType(1:nChoice));
        if choice == 0
            choice = 1;
        end
        newParam = myType{choice};
        newParamDisp = newParam;
        % Update the surface parameter and surface editor
        selectedSurface.(myName) = newParam;
    else
        if strcmpi('logical',myType)
            % type is choice of popmenu true or false
            trueFalse = {'0','1'};
            choice = menu(myName,'False','True');
            if choice == 0
                choice = 1;
            end
            newParam = trueFalse{choice};
            
            newParamDisp = newParam;
            
            % Update the surface parameter and surface editor
            selectedSurface.(myName) = newParam;
        elseif strcmpi('numeric',myType) || strcmpi('char',myType)
            
        else
            disp(['Error: Unknown parameter type: ',myType]);
            return;
        end
    end
    aodHandles.OpticalSystem.SurfaceArray(CURRENT_SELECTED_SURFACE) = selectedSurface;
    parentWindow.ParentHandles = aodHandles;
    updateSurfaceOrComponentEditorPanel( parentWindow , CURRENT_SELECTED_SURFACE)
    
end


function tblSurfaceBasicParameters_CellEditCallback(~, eventdata,parentWindow)
    % hObject    handle to aodHandles.tblSurfaceList (see GCBO)
    % eventdata  structure with the following fields (see UITABLE)
    %	Indices: row and column indices of the cell(s) edited
    %	PreviousData: previous data for the cell(s) edited
    %	EditData: string(s) entered by the user
    %	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
    %	Error: error string when failed to convert EditData to appropriate value for Data
    % parentWindow: object with structure with aodHandles and user data (see GUIDATA)
    
    global CURRENT_SELECTED_SURFACE
    aodHandles = parentWindow.ParentHandles;
    selectedSurface =  aodHandles.OpticalSystem.SurfaceArray(CURRENT_SELECTED_SURFACE);
    if CURRENT_SELECTED_SURFACE > 1
        previousSurface =  aodHandles.OpticalSystem.SurfaceArray(CURRENT_SELECTED_SURFACE-1);
    else
        previousSurface =  selectedSurface;
    end
    selectedCellIndex = eventdata.Indices;
    if isempty(selectedCellIndex)
        return;
    else
        editedRow = selectedCellIndex(1,1);
        editedCol = selectedCellIndex(1,2);
        
    end
    if editedCol == 2 % surface  value edited
        switch editedRow
            case 1
                newParam = str2num(eventdata.EditData);
                if isempty(newParam)
                    newParam = 0;
                end
                selectedSurface.Thickness = newParam;
            case 2
                glassName = (eventdata.EditData);
                newParam = Glass(glassName);
                selectedSurface.Glass = newParam;
            case 3
                coatingName = (eventdata.EditData);
                newParam = Coating(coatingName);
                selectedSurface.Coating = newParam;
            otherwise
                [ paramNames,paramTypes,paramValues,paramValuesDisp] = ...
                    getSurfaceParameters(selectedSurface,'Basic',editedRow);
                
                myType = paramTypes{1};
                myName = paramNames{1};
                
                if  (iscell(myType) && length(myType)>1)||(strcmpi('logical',myType))
                    % type is choice of popmenu and so already saved with popup menu
                else
                    if strcmpi('numeric',myType)
                        newParam = str2num(eventdata.EditData);
                        if isempty(newParam)
                            newParam = 0;
                        end
                        % Update the surface parameter and surface editor
                        selectedSurface.UniqueParameters.(myName) = newParam;
                    elseif strcmpi('char',myType)
                        newParam = (eventdata.EditData);
                        % Update the surface parameter and surface editor
                        selectedSurface.UniqueParameters.(myName) = newParam;
                    elseif strcmpi('Glass',myType)
                        glassName = (eventdata.EditData);
                        newParam = Glass(glassName);
                        % Update the surface parameter and surface editor
                        selectedSurface.UniqueParameters.(myName) = newParam;
                    elseif strcmpi('Coating',myType)
                        coatingName = (eventdata.EditData);
                        newParam = Glass(coatingName);
                        % Update the surface parameter and surface editor
                        selectedSurface.UniqueParameters.(myName) = newParam;
                    else
                        disp(['Error: Unknown parameter type: ',myType]);
                        return;
                    end
                end
                
        end
    else
    end
    aodHandles.OpticalSystem.SurfaceArray(CURRENT_SELECTED_SURFACE) = selectedSurface;
    parentWindow.ParentHandles = aodHandles;
    updateQuickLayoutPanel(parentWindow,CURRENT_SELECTED_SURFACE);
    updateSurfaceOrComponentEditorPanel( parentWindow );
end

function tblSurfaceApertureParameters_CellEditCallback(~, eventdata,parentWindow)
    % hObject    handle to aodHandles.tblSurfaceList (see GCBO)
    % eventdata  structure with the following fields (see UITABLE)
    %	Indices: row and column indices of the cell(s) edited
    %	PreviousData: previous data for the cell(s) edited
    %	EditData: string(s) entered by the user
    %	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
    %	Error: error string when failed to convert EditData to appropriate value for Data
    % parentWindow: object with structure with aodHandles and user data (see GUIDATA)
    
    global CURRENT_SELECTED_SURFACE
    aodHandles = parentWindow.ParentHandles;
    selectedSurface =  aodHandles.OpticalSystem.SurfaceArray(CURRENT_SELECTED_SURFACE);
    
    selectedCellIndex = eventdata.Indices;
    if isempty(selectedCellIndex) || isempty(eventdata.NewData)
        return;
    else
        editedRow = selectedCellIndex(1,1);
        editedCol = selectedCellIndex(1,2);
        
    end
    if editedCol == 2 % surface  value edited
        switch editedRow
            case 1
                % already saved by poup menu
            case 2
                newParam = str2num(eventdata.EditData);
                if isempty(newParam)
                    newParam = 0;
                end
                selectedSurface.Aperture.Decenter(1) = newParam;
            case 3
                newParam = str2num(eventdata.EditData);
                if isempty(newParam)
                    newParam = 0;
                end
                selectedSurface.Aperture.Decenter(2) = newParam;
            case 4
                newParam = str2num(eventdata.EditData);
                if isempty(newParam)
                    newParam = 0;
                end
                selectedSurface.Aperture.Rotation= newParam;
            case 5
                % already saved
            case 6
                % already saved
            case 7
                newParam = str2num(eventdata.EditData);
                if isempty(newParam)
                    newParam = 0;
                end
                selectedSurface.Aperture.AdditionalEdge = newParam;
            otherwise
                [ paramNames,paramTypes,paramValues,paramValuesDisp] = ...
                    getSurfaceParameters(selectedSurface,'Aperture',editedRow);
                
                myType = paramTypes{1};
                myName = paramNames{1};
                
                if  (iscell(myType) && length(myType)>1)||(strcmpi('logical',myType))
                    % type is choice of popmenu and so already saved with popup menu
                else
                    if strcmpi('numeric',myType)
                        newParam = str2num(eventdata.EditData);
                        if isempty(newParam)
                            newParam = 0;
                        end
                        % Update the surface parameter and surface editor
                        selectedSurface.Aperture.UniqueParameters.(myName) = newParam;
                    elseif strcmpi('char',myType)
                        newParam = (eventdata.EditData);
                        % Update the surface parameter and surface editor
                        selectedSurface.Aperture.UniqueParameters.(myName) = newParam;
                    elseif strcmpi('Glass',myType)
                        glassName = (eventdata.EditData);
                        newParam = Glass(glassName);
                        % Update the surface parameter and surface editor
                        selectedSurface.Aperture.UniqueParameters.(myName) = newParam;
                    elseif strcmpi('Coating',myType)
                        coatingName = (eventdata.EditData);
                        newParam = Glass(coatingName);
                        % Update the surface parameter and surface editor
                        selectedSurface.Aperture.UniqueParameters.(myName) = newParam;
                    else
                        disp(['Error: Unknown parameter type: ',myType]);
                        return;
                    end
                end
        end
    else
        
    end
    %     displaySurfaceDetail(selectedSurface,aodHandles);
    aodHandles.OpticalSystem.SurfaceArray(CURRENT_SELECTED_SURFACE) = selectedSurface;
    parentWindow.ParentHandles = aodHandles;
    updateSurfaceOrComponentEditorPanel( parentWindow , CURRENT_SELECTED_SURFACE)
    
end

function tblSurfaceTiltDecenterParameters_CellEditCallback(~, eventdata,parentWindow)
    % hObject    handle to aodHandles.tblSurfaceTiltDecenterParameters (see GCBO)
    % eventdata  structure with the following fields (see UITABLE)
    %	Indices: row and column indices of the cell(s) edited
    %	PreviousData: previous data for the cell(s) edited
    %	EditData: string(s) entered by the user
    %	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
    %	Error: error string when failed to convert EditData to appropriate value for Data
    % parentWindow: object with structure with aodHandles and user data (see GUIDATA)
    
    
    global CURRENT_SELECTED_SURFACE
    aodHandles = parentWindow.ParentHandles;
    selectedSurface =  aodHandles.OpticalSystem.SurfaceArray(CURRENT_SELECTED_SURFACE);
    
    selectedCellIndex = eventdata.Indices;
    if isempty(selectedCellIndex) || isempty(eventdata.NewData)
        return;
    else
        editedRow = selectedCellIndex(1,1);
        editedCol = selectedCellIndex(1,2);
        
    end
    if editedCol == 2 % surface  value edited
        switch editedRow
            case 1
                orderString = eventdata.EditData;
                if length(orderString) ~= 12
                    selectedSurface.TiltDecenterOrder = {'Dx','Dy','Dz','Tx','Ty','Tz'};
                else
                    selectedSurface.TiltDecenterOrder = {orderString(1:2),orderString(3:4),...
                        orderString(5:6),orderString(7:8),orderString(9:10),orderString(11:12)};
                end
                
            case 2
                newParam = str2num(eventdata.EditData);
                if isempty(newParam)
                    newParam = 0;
                end
                selectedSurface.Tilt(1) = newParam;
            case 3
                newParam = str2num(eventdata.EditData);
                if isempty(newParam)
                    newParam = 0;
                end
                selectedSurface.Tilt(2) = newParam;
            case 4
                newParam = str2num(eventdata.EditData);
                if isempty(newParam)
                    newParam = 0;
                end
                selectedSurface.Tilt(3) = newParam;
            case 5
                newParam = str2num(eventdata.EditData);
                if isempty(newParam)
                    newParam = 0;
                end
                selectedSurface.Decenter(1) = newParam;
            case 6
                newParam = str2num(eventdata.EditData);
                if isempty(newParam)
                    newParam = 0;
                end
                selectedSurface.Decenter(2) = newParam;
            case 7
                % already saved by pop up menu
            otherwise
        end
    else
    end
    aodHandles.OpticalSystem.SurfaceArray(CURRENT_SELECTED_SURFACE) = selectedSurface;
    parentWindow.ParentHandles = aodHandles;
    updateSurfaceOrComponentEditorPanel( parentWindow , CURRENT_SELECTED_SURFACE)
    
end


function InsertNewSurface(parentWindow,surfaceTypeDisp,surfaceType,insertPosition)
    %update surface list table
    aodHandles = parentWindow.ParentHandles;
    nSurface = getNumberOfSurfaces(aodHandles.OpticalSystem);
    % Update the surface array
    for kk = nSurface:-1:insertPosition
        aodHandles.OpticalSystem.SurfaceArray(kk+1) = aodHandles.OpticalSystem.SurfaceArray(kk);
    end
    aodHandles.OpticalSystem.SurfaceArray(insertPosition) = Surface(surfaceType);
    
    set(aodHandles.popStopSurfaceIndex,'String',[num2cell(1:nSurface+1)],...
        'Value',getStopSurfaceIndex(aodHandles.OpticalSystem));
    % If possible add here a code to select the first cell of newly added row
    % automatically
    parentWindow.ParentHandles = aodHandles;
    updateSurfaceOrComponentEditorPanel( parentWindow );
end

function RemoveSurface(parentWindow,removePosition)
   global CAN_REMOVE_SURFACE
    % check if it is stop surface
    if getStopSurfaceIndex(parentWindow.ParentHandles.OpticalSystem) == removePosition
        stopSurfaceRemoved = 1;
    else
        stopSurfaceRemoved = 0;
    end
    aodHandles = parentWindow.ParentHandles;
    
    % Update the surface array
    aodHandles.OpticalSystem.SurfaceArray = aodHandles.OpticalSystem.SurfaceArray([1:removePosition-1,removePosition+1:end]);
    if stopSurfaceRemoved
        if aodHandles.OpticalSystem.SurfaceArray(removePosition).ImageSurface
            aodHandles.OpticalSystem.SurfaceArray(removePosition-1).Stop = 1;
        else
            aodHandles.OpticalSystem.SurfaceArray(removePosition).Stop = 1;
        end
    end
    % The next selected row will be the one in the removed position, so if
    % it is image plane then dont let further removal
    if aodHandles.OpticalSystem.SurfaceArray(removePosition).ImageSurface
        CAN_REMOVE_SURFACE = 0;
    end
    parentWindow.ParentHandles = aodHandles;
    updateSurfaceOrComponentEditorPanel( parentWindow );
end
function ret = checkTheCurrentSystemDefinitionType(aodHandles)
    
    if get(aodHandles.popSystemDefinitionType,'Value') == 1
        ret = 1;
    else
        choice = questdlg(['Your system is not defined using Surface',...
            ' Based method. Editing in the surface editor',...
            ' window automatically converts your system to Surface',...
            ' Based definition. Do you want to continue editing?'], ...
            'Change System Definition Type', ...
            'Yes','No','No');
        % Handle response
        switch choice
            case 'Yes'
                set(aodHandles.popSystemDefinitionType,'Value', 2);
                ret = 1;
            case 'No'
                ret = 0;
                return;
        end
    end
end

