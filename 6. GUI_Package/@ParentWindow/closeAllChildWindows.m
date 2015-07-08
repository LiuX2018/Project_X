function [ childsclosed ] = closeAllChildWindows(parentWindow)
    %CloseAllChildrens: Closes all child windows opened except the main window
    % Member of ParentWindow class

    nChildren = size(parentWindow.ParentHandles.ChildWindows,2);
    for k = 1 : nChildren
        curentChildWindow = parentWindow.ParentHandles.ChildWindows(k);
        if isvalid(curentChildWindow) % Undeleted valid object
            closeChild(curentChildWindow);
        end
    end
    childsclosed = 1;
end