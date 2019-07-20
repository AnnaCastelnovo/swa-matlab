function [Info, SW] = swa_FindSWTravelling(Info, SW, indSW, flag_wait)
% function which calculate the travelling parameters for each slow wave given the
% delay maps previously calculated in swa_FindSWChannels

% check inputs
if nargin < 4
    flag_wait = 1;
end

if nargin < 3
    indSW = [];
    flag_wait = 0;
end

% check for empty structure
if length(SW) < 1
    fprintf(1, 'Warning: Wave structure is empty, you must find waves in the reference first \n');
    return
elseif length(SW) < 2
    if isempty(SW.Ref_Region)
        fprintf(1, 'Warning: Wave structure is empty, you must find waves in the reference first \n');
        return
    end
end

% check for sufficient parameter inputs
if ~isfield(Info.Parameters, 'Travelling_GS');
    Info.Parameters.Travelling_GS = 20; % size of interpolation grid
    fprintf(1,'Information: Interpolation grid set at 20x20 by default. \n');
    Info.Parameters.Travelling_MinDelay = 40; % minimum travel time (ms)
end

if ~isfield(Info.Parameters, 'Travelling_RecalculateDeley')
    Info.Parameters.Travelling_RecalculateDelay = true; % consistency with previous versions
end

%% Check Electrodes for 2D locations (match to grid)

if Info.Parameters.Travelling_RecalculateDelay
    Info.Electrodes = swa_add2dlocations(Info.Electrodes, Info.Parameters.Travelling_GS);
    fprintf(1,'Calculation: 2D electrode projections (Info.Electrodes). \n');
    
    xloc = [Info.Electrodes.x]; xloc=xloc(:);
    yloc = [Info.Electrodes.y]; yloc=yloc(:);
    
    % create the plotting mesh
    GS = Info.Parameters.Travelling_GS;
    XYrange = linspace(1, GS, GS);
    XYmesh = XYrange(ones(GS,1),:);
    
    % check Matlab version for interpolant...
    if exist('scatteredInterpolant', 'file')
        % If its available use the newest function
        F = scatteredInterpolant(xloc, yloc, SW(1).Travelling_Delays(:),...
            'natural', 'none');
        interp_version = 1;
    else
        % Use the old function
        F = TriScatteredInterp(xloc, yloc, SW(1).Travelling_Delays(:),...
            'natural');
        interp_version = 0;
    end
    
else
    
    % check if the delay map is already included in the SW structure
    if isempty(SW(1).Travelling_DelayMap)
        error('User requested no map calculation, but no maps were found in the structure');
    end
    
    % calculate necessary parameters
    XYrange = linspace(1, ...
        size(SW(1).Travelling_DelayMap, 1), size(SW(1).Travelling_DelayMap, 2));

end

    %% ~~~~~~~~~ ADDED BY ANNA
    
    % user inputs: technically to add outside this function
    C3_index = 59; % C3 in 185 chanlocs
    C4_index = 155; % C4 in 185 chanlocs
    measurement_input = 56;
    
    % find scaling factor for real head dimensions 
    switch measurement_input
        case 47
            measured_C3_C4_cm= 12;
        case 51
            measured_C3_C4_cm= 13;
        case 54
            measured_C3_C4_cm= 14;
        case 56
            measured_C3_C4_cm= 15;
        case 58
            measured_C3_C4_cm= 16;
    end
    C3_pos   = [Info.Electrodes(C3_index).X Info.Electrodes(C3_index).Y];
    C4_pos   = [Info.Electrodes(C4_index).X Info.Electrodes(C4_index).Y];
    loc_file_C3C4_distance = sqrt(sum((C3_pos - C4_pos).^2, 2));
    scaleF = measured_C3_C4_cm/loc_file_C3C4_distance;
    
    % chanlocs.XYZ are electrode positions in 'units'
    % x y z are scaled chanlocs (streateched to a head with circumference 
    % that we use to project the wave on a 'real' scalp (in cm))
    x = [Info.Electrodes.X]' * scaleF;
    y = [Info.Electrodes.Y]' * scaleF;
    z = [Info.Electrodes.Z]' * scaleF;
    
    % parameters for 2D transformation from egi chanlocs to 40*40 grid
    a  = min(x);
    x1 = x - a;
    b  = max(x1);
    c  = min(y);
    y1 = y - c;
    d  = max(y1);
    
    % find the function of the scaled surface
    Fspazio = scatteredInterpolant(x, y, z, 'linear');
    
    %% ~~~~~~~~~ 

%% Loop for each SW
if isempty(indSW)
    loopRange = 1:length(SW);
else
    loopRange = indSW;
end

if flag_wait
    h = waitbar(0,'Please wait...', 'Name', 'Finding Streams...');
end

for nSW = loopRange
    
    if Info.Parameters.Travelling_RecalculateDelay
        Delays      = SW(nSW).Travelling_Delays;
        
        % Check for minimum travel time...
        if max(Delays) < Info.Parameters.Travelling_MinDelay * Info.Recording.sRate/1000
            continue
        end
        
        % Interpolate delay map [zeros or nans above...]
        Delays = Delays(:); % ensure data is in column format
        
        % check interpolation function
        if interp_version
            F.Values = Delays;
        else
            F.V = Delays;
        end
        SW(nSW).Travelling_DelayMap = F(XYmesh, XYmesh'); % Delay map (with zeros)
        
        % Define Starting Point(s) on the GSxGS grid...
        sx = xloc(SW(nSW).Channels_Active);
        sy = yloc(SW(nSW).Channels_Active);
        
    else
        % get the starting x/y positions from the struct directly
        sx = Info.Electrodes.x(SW(nSW).Channels_Active);
        sy = Info.Electrodes.y(SW(nSW).Channels_Active);
    end
    
    % calculate the gradients over the delay map
    [u,v] = gradient(SW(nSW).Travelling_DelayMap);
    
    % Find Streamline(s)
    % ''''''''''''''''''
    % Use adstream2 
    % TODO: optimise by coding entire loop in C
    Streams         = cell(1,length(sx));
    Distances       = cell(1,length(sx));
    for n = 1:length(sx)
        % find streams backwards from current point
        [StreamsBack, DistancesBack,~] = adstream2b(...
            XYrange,XYrange,-u,-v,sx(n),sy(n), cosd(45), 0.1, 1000);
        % find streams forward from current point
        [StreamsForw, DistancesForw,~] = adstream2b(...
            XYrange,XYrange,u,v,sx(n),sy(n), cosd(45), 0.1, 1000);
        % combine the two directions for continuous stream
        Streams{n}      = [fliplr(StreamsBack), StreamsForw];
        Distances{n}    = [fliplr(DistancesBack), DistancesForw];
    end
       
    % Process and save streamlines...
    % '''''''''''''''''''''''''''''''
    Streams(cellfun(@isempty, Streams)) = []; %Remove empty streams
    Distances(cellfun(@isempty, Distances)) = []; %Remove empty streams
    
    if isempty(Streams) % continue if there were no streams found
        continue
    end
    
    % Minimum Distance Threshold (25% of longest path)
    tDist = cellfun(@(x) sum(x), Distances);    %% Plot Functions
    Streams(tDist < max(tDist)/4) = [];
    Distances(tDist < max(tDist)/4) = [];
    
        %% ~~~~~~~~~~ ADDED BY ANNA: LOCAL SPEED
    
    % create grid with delays (40x40)
    XYrange = linspace(1, Info.Parameters.Travelling_GS, Info.Parameters.Travelling_GS);
    XYmesh = XYrange(ones(Info.Parameters.Travelling_GS,1),:);
    x1 = reshape(XYmesh, 1, []);
    x2 = reshape(XYmesh', 1, []);
    ytime = reshape([SW(nSW).Travelling_DelayMap]', 1, []);
    
    % remove isnan
    x1(isnan(ytime)) = [];
    x2(isnan(ytime)) = [];
    ytime(isnan(ytime)) = [];    
    
    % interpolate delays on a gridded surface (and not on electrode surface)
    Ftime = scatteredInterpolant(x1',x2',ytime');                 
    % calculate the delays between each point of the traveling streams of
    % each wave and convert in seconds
    sLocalDelay = cellfun(@(x) abs(diff(Ftime(x')))/ Info.Recording.sRate, Streams, 'UniformOutput', false);  

%   put each stream back in 2d (or 3d) using converison factor from 40*40 grid to real electrode dimensions
    sLocalPosition = cellfun(@(x) 0.01 * ((x' * b)/(Info.Parameters.Travelling_GS-1) - 1 + a), Streams, 'UniformOutput', false);
%     TODO: 3dtest = cellfun(@(x) Fspazio(x(:,1), x(:,2)), sLocalPosition, 'UniformOutput', false); % calcolate the 3rd dimension

    % compute local speed (absolute travelled distance divided the delay)
    sLocalDistance = cellfun(@(x) sqrt(sum(diff(x).^2, 2)), sLocalPosition, 'UniformOutput', false);
    SW(nSW).LocalSpeed    = cellfun(@(x,y) x./y, sLocalDistance, sLocalDelay, 'UniformOutput',false);
%     SW(nSW).LocalSpeed    = cell2mat(cellfun(@(x) nanmedian(x), sLocalSpeed, 'UniformOutput', false));
    
    % compute total speed
    sTotalDistance = cell2mat(cellfun(@(x) sum(x), sLocalDistance, 'UniformOutput', false));
    sTotalTime = cell2mat(cellfun(@(x) sum(x), sLocalDelay, 'UniformOutput', false));
    SW(nSW).TotalSpeed = sTotalDistance./sTotalTime;
        
    %% ~~~~~~~~~~ 

   
    % Longest displacement
    tDisp = cellfun(@(x) (sum((x(:,1)-x(:,end)).^2))^0.5, Streams); % total displacement
    [~,maxDispId] = max(tDisp);
    SW(nSW).Travelling_Streams{1} = Streams{maxDispId};
    
    % Longest distance travelled (if different from displacement)
    tDist = cellfun(@(x) sum(x), Distances);    %% Plot Functions
    [~,maxDistId] = max(tDist);
    if maxDistId ~= maxDispId
        SW(nSW).Travelling_Streams{end+1} = Streams{maxDistId};
    end  

    % Most different displacement angle compared to longest stream (at least 45 degrees)
    streamAngle = cellfun(@(x) atan2d(x(1,end)- x(1,1),x(2,end)-x(2,1)), Streams);
    [maxAngle,maxAngleId] = max(streamAngle - streamAngle(maxDispId));
    if maxAngle > 45 || maxAngleId ~= maxDispId || maxAngleId ~= maxDistId
        SW(nSW).Travelling_Streams{end+1} = Streams{maxAngleId};
    end

    % Update waitbar
    if flag_wait
        waitbar(nSW/length(SW),h,sprintf('Slow Wave %d of %d',nSW, length(SW)))
    end
    
end

% DELETE the waitbar; don't try to CLOSE it.
if flag_wait
    delete(h);
end
