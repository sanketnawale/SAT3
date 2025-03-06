clear all; close all; clc;

%% 🌍 Simulation Parameters
tic
clear all; close all; clc;
Simulation_T = 110 * 60;  % Total simulation time (110 minutes in seconds)
Time_Step = 60;           % Time step (1 min = 60 sec)
MonteCarlo = 5000;        % Monte Carlo simulation (number of packets per node)
Nodes = 3;                % Number of ground nodes (Rome, Milan, NodeRM)
Pkct_Per_Hour = 100;      % Packets per hour for each node

%% 🌍 Ground Nodes (Rome, Milan, NodeRM)
Node_Coordinates = [ 
    41.9028, 12.4964;  % Rome (Node 1)
    45.4642, 9.1900;   % Milan (Node 2)
    41.9, 12.5         % NodeRM (Node 3, near Rome)
];




%% 🛰️ Satellite Constellation (Walker)
Sat_Per_Plane = 36;              % Satellites per plane
Num_Planes = 6;                 % Number of plan
Total_Sats = Sat_Per_Plane * Num_Planes;  % Total number of satellites in constellation
Orbital_Inclination = deg2rad(60);  % Inclination in radians
H = 1200e3;                    % Satellite altitude (meters)
Earth_Radius = 6378e3;           % Earth radius (meters)
Time_Vector = 0:Time_Step:Simulation_T;  % Simulation time vector


% 🛰️ Generate Walker Delta Constellation
oev = walker_delta(Sat_Per_Plane, Num_Planes, 1, pi, Earth_Radius + H, Orbital_Inclination);
Num_Satellites = size(oev, 1);
num_steps = length(Time_Vector);

%% 📡 Obtain Satellite Geometry
[Distances, Elevation_Angles, Ground_Distances, Visibility, Num_Visible_Sats, Sat_IDs, Latitudes, Longitudes] = ...
    Satellite_Geometry(H, Node_Coordinates, oev, Earth_Radius, Time_Vector);

%% 📊 Initialize Satellite Visibility Matrix
% Columns: [Time (min), Rome Visible Sats, Milan Visible Sats, NodeRM Visible Sats]
Visible_Sat_Matrix = zeros(length(Time_Vector), 4);

%% Additional storage for packet reception times info
% (For Rome and Milan, we build a string for each time step indicating
% the packet arrival (reception) times per visible satellite.
% For NodeRM, we store a simple message.)
Rome_PktReception = cell(length(Time_Vector), 1);
Milan_PktReception = cell(length(Time_Vector), 1);
NodeRM_PktReception = cell(length(Time_Vector), 1);

%% 📡 Random Access Logic for Packet Transmission, Reception & Collision Detection

% Initialize Success Rates and Collisions (for transmitting nodes)
SuccessRate = zeros(Nodes, length(Time_Vector));  % Success rate per node per time step
Collisions = zeros(Nodes, length(Time_Vector));     % Collision counts per node per time step

% Initialize Received Packets for NodeRM (relay node)
Received_Packets_NodeRM = zeros(1, length(Time_Vector));
num_nodes = size(Node_Coordinates, 1);  % Define number of ground nodes
Signal_Delay = zeros(num_nodes, Num_Satellites, num_steps);  % Initialize delay matrix


% Simulation loop over time steps
for t = 1:length(Time_Vector)
    current_time_min = Time_Vector(t) / 60;  % Convert time to minutes
    fprintf('\n⏳ Time %.2f min: \n', current_time_min);

    % Store visibility data in matrix (include NodeRM as 4th column)
    Visible_Sat_Matrix(t, :) = [current_time_min, Num_Visible_Sats(1, t), Num_Visible_Sats(2, t), Num_Visible_Sats(3, t)];

    % Print visible satellite info for each node
    fprintf('📡 Rome sees %d satellites: %s\n', Num_Visible_Sats(1, t), mat2str(Sat_IDs{1, t}));
    fprintf('📡 Milan sees %d satellites: %s\n', Num_Visible_Sats(2, t), mat2str(Sat_IDs{2, t}));
    fprintf('📡 NodeRM sees %d satellites: %s\n', Num_Visible_Sats(3, t), mat2str(Sat_IDs{3, t}));

    % Identify common satellites (Rome and Milan)
    Common_Sats = intersect(Sat_IDs{1, t}, Sat_IDs{2, t});
    if ~isempty(Common_Sats)
        fprintf('🔄 Common Satellites seen by Rome & Milan: %s\n', mat2str(Common_Sats));
    end

    % Container for NodeRM relay packet reception times
    NodeRM_Packet_Times = [];
    
    % For transmitting nodes (Rome = node 1; Milan = node 2)
    for n = 1:Nodes-1
        % Build a string to store packet reception times per satellite for this node
        pkt_str = '';
        if Num_Visible_Sats(n, t) == 0
            if n == 1
                pkt_str = 'None';
                Rome_PktReception{t} = pkt_str;
            elseif n == 2
                pkt_str = 'None';
                Milan_PktReception{t} = pkt_str;
            end
            continue;
        end

        % Each transmitting node sends a batch of packets during the time step
        Num_Packets = 1000;
        Transmission_Times = rand(1, Num_Packets) * Time_Step;  % Random transmission times (in sec)
        sorted_times = sort(Transmission_Times);

        % Get list of satellites visible for this node in the current time step
        Visible_Sats = Sat_IDs{n, t};
        if isempty(Visible_Sats)
            fprintf('🚫 No satellites visible for Node %d at %.2f min, skipping transmission.\n', n, current_time_min);
            if n == 1, Rome_PktReception{t} = 'None'; end
            if n == 2, Milan_PktReception{t} = 'None'; end
            continue;
        end

        % Initialize reception storage for each satellite over the whole constellation
        Sat_Receive_Times = cell(Total_Sats, 1);
        for s = 1:Total_Sats
            Sat_Receive_Times{s} = [];
        end

        % Assign each packet to every visible satellite (if the satellite ID is valid)
        for pkt = 1:Num_Packets
            for chosen_sat = Visible_Sats
                if chosen_sat <= Total_Sats
                    arrival_time = sorted_times(pkt) + Signal_Delay(n, chosen_sat, t);
                    Sat_Receive_Times{chosen_sat} = [Sat_Receive_Times{chosen_sat}, sorted_times(pkt)];
                end
            end
        end

        % Now, for each visible satellite, print and store the packet reception times
        for s = Visible_Sats
            if s <= Total_Sats
                if ~isempty(Sat_Receive_Times{s})
                    sat_arrivals = sort(Sat_Receive_Times{s});
                    formatted_arr = ['[', strtrim(num2str(sat_arrivals, '%.2f ')), ']'];
                    fprintf('⏰ Node %d, Satellite %d arrival packet timings (within %.2f sec): %s\n', ...
                    n, s, Time_Step, formatted_arr);
                    pkt_str = [pkt_str sprintf('Sat %d: %s; ', s, formatted_arr)];

                else
                    fprintf('⏰ Node %d, Satellite %d: No packet arrivals during this time step.\n', n, s);
                    pkt_str = [pkt_str sprintf('Sat %d: []; ', s)];
                end
            else
                fprintf('⏰ Node %d: Satellite %d ID exceeds the allocated satellite count (%d).\n', n, s, Total_Sats);
            end
        end
        
        % Save the packet reception string for the transmitting node
        if n == 1
            Rome_PktReception{t} = pkt_str;
        elseif n == 2
            Milan_PktReception{t} = pkt_str;
        end

        % Now check for collisions and update SuccessRate and Collisions
        for s = 1:Total_Sats
    if ~isempty(Sat_Receive_Times{s})
        % Sort arrival times (not just transmission times)
        sat_arrival_times = sort(Sat_Receive_Times{s}); 

        % ✅ Check for collisions based on 10 ms threshold
        collisions = sum(diff(sat_arrival_times) < 0.01);  
        total_packets = length(sat_arrival_times);

        % ✅ Store results
        Collisions(n, t) = collisions;
        SuccessRate(n, t) = total_packets - collisions;

        % ✅ Relay Logic for NodeRM (Rome to NodeRM)
        if n == 1
            if SuccessRate(n, t) > 0
                % **Only relay successful (non-collided) packets**
                non_collided_packets = sat_arrival_times(collisions+1:end);
                NodeRM_Packet_Times = [NodeRM_Packet_Times, non_collided_packets]; 
            end
        end
    end
end

% ✅ Debugging Output
fprintf('📊 Node %d transmitted %d packets, %d collisions\n', n, Num_Packets, Collisions(n, t));
    end

    % 🚀 NodeRM Packet Reception (Relay) Logic
    if ~isempty(NodeRM_Packet_Times)
        NodeRM_Packet_Times = sort(NodeRM_Packet_Times);
        NodeRM_Visible_Sats = Sat_IDs{3, t};  
        if ~isempty(intersect(NodeRM_Visible_Sats, Sat_IDs{1, t}))
            Received_Packets_NodeRM(t) = 1;  % Accept first arriving packet as successful relay
            first_packet_time = NodeRM_Packet_Times(1);
            if length(NodeRM_Packet_Times) > 1
                fprintf('🚨 COLLISION at NodeRM: %d packets received, but only 1 is accepted (Time: %.2f min)\n', ...
                    length(NodeRM_Packet_Times), first_packet_time/60);
            else
                fprintf('📡 NodeRM successfully received a packet at %.2f min\n', first_packet_time/60);
            end
            NodeRM_PktReception{t} = sprintf('Received: %s', mat2str(NodeRM_Packet_Times));
        else
            fprintf('📡 NodeRM Received No Packets at %.2f min (No overlapping visible satellites).\n', current_time_min);
            NodeRM_PktReception{t} = 'None';
        end
    else
        fprintf('📡 NodeRM Received No Packets at %.2f min (No successful relay).\n', current_time_min);
        NodeRM_PktReception{t} = 'None';
    end

end  % End simulation loop

%% Display the Visibility Matrix table
disp(array2table(Visible_Sat_Matrix, 'VariableNames', {'Time_Min','Rome_Sats','Milan_Sats','NodeRM_Sats'}));

%% Build a detailed table containing visible satellite IDs and packet reception times
Time_Min = round((Time_Vector)' / 60, 1);

% Convert the raw Sat_IDs cells into readable strings
Rome_Sat_IDs_str = cellfun(@(x) mat2str(x), Sat_IDs(1,:)', 'UniformOutput', false);
Milan_Sat_IDs_str = cellfun(@(x) mat2str(x), Sat_IDs(2,:)', 'UniformOutput', false);
NodeRM_Sat_IDs_str = cellfun(@(x) mat2str(x), Sat_IDs(3,:)', 'UniformOutput', false);

DetailedTable = table(Time_Min, Rome_Sat_IDs_str, Milan_Sat_IDs_str, NodeRM_Sat_IDs_str, ...
    Rome_PktReception, Milan_PktReception, NodeRM_PktReception, ...
    'VariableNames', {'Time_Min','Rome_Sat_IDs','Milan_Sat_IDs','NodeRM_Sat_IDs', ...
    'Rome_PktReception','Milan_PktReception','NodeRM_PktReception'});
disp(DetailedTable);

%% Save DetailedTable to Excel with a specified path
% Specify the folder path (change this to your desired folder)
folderPath = 'D:\thesis\walker\Analysis-and-Simulation-of-LoRaWAN-LR-FHSS-main (1)\Analysis-and-Simulation-of-LoRaWAN-LR-FHSS-main';
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
filename = fullfile(folderPath, 'DetailedResults.xlsx');
writetable(DetailedTable, filename, 'Sheet', 'DetailedResults');
fprintf('Detailed results saved to: %s\n', filename);

%% 📊 Final Summary and Graphs
fprintf('\n==== Final Results ====\n');
fprintf('✅ Overall Success Rate for Rome: %.2f%%\n', mean(SuccessRate(1, :)) / Num_Packets * 100);
fprintf('✅ Overall Success Rate for Milan: %.2f%%\n', mean(SuccessRate(2, :)) / Num_Packets * 100);
fprintf('📡 Total Packets Successfully Received by NodeRM: %d\n', sum(Received_Packets_NodeRM));
toc;

%% GRAPH PLOTTING

% Create a common time axis (in minutes)
time_minutes = Visible_Sat_Matrix(:, 1);

% Graph 1: Collisions Over Time for Rome (Node 1) and Milan (Node 2)
figure;
plot(time_minutes, Collisions(1, :), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(time_minutes, Collisions(2, :), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Time (min)'); ylabel('Number of Collisions');
title('Collisions over Time for Rome and Milan');
legend('Rome (Node 1)', 'Milan (Node 2)'); grid on;

% Graph 2: Successful Transmissions Over Time for Rome and Milan
figure;
plot(time_minutes, SuccessRate(1, :), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(time_minutes, SuccessRate(2, :), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Time (min)'); ylabel('Successful Transmissions');
title('Successful Transmissions over Time for Rome and Milan');
legend('Rome (Node 1)', 'Milan (Node 2)'); grid on;

% Graph 3: NodeRM Packet Reception over Time
figure;
stem(time_minutes, Received_Packets_NodeRM, 'g', 'LineWidth', 1.5, 'Marker', 'o');
xlabel('Time (min)'); ylabel('NodeRM Reception (1 = Received)');
title('NodeRM Packet Reception over Time'); grid on;



% ✅ Check if Latitudes and Longitudes exist before plotting
if exist('Latitudes', 'var') && exist('Longitudes', 'var')
    figure;
    hold on;
    grid on;
    
    % 🌍 **Plot Ground Stations**
    scatter(Node_Coordinates(:, 2), Node_Coordinates(:, 1), 100, 'ro', 'filled');  % Red markers for nodes
    text(Node_Coordinates(:, 2) + 1, Node_Coordinates(:, 1), {'Rome', 'Milan', 'NodeRM'}, 'FontSize', 12);

    xlabel('Longitude (°)');
    ylabel('Latitude (°)');
    title('2D Animated Ground Tracks of Satellites');
    xlim([-180, 180]);
    ylim([-90, 90]);

    % 🌟 **Initialize Satellite Plot Objects**
    sat_plots = gobjects(Num_Satellites, 1);
    ground_tracks = gobjects(Num_Satellites, 1);

    for s = 1:Num_Satellites
        % **Plot empty ground track (will update over time)**
        ground_tracks(s) = plot(Longitudes(s, 1), Latitudes(s, 1), 'b--', 'LineWidth', 1); % Dashed line for ground track
        sat_plots(s) = plot(Longitudes(s, 1), Latitudes(s, 1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Blue circles for satellites
    end

    % 🛰️ **Animate the Satellite Movement**
    for t = 1:num_steps
        for s = 1:Num_Satellites
            % Update satellite positions in the plot
            set(sat_plots(s), 'XData', Longitudes(s, t), 'YData', Latitudes(s, t));

            % **Update ground track by plotting past positions**
            set(ground_tracks(s), 'XData', Longitudes(s, 1:t), 'YData', Latitudes(s, 1:t));
        end

        % 📌 **Update Plot Title with Time**
        title(sprintf('2D Animated Ground Tracks of Satellites (Time: %.2f min)', Time_Vector(t) / 60));

        pause(0.1);  % Small pause for animation effect
    end

    hold off;
else
    fprintf('⚠️ Warning: Latitudes and Longitudes are not available for plotting.\n');
end
