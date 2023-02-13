clc
clear
net = WirelessNetwork;
net.channel_sets = [4 4 5 5 5 5 5];
total_channels = sum(net.channel_sets); % Should be 33

% 25 is our maximum number of base stations we can afford 
location1 = 375+1i*225;
net.AddBase(location1, 1, 39);

location2 = 321+1i*388;
net.AddBase(location2, 4, 40);

location3 = 326+1i*284;
net.AddBase(location3, 3, 39);

location4 = 440+1i*211;
net.AddBase(location4, 4, 39);

location5 = 505+1i*304;
net.AddBase(location5, 2, 40);

location6 = 120+1i*255;
net.AddBase(location6, 7, 40);

location7 = 200+1i*100;
net.AddBase(location7, 2, 39);

location8 = 320+1i*165;
net.AddBase(location8, 2, 39);

location9 = 100+1i*500;
net.AddBase(location9, 2, 39);

location10 = 375+1i*540;
net.AddBase(location10, 2, 40);

location11 = 530+1i*150;
net.AddBase(location11, 1, 34);

location12 = 230+1i*570;
net.AddBase(location12, 4, 39);

location13 = 260+1i*300;
net.AddBase(location13, 1, 35);

location14 = 100+1i*100;
net.AddBase(location14, 5, 39);

location15 = 440+1i*450;
net.AddBase(location15, 1, 40);

location16 = 290+1i*250;
net.AddBase(location16, 5, 36);

location17 = 230+1i*230;
net.AddBase(location17, 6, 39);

location18 = 100+1i*200;
net.AddBase(location18, 3, 40);

location19 = 180+1i*200;
net.AddBase(location19, 4, 40);

location20 = 230+1i*60;
net.AddBase(location20, 3, 40);

location21 = 245+1i*445;
net.AddBase(location21, 6, 40);

location22 = 165+1i*360;
net.AddBase(location22, 5, 40);

location23 = 412+1i*300;
net.AddBase(location23, 7, 40);

location24 = 375+1i*115;
net.AddBase(location24, 6, 40);

location25 = 258+1i*130;
net.AddBase(location25, 7, 40);


net.ComputeCost
cst = net.ComputeCost;
% The coverage map
figure
net.map_resolution = 600;
net.MapCoverage;
% The SINR outage map 
figure
net.MapSINR 
figure
net.MapCity
load('MobileLocationData');
% City map with mobiles at 19:00 superimposed
net.MapCity(ms_locations{19});
% Coverage map with mobiles at 19:00 superimposed
net.MapMobiles(ms_locations{19})
[num_connected, num_over_thresh] = net.AnalyzeNetwork(ms_locations);
% A plot of the number of connected users, the number of call attempts,
% and the number over threshold for each hour of the day.
figure
plot( 1:length(call_attempts), call_attempts, '-k', ...
    1:length(call_attempts), num_over_thresh, '-r', ...
    1:length(call_attempts), num_connected, '-b' );
xlabel( 'time (in hours)' );
ylabel( 'number calls' );
legend( 'call attempts', 'number covered', 'number connected' );
fraction_over_thresh = num_over_thresh./call_attempts;
fraction_connected = num_connected./call_attempts;
% A plot of the fraction of users that are covered and the fraction
% that are covered and connected for each hour of the day.
figure
plot( 1:length(call_attempts), fraction_over_thresh, '-r', ...
    1:length(call_attempts), fraction_connected, '-b');
xlabel( 'time (in hours)' );
ylabel( 'fraction' );
legend( 'fraction covered', 'fraction connected');
% The percentage connected.
sum(num_connected)
sum(call_attempts)
percent_connected = sum(num_connected)/sum(call_attempts);
fprintf( 'The percentage connected is %2.2f%%\n', 100*percent_connected)
% The cost of the network (in dollars).
cons = sum(num_connected);
att = sum(call_attempts);
fprintf( 'The cost is $%d\n', cst)
% The total number of users connected throughout the day (the main KPI). 
fprintf( 'The number connected is %d\n', cons)
fprintf( 'The call attempts are %d\n', att)
net.Save( "PreliminaryDesign.mat")