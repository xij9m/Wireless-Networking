classdef WirelessNetwork < handle
    % WirelessNetwork class, version 4.7.  Sept. 19, 2022.
    %   An object whose properties specify a wireless network and whose methods analyze the network
    %   usage: net = WirelessNetwork( [bs_locations] [bs_channels], [bs_power_dBm] )
    %          bs_locations = column vector containing locations of the base stations (as complex numbers)
    %          bs_channels = matrix of size num_bs by num_sec containing the channel sets used at each bs sector
    %          bs_power_dBm = matrix of size num_bs by num_sec containing the power (in dBm) transmitted by each bs sector
    %
    %   notes: Run with no argument to create an empty object: net = WirelessNetwork;
    %          Once created, can use the "CreateFirstTier" method to create a network with 1 co-channel tier of interferers:
    %          net.CreateFirstTier( N ), where N = 1, 3, 4, or 7.

    properties
        % these are common variables, but not defined as formal properties:        
        % num_bs is the number of base stations
        % num_channelsets is the number of disjoint channel sets
        % num_sectors it the number of sectors
        % num_ms is the number of mobile stations 
        % ms_locations is a (row) vector containing the locations of the mobiles
                           
        bs_locations       % A length-num_bs (column) vector containing the locations of the base stations
        bs_channels        % A num_bs by num_sectors matrix containing the channel set used by each sector 
        bs_power_dBm       % A num_bs by num_sectors matrix containing the power of the sectors (in dBm)

        channel_sets       % A length-num_channelsets (row) vector containing number of RF channels
                           % in each channel set
                           
        sector_offset_angle = 0;   % offset angle for the sectors (in radians)
        map_resolution = 250;      % resolution in each of the x- and y-directions when plotting
        
        noise_dBm = -80;     % the noise floor in dBm (set to -Inf to turn off noise)
        
        pl_exponent = 3;      % path loss exponent, which can be overridden by most functions
        SINR_thresh = 10;     % SINR Threshold (in dB) 
        ms_per_channel = 10;  % number mobiles supported per RF channel 
        max_users;            % A num_bs by num_sectors matrix containing the maximum number of users supported by each sector
        map_filename = 'Morgantown';   % filename for the map
        save_filename = 'MyProject';   % default filename for the save file
        cost;

    end
    
    properties ( Hidden, SetAccess = private)
        new_fig = 1;
    end
    
    methods
        
        function obj = WirelessNetwork( varargin )
            % Constructor: Constructs the object
            
            if (nargin > 0)
                % first argument is the base station locations
                obj.bs_locations = varargin{1};
            
                % default channel allocation -- univeral frequency reuse and no sectorization
                obj.bs_channels = ones( size(obj.bs_locations) );
                
                % default power allocation
                obj.bs_power_dBm = ones( size(obj.bs_locations) );

            end
            
            if (nargin > 1)
                % second argument is the channel allocation
                % (if none specified assume universal frequency reuse)
                obj.bs_channels = varargin{2};
            end 
            
            if (nargin > 2)
                % third argument is the power allocation
                % if none specified, assume equal power
                obj.bs_power_dBm = varargin{3};
            end
            
        end
        
        function [bs_locations_cluster bs_channels_cluster bs_power_dBm_cluster] = CreateCluster( obj, N, center )
            % method CreateCluster creates a cluster of N cells centered at location "center"
            % usage: [bs_locations_cluster bs_channels_cluster bs_power_dBm_cluster]  = obj.CreateCluster( N, center )
            %        N = cluster size: number of cell/cluster (can be between 1 and 7, inclusive)
            %        center = center of the cluster
            %        bs_locations_cluster = locations of the base stations in the cluster
            %        bs_channels_cluster = channel assignment for the cluster
            %        bs_power_dBm_cluster = power of each cell in the cluster
            %
            
            % creates a prototype cluster
            H = sqrt(3);            % the height of a cell (assumes radius R=1)
            bs_channels_cluster = (1:N)';
            bs_power_dBm_cluster = zeros(N,1);
            bs_locations_cluster = zeros(N,1);
            
            % the centers of the cells in the cluster
            centers = [ 0;
                1.5 + 0.5*H*1i;
                H*1i;
                -1.5 + 0.5*H*1i;
                -1.5 - 0.5*H*1i;
                -H*1i;
                1.5 - 0.5*H*1i];
                    
            if (N >= 1)&&(N <= 7 )
                bs_locations_cluster = centers( 1:N );
            else
                printf( 'invalid cluster size' );
            end            
           
            center_real = ( max( real( bs_locations_cluster ) ) + min( real( bs_locations_cluster ) ) )/2;
            center_imag = ( max( imag( bs_locations_cluster ) ) + min( imag( bs_locations_cluster ) ) )/2;
            shift = center_real + center_imag*1i;
            
            % shift the center of the cluster to be at "center"
            bs_locations_cluster = bs_locations_cluster -shift + center;
            
        end
        
        function obj = CreateCentralCluster( obj, N, center )
            % method CreateCentralCluster creates a cluster of N cells centered at location "center"
            % usage: obj.CreateCluster( N, center )
            %        N = cluster size: number of cell/cluster (can be between 1 and 7, inclusive)
            %        center = center of the cluster
            [bs, channels, power_dBm] = obj.CreateCluster( N, center );
            obj.bs_locations( 1:N, 1 ) = bs;
            obj.bs_channels( 1:N, 1 ) = channels;
            obj.bs_power_dBm( 1:N, 1 ) = power_dBm;           
        end
        
        function obj = CreateFirstTier( obj, N )
            % method CreateFirstTier creates a network of hexagonal cells
            % containing a central cluster and one tier of interfering
            % clusters.  The total number of clusters created is 7.
            % N cells
            % usage: obj.CreateTiers( N )
            %        N = cluster size (can be 1, 3, 4, or 7)
                       
            % number of co-channel tiers
            tiers = 1; 
            
            % the height of a cell
            H = sqrt(3);
            
            % place bs into multiple tiers
            if (N==1)
                centers = [ 0;
                    +H*1i;
                    -H*1i;
                    +1.5 + 0.5*H*1i;
                    +1.5 - 0.5*H*1i;
                    -1.5 + 0.5*H*1i;
                    -1.5 - 0.5*H*1i];
            elseif (N == 3)
                centers = [ 0;
                    3
                    1.5 + 1.5*H*1i;
                    -1.5 + 1.5*H*1i;
                    -3
                    -1.5 - 1.5*H*1i;
                    1.5 - 1.5*H*1i ];
            elseif (N == 4)
                centers = [0;
                    3+H*1i;
                    2*H*1i;
                    -3+H*1i;
                    -3-H*1i;
                    -2*H*1i;
                    3-H*1i];
            elseif (N==7)
                centers = [0
                    -3-2*H*1i
                    -4.5+0.5*H*1i
                    -1.5+2.5*H*1i
                    +3+2*H*1i
                    +4.5-0.5*H*1i
                    +1.5-2.5*H*1i];
            end
           
            % clear values
            obj.bs_locations = [];
            obj.bs_channels = [];
            obj.bs_power_dBm = [];

            % place each cluster
            for i=1:length(centers)
                [bs, channels, power_dBm] = obj.CreateCluster( N, centers(i) );                
                obj.bs_locations( ((i-1)*N+1):i*N, 1 ) = bs;
                obj.bs_channels( ((i-1)*N+1):i*N, 1 ) = channels;
                obj.bs_power_dBm( ((i-1)*N+1):i*N, 1 ) = power_dBm;
            end        
        end
        
        function obj = Sectorize( obj, num_sec, varargin )
            % method Sectorize creates sectors
            % usage: obj.Sectorize( num_sec, [sector_offset_angle] )
            %        num_sec = number of sectors (any positive integer)
            %        [sector_offset_angle] = offset angle for the sectors, in radians (optional)
            
            % set the sector offset angle
            if (nargin > 2)
                obj.sector_offset_angle = varargin{1};
            else
                obj.sector_offset_angle = 0;
            end             
             
            % eliminate previous sector information
            obj.bs_power_dBm = obj.bs_power_dBm(:,1);
            obj.bs_channels = obj.bs_channels(:,1);
            
            % determine the number of channels used by sector 1
            num_channels = max( obj.bs_channels );
            
            % set the power and channels of the other sectors
            for i=2:num_sec
                obj.bs_power_dBm(:,i) = obj.bs_power_dBm(:,1);
                obj.bs_channels(:,i) = obj.bs_channels(:,1)+(i-1)*num_channels;
            end
            
            % clear the channel sets (need to ask users to re-enter)
            obj.channel_sets = [];
            
            % compute the cost (which also sets the maximum number of users)
            obj.ComputeCost;
        
        end
        
        function [SINR_dB,num_ms_per_sec] = ComputeSINR( obj, ms_locations,   varargin )
            % method ComputeSINR computes the SINR of mobiles at locations "ms_locations"
            % usage: SINRdB = obj.ComputeSINR( ms_locations, [pl_exponent], [channel] )
            %        ms_locations = row vector containing the locations of the mobiles (as complex numbers) 
            %        [pl_exponent] = path loss exponent (optional)
            %        [channel] = row vector containing the channel set to use for each mobile (optional)
            
            % if the pl_exponent is specified, then temporarily use it
            temp_pl_exponent = obj.pl_exponent;
            if nargin>2    
                obj.pl_exponent = varargin{1};                
            end
            
            % find powers and sectors
            [power_dBm, sector] = obj.FindPower( ms_locations );
            
            % if channel specified, then use it
            if nargin>3
                channel = varargin{2};                
            else
                % find channel at each location
                [channel,best_bs,best_bs_power_dBm] = obj.FindChannel( power_dBm, sector );
            end
            
            % determine number of bs and sectors
            [num_bs,num_sec] = size( obj.bs_channels );
            
            % initialize the CCI mask (used to zero out non-interferers)
            CCI_mask = zeros( size( power_dBm) );
            
            for j=1:num_bs
                % determine channels used by this bs in direction of each mobile
                channel_toward_ms = obj.bs_channels(j, sector(j,:) );
                
                % determine which ms are using the same channel
                CCI_mask(j,:) = logical( channel == channel_toward_ms );
            end
            
            % convert power to mW
            power_mW = 10.^(power_dBm/10);         
             
            % if the channel was specified, need to identify the strongest bs at that channel
            if nargin>3
                [best_bs_power_mW,best_bs] = max( power_mW.*CCI_mask );
                best_bs_power_dBm = 10*log10( best_bs_power_mW );
            end
            
            for j=1:num_bs
                % delete the ones served by this bs
                CCI_mask( j, best_bs == j ) = 0;
            end
            
            % add noise
            CCI_mW = power_mW.*CCI_mask + 10^(obj.noise_dBm/10);
            CCI_total_mW = sum( CCI_mW );
            
            % compute the SINR, in dB
            SINR_dB = best_bs_power_dBm -10*log10( CCI_total_mW ); 
            
            % detemine number mobiles in each sector that are above
            % threshold
            num_ms_per_sec = zeros( num_bs, num_sec );
            for j=1:num_bs
                for i=1:num_sec
                    num_ms_per_sec(j,i) = sum( and( SINR_dB > obj.SINR_thresh, and( best_bs == j, sector(j,:) == i ) ) );
                end
            end
            
            % restore the value of pl_exponent
            obj.pl_exponent = temp_pl_exponent;
 
           
        end        
        
        function [channel, best_bs, best_bs_power_dBm ] = FindChannel( obj, power_dBm, sector )
            % method FindChannel finds the channel used between each bs and ms, and finds the best channel for each ms
            % usage: [channel, best_bs, best_bs_power_dBm ] = obj.FindChannel( power_dBm, sector )
            %        power_dBm = a num_bs by num_ms matrix containing the received power at each ms due to each bs
            %        sector = a numb_bs by num_ms matrix containing the sector used by each bs to reach each ms
            %        channel = a row vector containing the best channel for each ms
            %        best_bs = a row vevctor containing the best bs for each ms
            %        best_bs_power_dBm = a row vector containing the power of the best bs for each ms (in dBm)
            
            % determine the best_bs
            [best_bs_power_dBm, best_bs] = max( power_dBm );
            [num_bs, num_ms] = size( power_dBm );
            
            % find the channel for each ms (for future improvement, should vectorize this part)
            channel = zeros(1,num_ms); 
            for i=1:num_ms
                channel(i) = obj.bs_channels( best_bs(i), sector(best_bs(i), i ) );
            end

            
        end
        
        function [power_dBm, sector] = FindPower( obj, ms_locations, varargin )
            % method FindPower finds a matrix specifying the recieved power at each ms due to each bs,
            % and the sector of each bs facing each ms
            % usage: [power_dBm, sector] = obj.FindPower( ms_locations,[pl_exponent] )
            %        ms_locations = row vector containing the locations of the mobiles (as complex numbers) 
            %        [pl_exponent] = path loss exponent (optional argument)  
            %        power_dBm = a num_bs by num_ms matrix containing the received power at each ms due to each bs
            %        sector = a numb_bs by num_ms matrix containing the sector used by each bs to reach each ms
            
            % if the pl_exponent is specified, then use it
            if nargin>2
                obj.pl_exponent = varargin{1};
            end
            
            % determine numer of base stations and mobile stations
            [num_bs,num_sec] = size( obj.bs_channels );    
            num_ms = length( ms_locations );
           
            % sector division
            sector_angles =  (0:(num_sec-1))*2*pi/num_sec + obj.sector_offset_angle;
            sector_dividers_complex = exp( 1i*sector_angles );
            sector_dividers = sort( angle( sector_dividers_complex ) );
            
            % initialize
            sector = num_sec*ones( num_bs, num_ms );
            power_dBm = zeros( num_bs, num_ms );          
         
            for i=1:num_bs
                % shift the origin to be centered at the bs
                shifted_bs = ms_locations - obj.bs_locations(i);
                
                % determine angle from this base to each mobile
                bs_ms_angle = angle( shifted_bs );
                
                % find the sector it belongs to
                for j=1:num_sec %num_sectors
                    sector(i, bs_ms_angle >= sector_dividers(j) ) = j;
                end
                
                distance = abs( shifted_bs );  
                
                % the power transmitted by this base station to each mobile
                Pt_dBm = obj.bs_power_dBm( i, sector(i,:) );
                
                % the received power at the mobile
                power_dBm(i,:) = Pt_dBm - 10*obj.pl_exponent*log10( distance );
            end            

        end
        
        function obj = FindMaxUsers( obj, varargin )
            % method FindMaxUsers determines the maximum number of users supported by each sector
            % usage: obj.FindMaxUsers( [channel_sets] )
            %        channel_sets = vector specifying how many RF channels are in each channel set (optional)
           
            % if channel_sets not correctly set, the set to a default of 1 RF channel per set
            num_channels = max( max( obj.bs_channels ) );
            
            % determine number of bs and sectors
            [num_bs,num_sec] = size( obj.bs_channels );
            
            if (nargin > 1)
                % set the channels
                obj.channel_sets = varargin{1};
            end
            
            if ( length( obj.channel_sets ) ~= num_channels )
                % fprintf( 'obj.channel_sets unset or wrong length, setting to default\n' );
                obj.channel_sets = ones(1,num_channels);
            end
            
            % now set the maximum number of users supported by each sector
            obj.max_users = zeros( num_bs, num_sec );
            for i=1:num_bs
                for j=1:num_sec
                    obj.max_users(i,j) = obj.ms_per_channel*obj.channel_sets( obj.bs_channels( i, j ) );
                end
            end
            
        end
        
        function [num_connected, num_over_thresh] = AnalyzeNetwork( obj, ms_locations, varargin )
            % method AnalyzeNetwork will determine the number of connected users and number over the threshold
            % usage: [num_connected, num_over_thresh] = obj.AnalyzeNetwork( ms_locations, [pl_exponent] )
            %        num_connected = number of mobiles that are above the threshold AND can obtain a channel
            %        num_over_thresh = number of mobiles over the SINR threshold
            %        ms_locations = mobile locations, in a vector or cell array
            %        [pl_exponent] = path loss exponent (optional argument)
            
            % if the pl_exponent is specified, then use it
            if nargin>2
                obj.pl_exponent = varargin{1};
            end
            
            % if the input ms_locations is not already a cell array, make it one
            if ~iscell( ms_locations )
                ms_in = ms_locations;
                ms_locations = cell(1);
                ms_locations{1} = ms_in;
            end
            
            % update the cost (moved here on 10-11-12; was previously below the following for loop)
            obj.ComputeCost;
                       
            for i=1:length( ms_locations )
                % compute the SINR for this realization of the mobile locations
                [SINR_dB, num_ms_per_sec] = obj.ComputeSINR( ms_locations{i} );
                
                % determine how many are over the threshold
                num_over_thresh(i) = sum (SINR_dB > obj.SINR_thresh );
                
                % determine how many can connected
                num_connected(i) = sum( sum( min( num_ms_per_sec, obj.max_users ) ) );
            end
             
                 
        end
        
        function cost_out = ComputeCost( obj )
            % method ComputeCost computes the cost (in thousands of dollars) of the network
            % usage cost = obj.ComputeCost
            
            % determine number of bs and sectors
            [num_bs,num_sec] = size( obj.bs_channels );
            
            % make sure to set the maximum number of users
            obj.FindMaxUsers;
            
            % determine number of channels
            num_channels = 0;
            for i=1:num_bs
                for j=1:num_sec
                    num_channels = num_channels + obj.channel_sets( obj.bs_channels(i,j) );
                end
            end
            
            cost_per_channel = 25;            
            cost_per_bs = 350;
            
            cost_out = num_channels*cost_per_channel + num_bs*cost_per_bs;
            
            obj.cost = cost_out;
            
        
        end
        
        function MapCoverage( obj, varargin )
            % method MapCoverage shows a coverage map for a cellular wireless network.
            %
            % usage: obj.MapCoverage( [pl_exponent], [w e s n] )
            %        [pl_exponent] = path loss exponent (optional argument)  
            %        [w e s n] = coorinates bounding the map (optional argument)
                       
            number_channels = max( max( obj.bs_channels) );
            
            % if the pl_exponent is specified, then use it
            if (nargin>1)
                obj.pl_exponent = varargin{1};
            end           
            
            % compute the minimum separation between base stations (assuming more than 1 bs)
            num_bs = length( obj.bs_locations );
            if (num_bs > 1)           
                min_dist = min( abs( obj.bs_locations(2:end) - obj.bs_locations(1) ) );
            else
                min_dist = 1;
            end
            
            if (nargin > 2)
                west = varargin{2}(1);
                east = varargin{2}(2);
                south = varargin{2}(3);
                north = varargin{2}(4);
            else
                west = min( real( obj.bs_locations ) )-min_dist/2;
                east = max( real( obj.bs_locations ) )+min_dist/2;
                south = min( imag( obj.bs_locations ) )-min_dist/2;
                north = max( imag( obj.bs_locations ) )+min_dist/2;
            end
            
            grid_width = east - west;
            grid_height = north - south;
            vertical_spacing = grid_height/obj.map_resolution;
            horizontal_spacing = grid_width/obj.map_resolution;
            
            x_temp = west:horizontal_spacing:east;
            y_temp = south:vertical_spacing:north;
            
            % create an array of coordinates
            x_coord = ones(length(x_temp),1)*x_temp;
            y_coord = y_temp'*ones(1,length(y_temp));
            coord = x_coord + 1i*y_coord;
            ms_locations = reshape( coord, 1, length(x_temp)*length(y_temp) );
           
            % find powers
            [power_dBm, sector] = obj.FindPower( ms_locations );
            
            % find channel at each location
            channel = obj.FindChannel( power_dBm, sector );
            channel_matrix = reshape( channel, length(y_temp), length(x_temp) );
            
            % figure
            cmap = colormap( jet(number_channels) );
            colormap( cmap );
            image( x_temp, y_temp, channel_matrix );
            hold on
            plot( obj.bs_locations, 'k*' );
            axis xy % xy mode required
            
            for k=1:num_bs
                text( real( obj.bs_locations(k) ), imag( obj.bs_locations(k) ), num2str( k ), ...
                    'FontSize', 10, ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'top' );
            end
            
            hold off
                       
        end

  
        function MapSINR( obj, varargin )
            % method MapSINR shows an SINR map for a cellular wireless network.
            %
            % usage: obj.MapSINR( [pl_exponent], [w e s n], [thresh] )
            %        [pl_exponent] = path loss exponent (optional argument)  
            %        [w e s n] = coordinates bounding the map (optional argument)  
            %        [thresh]    = SINR threshold
           
            % if the pl_exponent is specified, then use it
            if (nargin>1)
                obj.pl_exponent = varargin{1};
            end           
                       
            % determine the boundaries of the map
            if (nargin > 2)
                temp = varargin{2};
                west = temp(1);
                east = temp(2);
                south = temp(3);
                north = temp(4);
            else        
                % compute the minimum separation between base stations (assuming more than 1 bs)
                num_bs = length( obj.bs_locations );
                if (num_bs > 1)
                    min_dist = min( abs( obj.bs_locations(2:end) - obj.bs_locations(1) ) );
                else
                    min_dist = 1;
                end
                
                west = min( real( obj.bs_locations ) )-min_dist/2;
                east = max( real( obj.bs_locations ) )+min_dist/2;
                south = min( imag( obj.bs_locations ) )-min_dist/2;
                north = max( imag( obj.bs_locations ) )+min_dist/2;
            end
            
            grid_width = east - west;
            grid_height = north - south;
            vertical_spacing = grid_height/obj.map_resolution;
            horizontal_spacing = grid_width/obj.map_resolution;
            
            x_temp = west:horizontal_spacing:east;
            y_temp = south:vertical_spacing:north;
            
            % create an array of coordinates
            x_coord = ones(length(x_temp),1)*x_temp;
            y_coord = y_temp'*ones(1,length(y_temp));
            coord = x_coord + 1i*y_coord;
            ms_locations = reshape( coord, 1, length(x_temp)*length(y_temp) );
           
            % find SINR
            SINR_dB = obj.ComputeSINR( ms_locations );
            
            % reshape so it can be an image
            SINR_matrix = reshape( min( SINR_dB, 128 ), length(y_temp), length(x_temp) );
            
            % the threshold is either an input argument or a property
            if (nargin > 3)
                thresh = varargin{3};
            else
                thresh = obj.SINR_thresh;
            end
            
            mymap = [1 0 0
                0 0 1];
            
            SINR_matrix = (SINR_matrix > thresh);
            colormap( mymap );
            image( x_temp, y_temp, SINR_matrix );
            axis xy % xy mode required
            covered_area = sum( sum( SINR_matrix ) )/prod( size( SINR_matrix ) );
            fprintf( 'Covering %2.2f percent of the area with sufficient SINR\n', 100*covered_area );
            fprintf( 'worst case SINR is %f\n', min( min( SINR_dB ) ) );
            
        end
        
        function Recording = MapCity( obj, varargin )
            % method MapCity shows a map of the city, with the base station locations overlaid
            % and (optionally), the locations of the mobiles
            %
            % usage: Recording = obj.MapCity( [ms_locations], [pl_exponent] )
            %        Recording = a Matlab movie showing the locations of the mobiles
            %        [ms_locations] = mobile locations, in a vector or cell array
            %        [pl_exponent] = path loss exponent (optional argument)
            
            
            % determine the number of base stations
            [num_bs,num_sec] = size( obj.bs_channels );
            
            % load the morgantown map
            MapImage = imread( obj.map_filename, 'jpg');
            MapSize = length(MapImage);
            
            % plot the locations of the base stations (as astericks)
            % figure;
            image( MapImage );
            hold on
            
            % plot the base stations
            for k=1:num_bs
                marker = ['*k'];
                % hack -- need to flip the imaginary axis (because image plots upside down)
                plot( real( obj.bs_locations(k) ), MapSize-imag( obj.bs_locations(k) ), marker );
                % text( real( obj.bs_locations(k) ), imag( obj.bs_locations(k) ), num2str( obj.bs_channels(k) ), ...
                text( real( obj.bs_locations(k) ), MapSize-imag( obj.bs_locations(k) ), num2str( k ), ...
                    'FontSize', 10, ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'top' );
            end
            
            Recording(1) = getframe;
            
            hold off
            
            % if MS are present, then plot them.
            
            % if  ms_locations are in the argument list, then use them
            if (nargin>1)
                ms_locations = varargin{1};
                
                % if the pl_exponent is specified, then use it
                if (nargin>2)
                    obj.pl_exponent = varargin{1};
                end
                
                % if the input ms_locations is not already a cell array, make it one
                if ~iscell( ms_locations )
                    ms_in = ms_locations;
                    ms_locations = cell(1);
                    ms_locations{1} = ms_in;
                end
                
                for i=1:length( ms_locations )
                    
                    % load the map
                    image( MapImage );
                    hold on
                    
                    % plot the base stations
                    for k=1:num_bs
                        marker = ['*k'];
                        % need to flip the imaginary axis (because image plots upside down)
                        plot( real( obj.bs_locations(k) ), MapSize-imag( obj.bs_locations(k) ), marker );
                        % text( real( obj.bs_locations(k) ), imag( obj.bs_locations(k) ), num2str( obj.bs_channels(k) ), ...
                        text( real( obj.bs_locations(k) ), MapSize-imag( obj.bs_locations(k) ), num2str( k ), ...
                            'FontSize', 10, ...
                            'HorizontalAlignment', 'left', ...
                            'VerticalAlignment', 'top' );
                    end
                    
                    % compute the SINR for this realization of the mobile locations
                    SINR_dB = obj.ComputeSINR( ms_locations{i} );
                    
                    % mobiles over the threshold are indicated in blue
                    over_thresh = ms_locations{i}( find( SINR_dB >= obj.SINR_thresh) );
                    plot( real( over_thresh), MapSize-imag(over_thresh), 'b.', 'MarkerSize', 7 );
                    
                    % mobiles under the threshold are in red
                    under_thresh = ms_locations{i}( find( SINR_dB < obj.SINR_thresh) );
                    plot( real(under_thresh), MapSize-imag(under_thresh), 'r.', 'MarkerSize', 7  );
                    
                    %title( sprintf( 'Time = %d:00\n', i ) );                
                    text( 600, 30, sprintf( 'Time = %d:00\n', i ) );
                    
                    % capture the frame
                    Recording(i) = getframe;
                    
                    hold off
                end
                
            end
            
        end
        
        function Recording = MapMobiles( obj, ms_locations, varargin )
            % method MapMobiles shows the location of the mobiles
            % and produces a movie
            % 
            % usage: Recording obj.MapMobiles( ms_locations, [pl_exponent] )
            %        Recording = a Matlab movie showing the locations of the mobiles
            %        ms_locations = mobile locations, in a vector or cell array
            %        [pl_exponent] = path loss exponent (optional argument)    
           

            % determine the number of base stations
            [num_bs,num_sec] = size( obj.bs_channels );
            
            % if the pl_exponent is specified, then use it
            if (nargin>2)
                obj.pl_exponent = varargin{1};
            end
            
            
            % if the input ms_locations is not already a cell array, make it one
            if ~iscell( ms_locations )
                ms_in = ms_locations;
                ms_locations = cell(1);
                ms_locations{1} = ms_in;
            end
            
            % open new figure window
            figure
            
            % supress new figure windows
            obj.new_fig = 0;
            
            for i=1:length( ms_locations )
                
                % plot a coverage map
                % Hard coded bounding box
                obj.MapCoverage( obj.pl_exponent, [1 601 1 601] );
                hold on;
                
                % compute the SINR for this realization of the mobile locations
                SINR_dB = obj.ComputeSINR( ms_locations{i} );
                
                % mobiles over the threshold are indicated in blue
                over_thresh = ms_locations{i}( find( SINR_dB >= obj.SINR_thresh) );
                plot( real( over_thresh), imag(over_thresh), ...
                    'k.', 'MarkerSize', 7, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'LineWidth', 1.5 );
                
                % mobiles under the threshold are in red
                under_thresh = ms_locations{i}( find( SINR_dB < obj.SINR_thresh) );
                plot( real(under_thresh), imag(under_thresh), ...
                    'kx', 'MarkerSize', 8, 'LineWidth', 2 );
                
                % title( sprintf( 'Time = %d:00\n', i ) );
                text( 10, 725, sprintf( 'Time = %d:00\n', i ) );
                
                % capture the frame
                Recording(i) = getframe;
                
                hold off
                
            end
            
            % allow new figure windows
            obj.new_fig = 1;
        end
        
        
        function SetPower( obj, bs_power_dBm_input   )
            % method SetPower sets the power of each sector
            % usage obj.SetPower( bs_power_dBm )
            %       bs_power_dBm = a num_sectors by num_bs matrix containing the power of the sectors (in dBm)

            obj.bs_power_dBm = bs_power_dBm_input;            

        end
        
        function SetChannels( obj, bs_channels_input   )
            % method SetChannels sets the channel set assigned to each sector
            % usage obj.SetChannels( bs_channels )
            %       bs_channels = a num_sectors by num_bs matrix containing the channel set used by each sector 
            
            obj.bs_channels = bs_channels_input;            

        end     
        
        function AddBase( obj, location, channels, power_dBm   )
            % method AddBase adds one (or more) base station(s)
            % usage obj.AddBase( location, channels, power_dBm   )
            %       location = location(s) of the new base station(s)
            %       channels = channel set used by the sectors of the new base station(s)
            %       power_dBm = power used by the sector(s) of the new base stations(s), in dBm            
            
            obj.bs_locations = [obj.bs_locations;
                location];
            obj.bs_channels = [obj.bs_channels;
                channels];
            obj.bs_power_dBm = [obj.bs_power_dBm;
                power_dBm];  
        end
        
        function RemoveBase( obj, bs_to_remove )
            % method RemoveBase removes one (or more) base station(s)
            % usage obj.RemoveBase( bs_to_remove )
            %       bs_to_remove = list of the indices of the base stations to remove
            
            % determine the number of base stations
            [num_bs,num_sec] = size( obj.bs_channels ); 
            
            % determine which base stations to keep
            bs_set = 1:num_bs;
            bs_to_keep = setdiff( bs_set, bs_to_remove )
            
            % keep the base stations, the channels, and the powers
            obj.bs_locations = obj.bs_locations( bs_to_keep, : );
            obj.bs_power_dBm = obj.bs_power_dBm( bs_to_keep, : );
            obj.bs_channels = obj.bs_channels( bs_to_keep, : );            
            
        end
        
        function Save( obj, varargin )
            % method Save saves the network design to a file
            % usage obj.Save( [save_filename] )
            %       [save_filename] = name of the file (in single quotes)
            
            if nargin > 1
                obj.save_filename = varargin{1};
            end
            
            save_struct.bs_locations = obj.bs_locations;
            save_struct.bs_channels = obj.bs_channels;
            save_struct.bs_power_dBm = obj.bs_power_dBm;
            save_struct.channel_sets = obj.channel_sets;
            save_struct.sector_offset_angle = obj.sector_offset_angle;
            save_struct.SINR_thresh = obj.SINR_thresh;
            save_struct.ms_per_channel = obj.ms_per_channel;
            
            save( obj.save_filename, 'save_struct' )
            
        end
        
        function Load( obj, varargin )
            % method Load loads the network design from a file
            % usage obj.Load( [save_filename] )
            %       [save_filename] = name of the file (in single quotes)
            
            if nargin > 1
                obj.save_filename = varargin{1};
            end
            
            load( obj.save_filename, 'save_struct' )
            
            obj.bs_locations = save_struct.bs_locations;
            obj.bs_channels = save_struct.bs_channels;
            obj.bs_power_dBm = save_struct.bs_power_dBm;
            obj.channel_sets = save_struct.channel_sets;
            obj.sector_offset_angle = save_struct.sector_offset_angle;
            obj.SINR_thresh = save_struct.SINR_thresh;
            obj.ms_per_channel = save_struct.ms_per_channel;
            
        end      
                                   
    end
    
end

