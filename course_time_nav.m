%clear all;
% rng(777);
% file_idx = 2;
N_trials = 10;

scenario_n_files = {'SITE247J.01N', ...
    'rinex_download\abpo_1jan2015_0000\abpo001a00.15n',...
    'rinex_download\amc\amc2001a00.15n', ...
    'rinex_download\crao\crao001a00.15n'};
scenario_o_files = {'SITE247J.01O',...
    'rinex_download\abpo_1jan2015_0000\abpo001a00.15o',...
    'rinex_download\amc\amc2001a00.15o', ...
    'rinex_download\crao\crao001a00.15o'};
gt_per_file_lla = ...
    [56.997739, 9.993229, 48.8;
    -19.0183060	47.2292117	1552.9923;
    38.8031250	-104.5245972	1911.3941;
    44.4132611	33.9909833	365.8];

for file_idx = 1:size(gt_per_file_lla,1)
    % Read RINEX ephemerides file and convert to internal Matlab format
    rinexe(scenario_n_files{file_idx},'eph.dat');
    Eph = get_eph('eph.dat');
    
    OBS_TYPE = 'C1'; % Pseudorange using C/A code on L1
    light_ms = 299792458 * 0.001; % meters
    TOW_MAX = 604800;
    
    % We identify the observation file and open it
    ofile1 = scenario_o_files{file_idx};
    fid1 = fopen(ofile1,'rt');
    [Obs_types1, rec_loc_assist, ~] = my_anheader(ofile1);
    
    gt_lla = gt_per_file_lla(file_idx, :);
    gt_ecef = lla2ecef(gt_lla);
    
    % add error:
    az = rand(N_trials,1)*360;
    dist = rand(N_trials,1)*3000; %meters
    
    distUnits = 'm';
    % Convert input distance to earth degrees (Lat, Lon are typicaly given in degrees)
    arclen = rad2deg(dist/earthRadius(distUnits));
    [latOut,lonOut] = reckon(gt_lla(1), gt_lla(2), arclen, az);
    
    assistance_lla = [latOut,lonOut, gt_lla(3)*ones(N_trials,1)];
    assistance_ecef = lla2ecef(assistance_lla);
    
%     fprintf('gt_lla: %f, %f, %f\n', gt_lla(1), gt_lla(2), gt_lla(3));
%     fprintf('loc assistance err = %f [m], %f deg\n', dist, az);
%     fprintf('assistance_lla: %f, %f, %f\n', assistance_lla(1), assistance_lla(2), assistance_lla(3));
    
    NoObs_types1 = size(Obs_types1,2)/2;
    Pos = [];
    PosDgln = [];
    
    random_bias = (2*rand(N_trials,1)-1)*1000; %in ms
    
    % Generate assistance time:
    two_sec_worth_of_ms = 2*1000;
    time_error_ms = randi([-two_sec_worth_of_ms, two_sec_worth_of_ms], N_trials, 1);
    
    q = 1; eof1 = 0;
    while 1
        [time1, dt1, sats1, eof1] = fepoch_0(fid1); %fetches epoch data with flag 0
        if eof1
            break;
        end
        
        NoSv1 = size(sats1,1);
        % We pick the observed C1 pseudoranges
        obs1 = grabdata(fid1, NoSv1, NoObs_types1);
        i = fobs_typ(Obs_types1, OBS_TYPE);
        
        tpos = zeros(5, N_trials);
        dgln_tpos = zeros(5, N_trials);
        
        for t = 1:N_trials
            rb = random_bias(t);
            te_ms = time_error_ms(t);
            loc_assist = assistance_ecef(t,:)';
            
            orig_PRs = obs1(:,i) / light_ms + rb; % observations Pseudoranges (milliseconds)
            % inherent bias is 0.379229019830779 ms

            % Construct synthetic fractional PRs:
            ctn_Ns = floor(orig_PRs); % full ms
            obs_frac_PRs = mod(orig_PRs, 1); %sub-ms

            TOW_assist_ms = floor(time1*1000 + te_ms); % some integer millisecond

            pos = my_recPosCoarseTime_ls(obs_frac_PRs, sats1, Eph, TOW_assist_ms, loc_assist, 0);
            tpos(:, t) = pos;
%             Pos = [Pos pos];

            pos = my_recPosCoarseTime_ls(obs_frac_PRs, sats1, Eph, TOW_assist_ms, loc_assist, 1);
            dgln_tpos(:, t) = pos;
%             PosDgln = [PosDgln pos];
        end
        Pos(q, 1:5, 1:N_trials) = tpos;
        PosDgln(q, 1:5, 1:N_trials) = dgln_tpos;
        q = q+1;
    end
    q = q-1;
    fclose(fid1);
    
    %Location Errors Distribution
    mlocs = Pos(:,1:3,:);
    mlocs = permute(mlocs, [2, 1, 3]); %loc - epoch - trial
    mlocs = reshape(mlocs, 3, []);
    mlocs_errs = sqrt(sum((mlocs - repmat(gt_ecef', 1, q*N_trials)).^2, 1));
    [f_mlocs_errs, x_mlocs_errs] = ksdensity(mlocs_errs);
    
    dlocs = PosDgln(:,1:3,:);
    dlocs = permute(dlocs, [2, 1, 3]); %loc - epoch - trial
    dlocs = reshape(dlocs, 3, []);
    dlocs_errs = sqrt(sum((dlocs - repmat(gt_ecef', 1, q*N_trials)).^2, 1));
    [f_dlocs_errs, x_dlocs_errs] = ksdensity(dlocs_errs);
    
    figure;
    plot(x_mlocs_errs, f_mlocs_errs, x_dlocs_errs, f_dlocs_errs);
    legend('my' ,' digglen');
    xlabel('pos error [m]'); ylabel('pdf');
    title(sprintf('file %d - distribution of position errors', file_idx));
    
    %Assistance-Time Errors Distribution
    mtimes = Pos(:,5,:);
    mtimes = permute(mtimes, [3, 1, 2]); % trial-epoch-time_err_est
    mtimes_errs = mtimes - repmat(time_error_ms*1e-3, 1, q);
    [f_mtimes_errs, x_mtimes_errs] = ksdensity(mtimes_errs(:));
    
    dtimes = PosDgln(:,5,:);
    dtimes = permute(dtimes, [3, 1, 2]); % trial-epoch-time_err_est
    dtimes_errs = dtimes - repmat(time_error_ms*1e-3, 1, q);
    [f_dtimes_errs, x_dtimes_errs] = ksdensity(dtimes_errs(:));
    
    figure;
    plot(x_mtimes_errs, f_mtimes_errs, x_dtimes_errs, f_dtimes_errs);
    legend('my' ,' digglen');
    xlabel('Assistance-Time Error [sec]'); ylabel('pdf');
    title(sprintf('file %d - distribution of assistance time errors', file_idx));
    
%     fprintf('My:\n');
%     me = mean(Pos,2);
%     computed_lla = ecef2lla(me(1:3)');
%     fprintf('computed_lla: %f, %f, %f\n', computed_lla(1), computed_lla(2), computed_lla(3));
%     fprintf('fix distance from assist approx = %f [meters]\n', norm(me(1:3)-assistance_ecef'));
%     fprintf('fix distance from expected = %f [meters]\n', norm(me(1:3)-lla2ecef(gt_lla)'));
%     fprintf('true time error = %f [sec], estimated time error = %f [sec]\n', time_error_ms*1e-3, me(5));
%     fprintf('ignore-> added bias = %f [ms], estimated = %f [ms]\n', random_bias + 0.379229019830779, me(4)/light_ms); %113690 is in the data...
% 
%     figure;
%     plot((Pos(1:3,:)-Pos(1:3,1)*ones(1,q))','linewidth',2)
%     title('Positions Over Time','fontsize',16)
%     legend('X','Y','Z')
%     xlabel('Epochs [1 s interval]','fontsize',16)
%     ylabel('Variation in Coordinates, Relative to the First Epoch [m]','fontsize',16)
%     set(gca,'fontsize',16)
%     
%     fprintf('\nDigglen:\n');
%     me = mean(PosDgln,2);
%     computed_lla = ecef2lla(me(1:3)');
%     fprintf('computed_lla: %f, %f, %f\n', computed_lla(1), computed_lla(2), computed_lla(3));
%     fprintf('fix distance from assist approx = %f [meters]\n', norm(me(1:3)-assistance_ecef'));
%     fprintf('fix distance from expected = %f [meters]\n', norm(me(1:3)-lla2ecef(gt_lla)'));
%     fprintf('true time error = %f [sec], estimated time error = %f [sec]\n', time_error_ms*1e-3, me(5));
%     fprintf('ignore-> added bias = %f [ms], estimated = %f [ms]\n', random_bias + 0.379229019830779, me(4)/light_ms); %113690 is in the data...
%     
%     figure;
%     plot((PosDgln(1:3,:)-PosDgln(1:3,1)*ones(1,q))','linewidth',2)
%     title('Positions Over Time','fontsize',16)
%     legend('X','Y','Z')
%     xlabel('Epochs [1 s interval]','fontsize',16)
%     ylabel('Variation in Coordinates, Relative to the First Epoch [m]','fontsize',16)
%     set(gca,'fontsize',16)
    
    drawnow;
end
%
% print -deps easy3
%%%%%%%%%%%%%%%%%%%%% end easy3.m %%%%%%%%%%%%%%%





