%clear all;
rng(777);
NEPOCHS = 300;
% file_idx = 2;
SAVE_FIGS = 0;

tm = datestr(now,'dd_mm_yyyy__HH_MM_SS');
if SAVE_FIGS mkdir('results', tm), end;

N_trials = 1;

scenario_n_files = {'SITE247J.01N', ...
    'rinex_download\abpo_1jan2015_0000\abpo001a00.15n',...
    'rinex_download\amc\amc2001a00.15n', ...
    'rinex_download\crao\crao001a00.15n',...
    'rinex_download\fair\fair001a00.15n'};
scenario_o_files = {'SITE247J.01O',...
    'rinex_download\abpo_1jan2015_0000\abpo001a00.15o',...
    'rinex_download\amc\amc2001a00.15o', ...
    'rinex_download\crao\crao001a00.15o',...
    'rinex_download\fair\fair001a00.15o'};
gt_per_file_lla = ...
    [57.015501  9.987793      48.8;        %[56.997739, 9.993229, 48.8;
    -19.0183060	47.2292117	  1552.9923;
    38.8031250	-104.5245972  1911.3941;
    44.4132611	33.9909833	  365.8;
    64.9779983	-147.4992396  319.1771];

N_files = size(gt_per_file_lla,1);
for file_idx = 4% 1:N_files
    fprintf('file %d/%d\n', file_idx, N_files);
    % Read RINEX ephemerides file and convert to internal Matlab format
    rinexe(scenario_n_files{file_idx},'eph.dat');
    Eph = get_eph('eph.dat');
    
    OBS_TYPE1 = 'P1';
    OBS_TYPE2 = 'P2';
    OBS_TYPE11 = 'L1';
    OBS_TYPE21 = 'L2';
    
    
    light_ms = 299792458 * 0.001; % meters
    TOW_MAX = 604800;
    
    % We identify the observation file and open it
    ofile1 = scenario_o_files{file_idx};
    fid1 = fopen(ofile1,'rt');
    [Obs_types1, rec_loc_assist, ~] = my_anheader(ofile1);
    
    gt_lla = gt_per_file_lla(file_idx, :);
    gt_ecef = lla2ecef(gt_lla);
    
    % add position error:
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
    Pos0 = [];
    
    random_bias = (2*rand(N_trials,1)-1)*1000; %in ms
    
    % Generate assistance time:
    two_sec_worth_of_ms = 2*1000;
    
    time_error_ms = randi([-two_sec_worth_of_ms, two_sec_worth_of_ms], N_trials, 1);
    
    tion = [];
    
    q = 1; eof1 = 0;
    while q < NEPOCHS
        [time1, dt1, sats1, eof1] = fepoch_0(fid1); %fetches epoch data with flag 0
        if eof1
            break;
        end
        
        NoSv1 = size(sats1,1);
        % We pick the observed C1 pseudoranges
        obs1 = grabdata(fid1, NoSv1, NoObs_types1);
        if sum(sum(isnan(obs1),2),1) > 0
            continue;
        end
        i1 = fobs_typ(Obs_types1, OBS_TYPE1);
        i2 = fobs_typ(Obs_types1, OBS_TYPE2);
        
        i11 = fobs_typ(Obs_types1, OBS_TYPE11);
        i21 = fobs_typ(Obs_types1, OBS_TYPE21);
        
        tpos = zeros(5, N_trials); %state vectors for t trials (Error experiments)
        tpos0 = zeros(5, N_trials);

%         tms_positions = zeros(3, N_trials);
        
        for t = 1:N_trials
            rb = random_bias(t);
            te_ms = time_error_ms(t);
            loc_assist = assistance_ecef(t,:)';
            
            orig_PRs1 = obs1(:,i1) / light_ms + rb; % observations Pseudoranges (milliseconds)
            orig_PRs2 = obs1(:,i2) / light_ms + rb;

            % Construct synthetic fractional PRs - P1
            ctn_Ns1 = floor(orig_PRs1); % full ms
            obs_frac_PRs1 = mod(orig_PRs1, 1); %sub-ms
            
            % Construct synthetic fractional PRs - P1
            ctn_Ns2 = floor(orig_PRs2); % full ms
            obs_frac_PRs2 = mod(orig_PRs2, 1); %sub-ms
            
            % L1/L2 phase obs:
            L1phase = obs1(:,i11);
            L2phase = obs1(:,i21);
            
            TOW_assist_ms = floor(time1*1000 + te_ms); % some integer millisecond
            
            curr_obs = [obs_frac_PRs1 obs_frac_PRs2 L1phase L2phase];
            
            mat_Pion_mean = load(sprintf('file%d_Pion_mean.mat', file_idx));
            Pion_mean = mat_Pion_mean.Pion_mean;
            
            [pos, dion] = my_recPosCoarseTime_ls_dualF(curr_obs , sats1, Eph, TOW_assist_ms, loc_assist, 1, Pion_mean);
            tpos(:, t) = pos;
            tion = [tion dion];
%             tms_positions(:, t) = int_ms_pos;
%             Pos = [Pos pos];

            pos = my_recPosCoarseTime_ls(obs_frac_PRs1, sats1, Eph, TOW_assist_ms, loc_assist);
            tpos0(:, t) = pos;
%             PosDgln = [PosDgln pos];
        end
        Pos(q, 1:5, 1:N_trials) = tpos;
        Pos0(q, 1:5, 1:N_trials) = tpos0;
        
        q = q+1;
        if mod(q,10) == 0
            fprintf('epoch %d...\n', q);
        end
    end
    q = q-1;
    fclose(fid1);
    
    %Location Errors Distribution
    mlocs = Pos(:,1:3,:);
    mlocs = permute(mlocs, [2, 1, 3]); %loc - epoch - trial
    mlocs = reshape(mlocs, 3, []);
    mlocs_errs = sqrt(sum((mlocs - repmat(gt_ecef', 1, q*N_trials)).^2, 1));
    [f_mlocs_errs, x_mlocs_errs] = ksdensity(mlocs_errs);
    
    dlocs = Pos0(:,1:3,:);
    dlocs = permute(dlocs, [2, 1, 3]); %loc - epoch - trial
    dlocs = reshape(dlocs, 3, []);
    dlocs_errs = sqrt(sum((dlocs - repmat(gt_ecef', 1, q*N_trials)).^2, 1));
    [f_dlocs_errs, x_dlocs_errs] = ksdensity(dlocs_errs);
    
    figure;
    plot(x_mlocs_errs, f_mlocs_errs, x_dlocs_errs, f_dlocs_errs);
    legend('dual freq' ,'ref');
    xlabel('pos error [m]'); ylabel('pdf');
    title(sprintf('file %d - distribution of position errors', file_idx));
    if SAVE_FIGS savefig(sprintf('results\\%s\\pos_pdf_file_%d.fig', tm, file_idx)), end;
    
    %Assistance-Time Errors Distribution
    mtimes = Pos(:,5,:);
    mtimes = permute(mtimes, [3, 1, 2]); % trial-epoch-time_err_est
    mtimes_errs = mtimes - repmat(time_error_ms*1e-3, 1, q);
    [f_mtimes_errs, x_mtimes_errs] = ksdensity(mtimes_errs(:));
    
    dtimes = Pos0(:,5,:);
    dtimes = permute(dtimes, [3, 1, 2]); % trial-epoch-time_err_est
    dtimes_errs = dtimes - repmat(time_error_ms*1e-3, 1, q);
    [f_dtimes_errs, x_dtimes_errs] = ksdensity(dtimes_errs(:));
    
    figure;
    plot(x_mtimes_errs, f_mtimes_errs, x_dtimes_errs, f_dtimes_errs);
    legend('dual freq' ,'ref');
    xlabel('Assistance-Time Error [sec]'); ylabel('pdf');
    title(sprintf('file %d - distribution of assistance time errors', file_idx));
    if SAVE_FIGS savefig(sprintf('results\\%s\\time_pdf_file_%d.fig', tm, file_idx)), end;

    drawnow;
end
%
% print -deps easy3
%%%%%%%%%%%%%%%%%%%%% end easy3.m %%%%%%%%%%%%%%%





