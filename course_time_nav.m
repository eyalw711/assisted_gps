%clear all;
% rng(777);
% file_idx = 2;
SAVE_FIGS = 0;

d=3;
bBar = 0;
tcode = 1e-3;
c = physconst('LightSpeed');
d_coarse = 76.5e-3;

tm = datestr(now,'dd_mm_yyyy__HH_MM_SS');
if SAVE_FIGS mkdir('results', tm), end;

N_trials = 1;

scenario_n_files = {'SITE247J.01N', ...
    'rinex_download/abpo_1jan2015_0000/abpo001a00.15n',...
    'rinex_download/amc/amc2001a00.15n', ...
    'rinex_download/crao/crao001a00.15n',...
    'rinex_download/fair/fair001a00.15n'};
scenario_o_files = {'SITE247J.01O',...
    'rinex_download/abpo_1jan2015_0000/abpo001a00.15o',...
    'rinex_download/amc/amc2001a00.15o', ...
    'rinex_download/crao/crao001a00.15o',...
    'rinex_download/fair/fair001a00.15o'};
gt_per_file_lla = ...
    [57.015501  9.987793      48.8;        %[56.997739, 9.993229, 48.8;
    -19.0183060	47.2292117	  1552.9923;
    38.8031250	-104.5245972  1911.3941;
    44.4132611	33.9909833	  365.8;
    64.9779983	-147.4992396  319.1771];

N_files = size(gt_per_file_lla,1);
for file_idx = 5%1:N_files
    fprintf('file %d/%d\n', file_idx, N_files);
    
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
    gt_ecef = lla2ecef(gt_lla)';
    
%     % add position error:
%     az = rand(N_trials,1)*360;
%     dist = rand(N_trials,1)*3000; %meters
%     
%     distUnits = 'm';
%     % Convert input distance to earth degrees (Lat, Lon are typicaly given in degrees)
%     arclen = rad2deg(dist/earthRadius(distUnits));
%     [latOut,lonOut] = reckon(gt_lla(1), gt_lla(2), arclen, az);
    
    assistance_lla = repmat(gt_lla, N_trials,1); % + 1e-4*randn(N_trials, 3); %[latOut,lonOut, gt_lla(3)*ones(N_trials,1)];
    assistance_ecef = lla2ecef(assistance_lla) + randn(N_trials, 3);
    ellBar = assistance_ecef(1,:)';
            
%     fprintf('gt_lla: %f, %f, %f\n', gt_lla(1), gt_lla(2), gt_lla(3));
%     fprintf('loc assistance err = %f [m], %f deg\n', dist, az);
%     fprintf('assistance_lla: %f, %f, %f\n', assistance_lla(1), assistance_lla(2), assistance_lla(3));
    
    NoObs_types1 = size(Obs_types1,2)/2;
    Pos = [];
    PosDgln = [];
    PosMS_UPGRD = [];
    
    clockBiasMag = 120; % seconds
    clock_bias = (2*rand(N_trials,1)-1)*1000*clockBiasMag; %[f] in ms
    ctn_clock_bias = clock_bias(1) * 1e-3; % [sec]
%     ctn_clock_bias = abs(ctn_clock_bias); %NOTICE REMOVE LATER
    
    % Generate assistance time:
    two_sec_worth_of_ms = 2*1000;
    
    time_error_ms = randi([-two_sec_worth_of_ms, two_sec_worth_of_ms], N_trials, 1);
    
    q = 1; eof1 = 0;
    while 1
        [gps_time, dt1, sats1, eof1] = fepoch_0(fid1); %fetches epoch data with flag 0
        if eof1
            break;
        end
        
        NoSv1 = size(sats1,1);
        % We pick the observed C1 pseudoranges
        obs1 = grabdata(fid1, NoSv1, NoObs_types1);
        obs_type_inx = fobs_typ(Obs_types1, OBS_TYPE);
        
        tpos = zeros(5, N_trials); %state vectors for t trials (Error experiments)
        dgln_tpos = zeros(5, N_trials);
        tms_positions = zeros(3, N_trials);
        
        ranges = obs1(:,obs_type_inx); % assume obs1 is range (almost true in these recordings)
        truNs = floor(ranges / (c*tcode));
        gps_codephases = mod(ranges / (c*tcode), 1); % [0,1]
        gps_arrival_times = gps_time + gps_codephases * tcode; % [sec]
        
        fprintf('clock offset by %.3f milliseconds\n', ctn_clock_bias / tcode);
        
        axis_beta_rev = ctn_clock_bias/tcode - floor(ctn_clock_bias/tcode); % [0,1] this worked +- in reversed b
        axis_beta = mod(-ctn_clock_bias, tcode)/tcode;
        
        ctn_codephases_rev = mod(gps_arrival_times - ctn_clock_bias, tcode)/tcode; % [0,1] this worked +- in reversed b
        ctn_codephases = mod(gps_codephases - axis_beta, 1);
        
        presumed_time = ceil((gps_time + ctn_clock_bias)/tcode)*tcode;
        presumed_time_ms = presumed_time / tcode;
        
        [trypos, H, omc, reconsNs] = my_recPosCoarseTime_ls(ctn_codephases, sats1, Eph, ...
            presumed_time_ms, ellBar, 0, 0);
        
        lastAlgoError = norm(gt_ecef - trypos(1:3))
        
        % our new method
        N_sats = numel(sats1);
        presumed_arrival_times = presumed_time + ctn_codephases*tcode;
        tDhat = presumed_arrival_times - d_coarse - bBar;
        
        [distances, J] = model(ellBar, tDhat, sats1, Eph);
        [~,j] = min(abs(J(:,4)));
        
        correction_times = zeros(N_sats, 1);
        delta_t = zeros(N_sats, 1);
        for i=1:N_sats
            eph_col_inx = find_eph( Eph, sats1(i), tDhat(i));
            eph_col = Eph(:, eph_col_inx);
            [tx_GPS, tcorr] = tx_RAW2tx_GPS(tDhat(i), eph_col);
            correction_times(i) = tcorr;
            delta_t(i) = Eph(19, eph_col_inx);
        end
        
        nu = zeros(N_sats, 1);
        nu(j) = ceil( (distances(j)/c - ctn_codephases(j)*tcode - delta_t(j)) / tcode);
        beta  = (nu(j)+ctn_codephases(j))*tcode - distances(j)/c ;
        nu    = round( (distances/c - ctn_codephases*tcode + beta - delta_t) / tcode );

        %if consistent then phi_i*tcode - distances/c + nu_i*tcode + beta =
        % eps
        
        
        
        %sivan code
        %[~,j] = min(abs( J(:,d+1) ));
        %nu(j) = ceil( (distances(j)/c - phi(j)*tcode) / tcode );
        %beta  = (nu(j)+phi(j))*tcode - distances(j)/c
        
        %delta = [ J/c ones(N_sats,1)/c ] \ (ctn_codephases*tcode + nu*tcode - distances/c);
        
        
        [truNs reconsNs  nu]
        [(omc) (((nu+ctn_codephases)*tcode*c - distances + correction_times*c))]
        
        delta = [ J ones(N_sats,1) ] \ ((nu+ctn_codephases)*tcode*c - distances + correction_times*c);
        
        ellHat = ellBar + delta(1:d,1)
        
        bHat = bBar + delta(d+1);
        betaHat = delta(d+2)/c;
        
        ellErr = ellHat - gt_ecef; %target;
        bErr = bHat - ctn_clock_bias;
        
        results.VanDiggelen_locationErr = norm(ellErr);
        results.VanDiggelen_clockErr    = bErr
        fprintf('done iter\n');
        
%         return;
%         %improve:
%         
%         
%         
%         for t = 1:N_trials
%             rb = clock_bias(t); %[f] in ms
%             
%             te_ms = time_error_ms(t);
%             loc_assist = assistance_ecef(t,:)';
%             
%             orig_PRs = obs1(:,obs_type_inx) / light_ms + rb; % observations Pseudoranges (milliseconds)
% 
%             % Construct synthetic fractional PRs:
%             ctn_Ns = floor(orig_PRs); % full ms
%             obs_frac_PRs = mod(orig_PRs, 1); %sub-ms
% 
%             TOW_assist_ms = floor(gps_time*1000 + te_ms); % some integer millisecond
% 
%             pos = my_recPosCoarseTime_ls(obs_frac_PRs, sats1, Eph, TOW_assist_ms, loc_assist, 0, 0);
%             tpos(:, t) = pos;
% %             tms_positions(:, t) = int_ms_pos;
% %             Pos = [Pos pos];
% 
%             pos = my_recPosCoarseTime_ls(obs_frac_PRs, sats1, Eph, TOW_assist_ms, loc_assist, 1);
%             dgln_tpos(:, t) = pos;
% %             PosDgln = [PosDgln pos];
%         end
%         Pos(q, 1:5, 1:N_trials) = tpos;
%         PosDgln(q, 1:5, 1:N_trials) = dgln_tpos;
% %         PosMS_UPGRD(q, 1:3, 1:N_trials) = tms_positions;
        
        q = q+1;
        if mod(q,10) == 0
            fprintf('epoch %d...\n', q);
        end
    end
    q = q-1;
    fclose(fid1);
    
%     %Location Errors Distribution
%     mlocs = Pos(:,1:3,:);
%     mlocs = permute(mlocs, [2, 1, 3]); %loc - epoch - trial
%     mlocs = reshape(mlocs, 3, []);
%     mlocs_errs = sqrt(sum((mlocs - repmat(gt_ecef', 1, q*N_trials)).^2, 1));
%     [f_mlocs_errs, x_mlocs_errs] = ksdensity(mlocs_errs);
%     
%     dlocs = PosDgln(:,1:3,:);
%     dlocs = permute(dlocs, [2, 1, 3]); %loc - epoch - trial
%     dlocs = reshape(dlocs, 3, []);
%     dlocs_errs = sqrt(sum((dlocs - repmat(gt_ecef', 1, q*N_trials)).^2, 1));
%     [f_dlocs_errs, x_dlocs_errs] = ksdensity(dlocs_errs);
%     
%     figure;
%     plot(x_mlocs_errs, f_mlocs_errs, x_dlocs_errs, f_dlocs_errs);
%     legend('my' ,' digglen');
%     xlabel('pos error [m]'); ylabel('pdf');
%     title(sprintf('file %d - distribution of position errors', file_idx));
%     if SAVE_FIGS savefig(sprintf('results\\%s\\pos_pdf_file_%d.fig', tm, file_idx)), end;
%     
%     %Assistance-Time Errors Distribution
%     mtimes = Pos(:,5,:);
%     mtimes = permute(mtimes, [3, 1, 2]); % trial-epoch-time_err_est
%     mtimes_errs = mtimes - repmat(time_error_ms*1e-3, 1, q);
%     [f_mtimes_errs, x_mtimes_errs] = ksdensity(mtimes_errs(:));
%     
%     dtimes = PosDgln(:,5,:);
%     dtimes = permute(dtimes, [3, 1, 2]); % trial-epoch-time_err_est
%     dtimes_errs = dtimes - repmat(time_error_ms*1e-3, 1, q);
%     [f_dtimes_errs, x_dtimes_errs] = ksdensity(dtimes_errs(:));
%     
%     figure;
%     plot(x_mtimes_errs, f_mtimes_errs, x_dtimes_errs, f_dtimes_errs);
%     legend('my' ,' digglen');
%     xlabel('Assistance-Time Error [sec]'); ylabel('pdf');
%     title(sprintf('file %d - distribution of assistance time errors', file_idx));
%     if SAVE_FIGS savefig(sprintf('results\\%s\\time_pdf_file_%d.fig', tm, file_idx)), end;
    
%     %improved int ms version:
%     locms_errs = sqrt(sum((PosMS_UPGRD' - repmat(gt_ecef', 1, q*N_trials)).^2, 1));
%     [f_mslocerr, x_mslocerr] = ksdensity(locms_errs);
%     figure;
%     plot(x_mslocerr, f_mslocerr, x_mlocs_errs, f_mlocs_errs);
%     legend('int ms' ,'LS ms');
%     xlabel('pos err [m]'); ylabel('pdf');
%     title(sprintf('file %d - int ms - distribution of position errors', file_idx));
    
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


function[distances, J] = model(ell, td, sats, Eph)
    N_sats = numel(sats);
    J = zeros(N_sats, 4);
    distances = zeros(N_sats, 1);

    for r = 1:N_sats
        eph_col_inx = find_eph( Eph, sats(r), td(r));
        eph_col = Eph(:, eph_col_inx);

        rho = satpos(td(r), eph_col);
        distances(r) = norm(ell - rho);

        sat_vel_approx = satpos(td(r) + 0.5, eph_col) - ...
            satpos(td(r) - 0.5, eph_col);

        e_r = (ell - rho)';
        e_r = e_r / norm(e_r);
        J(r, 1:3) = e_r;
        J(r, 4)   = e_r * sat_vel_approx;
    end
end



