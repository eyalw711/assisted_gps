function [state, int_ms_pos] = my_recPosCoarseTime_ls(... % old return values are [pos, El, GDOP, basic_obs]
    obs, sats, Eph, TOW_assist_ms, rec_loc_assist, is_dgln, is_int_ms)
% MY_RECPOSCOARSETIME_LS Computation of receiver position from fractional
%          pseudoranges using coarse time navigation and least squares

% obs (Observations) are the fractional PRs (sub millisecond)
% sats are the SV numbers associated with each observation
% Eph is a table of ephemerides, each column associated with a satellite

v_light = 299792458;
dtr = pi/180;
numSVs = size(obs,1);  % number of svs
el = zeros(numSVs,1);

TOW_assist = TOW_assist_ms * 1e-3;

% identify ephemerides columns in Eph
for k = 1:numSVs
    col_Eph(k) = find_eph(Eph, sats(k), TOW_assist);
end
Eph = Eph(:, col_Eph); % sort according to sats argument
clear col_Eph; %shouldn't use it anymore after sorting

% preliminary guess for receiver position, common bias, and assistance
% error
state = [rec_loc_assist; 0; 0]; % [x y z b et]'
no_iterations = 20; 
ps_corr = [];
sat_pos = [];

T_tilde = TOW_assist;

% find stellite positions at T_tilde
satPos_at_T_tilde = zeros(3,numSVs);
for k = 1:numSVs
    tx_GPS = tx_RAW2tx_GPS(T_tilde, Eph(:, k));
    satPos_at_T_tilde(:,k) = satpos(tx_GPS, Eph(:,k));
end

% and then find closest one (later on figure out how to find highest)
approx_distances = sqrt(sum((repmat(rec_loc_assist, 1, numSVs) - satPos_at_T_tilde).^2, 1));
[~, svInxListByDistance] = sort(approx_distances);

% assign N numbers:
[Ns, N0_inx] = my_assignNs(sats, svInxListByDistance, obs, Eph, TOW_assist ,rec_loc_assist, approx_distances);
Ns = Ns + 711;

% now find K numbers:
Ks = arrayfun(@(x) x - Ns(N0_inx), Ns);
if exist('is_dgln', 'var') && is_dgln 
    Ks = zeros(size(Ks));
end

fullPRs = Ns + obs; % full pseudoranges reconstruction in ms
fullPRs = fullPRs * (v_light * 1e-3); % in meters

for iter = 1:no_iterations
    H = [];
    delta_z = []; % observed minus computed observation
    
    b = state(4); % in meters
    Et = state(5); % in seconds
    
    for k = svInxListByDistance
        Kk = Ks(k);
        
%         [tx_GPS, tcorr] = tx_RAW2tx_GPS(T_tilde, Eph(:, k));
        [tx_GPS, tcorr] = tx_RAW2tx_GPS(T_tilde - Et - Kk*1e-3, Eph(:, k));
%         tx_RAW = T_tilde - Et - Kk*1e-3;
%         [tx_GPS, tcorr] = tx_RAW2tx_GPS(tx_RAW, Eph(:, k));
        
        X = satpos(tx_GPS, Eph(:,k));
        X_fut = satpos(tx_GPS + 1, Eph(:,k));
        
        % this if case calculates something about trop, leave it like this
        if iter == 1
            traveltime = 0.072;
            Rot_X = X;
            Rot_X_fut = X_fut; 
            trop = 0;
        else
            rho2 = (X(1)-state(1))^2+(X(2)-state(2))^2+(X(3)-state(3))^2; %distance squared
            traveltime = sqrt(rho2)/v_light;
            Rot_X = e_r_corr(traveltime,X);
            
            rho2_fut = (X_fut - state(1:3)).^2;
            traveltime_fut = sqrt(rho2_fut)/v_light;
            Rot_X_fut = e_r_corr(traveltime,X_fut);
            
            rho2 = (Rot_X(1)-state(1))^2+(Rot_X(2)-state(2))^2+(Rot_X(3)-state(3))^2;          
            [az,el,dist] = topocent(state(1:3,:),Rot_X-state(1:3,:));                                                            
            if iter == no_iterations, El(k) = el; end
            trop = tropo(sin(el*dtr),0.0,1013.0,293.0,50.0,...
                0.0,0.0,0.0);    
        end
        
        % subtraction of pos(4) corrects for receiver clock offset and
        % v_light*tcorr is the satellite clock offset
        if iter == no_iterations
            ps_corr = [ps_corr; obs(k)+v_light*tcorr-trop];
            sat_pos = [sat_pos; X'];
        end
        
        %predictedPR = norm(Rot_X-state(1:3)) / (v_light * 0.001) - b*1e3 + tcorr*1e3 - trop / (v_light * 1e-3); % ms
        predictedPR = norm(Rot_X-state(1:3)) + b - tcorr*v_light + trop; % meters
        delta_z = [delta_z; fullPRs(k) - predictedPR];  % meters
        
        % now add row to matrix H according to:
        % -e_k 1 v_k
        % notice that it is easier to plug in the location of the satellite
        % at its T_dot estimation, i.e. Rot_X
        sat_vel_mps = Rot_X_fut - Rot_X;
        
        e_k = (Rot_X - (Et+Kk*1e-3)*sat_vel_mps - state(1:3));
        e_k = e_k/norm(e_k);   %    e_k = e_k / fullPRs(k);
        
        v_k = -sat_vel_mps'*e_k; %relative speed in mps
        H_row = [-e_k' 1 v_k];
        H = [H; H_row];
    end % i
    x = H\delta_z;
    state = state+x;

%     dzmag = [dzmag norm(delta_z)];
%     if iter == no_iterations,
%         % GDOP = sqrt(trace(inv(H'*H))); 
%         % two lines that solve an exercise on computing tdop
%         % invm = inv(A'*A);
%         % tdop = sqrt(invm(4,4))
%     end
end % iter

if exist('is_int_ms', 'var') && is_int_ms
    %construct time is seconds!
    assistance_time_err = state(5);
    concluded_time_sat0 = TOW_assist - assistance_time_err;
    sat0_hypotheses_sec = ((round(concluded_time_sat0*1000)-20):1:(round(concluded_time_sat0*1000)+20)) * 1e-3;
    rec_time_hypotheses = sat0_hypotheses_sec + fullPRs(N0_inx)/ v_light;
    int_ms_positions = zeros(4, length(rec_time_hypotheses));
    int_ms_resid_mags = zeros(size(rec_time_hypotheses));
    
    for iter = 1:length(rec_time_hypotheses)
        time = rec_time_hypotheses(iter);
        [pos, ~, ~, ~ , resid_mag] = recpo_ls(fullPRs,sats,time,Eph);
        int_ms_positions(:, iter) = pos;
        int_ms_resid_mags(iter) = resid_mag;
    end
    
    % choose least residual:
    [~, Imin] = min(int_ms_resid_mags);
    int_ms_pos = int_ms_positions(1:3, Imin);
    
    
end

% basic_obs = [sat_pos ps_corr];
% figure; plot(dzmag); xlabel('iterations'); ylabel('norm(dz)');
%%%%%%%%%%%%%%%%%%%%%  recpo_ls.m  %%%%%%%%%%%%%%%%%%%%%
end
