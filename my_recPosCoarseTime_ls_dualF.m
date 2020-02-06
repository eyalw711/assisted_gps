function [state, dPion] = my_recPosCoarseTime_ls_dualF(... % old return values are [pos, El, GDOP, basic_obs]
    obs, sats, Eph, TOW_assist_ms, rec_loc_assist, is_dualF, Pion_mean)
% MY_RECPOSCOARSETIME_LS Computation of receiver position from fractional
%          pseudoranges using coarse time navigation and least squares

% obs (Observations) are the fractional PRs (sub millisecond)
% sats are the SV numbers associated with each observation
% Eph is a table of ephemerides, each column associated with a satellite

obs2 = obs(:, 2);
phi1 = obs(:, 3);
phi2 = obs(:, 4);

obs = obs(:,1);

f1 = 1575.42e6;
f2 = 1227.60e6;

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
[Ns2, ~] = my_assignNs(sats, svInxListByDistance, obs2, Eph, TOW_assist ,rec_loc_assist, approx_distances);

% now find K numbers:
Ks = arrayfun(@(x) x - Ns(N0_inx), Ns);

fullPRs = Ns + obs; % full pseudoranges reconstruction in ms
fullPRs = fullPRs * (v_light * 1e-3); % in meters

fullPRs2 = Ns2 + obs2; % full pseudoranges reconstruction in ms
fullPRs2 = fullPRs2 * (v_light * 1e-3); % in meters

dPion = zeros(size(sats));

alpha = (f1/f2)^2;
A = [1 1 0 0;
    1 -1 v_light/f1 0;
    1 alpha 0 0;
    1 -alpha 0 v_light/f2];
invA = inv(A);

if exist('is_dualF','var') && is_dualF
    for kk = 1:length(sats)
        xf = invA*[fullPRs(kk); phi1(kk); fullPRs2(kk); phi2(kk)];
        dPion(kk) = (-v_light/f1*xf(3)+v_light/f2*xf(4) + phi1(kk)-phi2(kk))/(alpha-1);
    end
end
% dPion = ((f2^2) / (f2^2 - f1^2)) * (fullPRs - fullPRs2);

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
        predictedPR = norm(Rot_X-state(1:3)) + b - tcorr*v_light + trop + Pion_mean(k); % Pion_mean(k); %dPion(k); % meters
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


% basic_obs = [sat_pos ps_corr];
% figure; plot(dzmag); xlabel('iterations'); ylabel('norm(dz)');
%%%%%%%%%%%%%%%%%%%%%  recpo_ls.m  %%%%%%%%%%%%%%%%%%%%%
end
