function pos = my_recPosCoarseTime_ls(... % old return values are [pos, El, GDOP, basic_obs]
    obs, sats, Eph, TOW_assist, rec_loc_assist)
% MY_RECPOSCOARSETIME_LS Computation of receiver position from fractional
%          pseudoranges using coarse time navigation and least squares

% obs (Observations) are the fractional PRs (sub millisecond)
% sats are the SV numbers associated with each observation
% Eph is a table of ephemerides, each column associated with a satellite

v_light = 299792458;
dtr = pi/180;
numSVs = size(obs,1);  % number of svs
el = zeros(numSVs,1);


% identify ephemerides columns in Eph
for k = 1:numSVs
    col_Eph(k) = find_eph(Eph, sats(k), TOW_assist);
end
Eph = Eph(:, col_Eph); % sort according to sats argument
clear col_Eph; %shouldn't use it anymore after sorting

% preliminary guess for receiver position, common bias, and assistance
% error
state = [rec_loc_assist; 0; 0]; % [x y z b et]'
no_iterations = 6; 
ps_corr = [];
sat_pos = [];

T_tilde = TOW_assist;

% TODO:
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

% now find K numbers:
Ks = arrayfun(@(x) x - Ns(N0_inx), Ns);

fullPRs = Ns + obs; % full pseudoranges reconstruction

for iter = 1:no_iterations
    H = [];
    delta_z = []; % observed minus computed observation
    
    b = state(4);
    Et = state(5);
    
    for k = svInxListByDistance
        Kk = Ks(k);
        
        tx_RAW = T_tilde - Et - Kk;
        [tx_GPS, tcorr] = tx_RAW2tx_GPS(tx_RAW, Eph(:, k));
        
        X = satpos(tx_GPS, Eph(:,k));
        
        % this if case calculates something about trop, leave it like this
        if iter == 1
            traveltime = 0.072;
            Rot_X = X;
            trop = 0;
        else
            rho2 = (X(1)-state(1))^2+(X(2)-state(2))^2+(X(3)-state(3))^2; %distance squared
            traveltime = sqrt(rho2)/v_light;
            Rot_X = e_r_corr(traveltime,X);
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
        
        predictedPR = norm(Rot_X-state(1:3)) / (v_light * 0.001) - b + tcorr*1000 - trop / (v_light * 0.001);
        delta_z = [delta_z; fullPRs(k) - predictedPR]; 
        
        % now add row to matrix H according to:
        % -e_k 1 v_k
        % notice that it is easier to plug in the location of the satellite
        % at its T_dot estimation, i.e. Rot_X
        H_row = [(Rot_X - state(1:3))'/norm(Rot_X - state(1:3)) 1 0]; %TODO: set velocity
        H = [H; H_row];
    end % i
    x = H\delta_z;
    pos = pos+x;
    if iter == no_iterations, GDOP = sqrt(trace(inv(H'*H))); 
        % two lines that solve an exercise on computing tdop
        % invm = inv(A'*A);
        % tdop = sqrt(invm(4,4))
    end
end % iter
basic_obs = [sat_pos ps_corr];

%%%%%%%%%%%%%%%%%%%%%  recpo_ls.m  %%%%%%%%%%%%%%%%%%%%%
