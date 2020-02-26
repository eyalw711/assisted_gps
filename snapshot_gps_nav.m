function results = snapshot_gps_nav(gs)


%%%
%%% SETTINGS
%%%

ONE_SHOT = 1;

%%%
%%% CONSTANTS
%%%

d = 3;
d_coarse = 76.5e-3; %ms
c = physconst('LightSpeed'); %300e6;
tcode = 1e-3;


gs

results = struct();


% unknown
clockBiasMag = 5; %500
clock_bias = clockBiasMag*(rand-0.5);

iter = 1;
while 1
    obs_gps = gs.obs_reader.fetch();
    if obs_gps.eof
        break;
    end
    fprintf('%d: at time %f: %d observations\n', iter, obs_gps.time, numel(obs_gps.obs));
    obs_ctn = obs_gps;
    obs_ctn.obs = mod((obs_gps.obs)/(c*tcode), 1);   % code phase obs
    obs_ctn.time = obs_gps.time + clock_bias;
    
    % knowns / observations
    N_sats = numel(obs_ctn.sats);
    t = obs_ctn.time + tcode*obs_ctn.obs;       % arrival times
    tDhat = t - d_coarse - gs.bBar;             % estimated departure times
    
    %%% van diggelen's method %%%
    %%% shadowing      method %%%
    
    [distances, J] = model(gs.ellBar, tDhat);
    
    [~,j] = min(abs(J(:,4)));
    
    nu = zeros(N_sats, 1);
    nu(j) = ceil( (distances(j)/c - obs_ctn.obs(j)*tcode) / tcode);
    beta  = (nu(j)+obs_ctn.obs(j))*tcode - distances(j)/c;
    
    nu    = round( (distances/c - obs_ctn.obs*tcode + beta) / tcode );
    
    w(1:N_sats,1)   = (1/gs.sigmaCode) * ones(N_sats,1);
    W = diag(w);
    
    delta = (W*[ J/c ones(N_sats,1)/c ]) \ (W*(obs_ctn.obs*tcode + nu*tcode - distances/c));
    
    ellHat = gs.ellBar + delta(1:d,1)
    
    bHat = gs.bBar + delta(d+1);
    betaHat = delta(d+2)/c;
    
    ellErr = ellHat - gs.gt; %target;
    bErr = bHat - clock_bias;
    
    results.VanDiggelen_locationErr = norm(ellErr);
    results.VanDiggelen_clockErr    = bErr;
    
    %%% IMPROVE USING THE NEW ESTIMATES
    
    distances = model(ellHat, tDhat); % just for imroving transmit times
    tDhat = t - distances/c - bHat;   % improved transmit times
    [distances,J] = model(ellHat, tDhat);
    
    results.VanDiggelen_resNorm = norm( obs_ctn.obs*tcode + nu*tcode - betaHat - distances/c )
    
    delta = (W*[ J/c ones(N_sats,1)/c ]) \ (W*[obs_ctn.obs*tcode + nu*tcode - distances/c])
    
    ellHat = ellHat + delta(1:d,1)
    
    bHat = bHat + delta(4);
    betaHat = delta(5)/c;
    
    ellErr = ellHat - gs.gt;
    bErr = bHat - clock_bias;
    
    results.VanDiggelen_locationErr2 = norm(ellErr);
    results.VanDiggelen_clockErr2    = bErr;
    
    distances = model(ellHat, tDhat); % just for imroving transmit times
    results.VanDiggelen_resNorm2 = norm( obs_ctn.obs*tcode + nu*tcode - betaHat - distances/c );
    
    %return
    %%% end of Van Diggelen
    
    if ONE_SHOT
        break;
    else
        iter = iter + 1;
    end
end


    function[distances, J] = model(ell, td)
        J = zeros(N_sats, 4);
        distances = zeros(N_sats, 1);
        
        for r = 1:N_sats
            assert(issorted(obs_ctn.sats), "sats are not sorted here!");
            
            eph_col_inx = find_eph( gs.ephemeris, obs_ctn.sats(r), td(r));
            eph_col = gs.ephemeris(:, eph_col_inx);
            
            rho = satpos(td(r), eph_col);
            distances(r) = norm(ell - rho);
            
            sat_vel_approx = satpos(td(r) + 0.5, eph_col) - ...
                satpos(td(r) - 0.5, eph_col);
            
            e_r = (gs.ellBar - rho)';
            e_r = e_r / norm(e_r);
            J(r, 1:3) = e_r;
            J(r, 4)   = e_r * sat_vel_approx;
        end 
    end
end