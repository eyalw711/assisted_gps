function gs = snapshot_gps_stats()

gs = struct();

gs.tm = datestr(now,'dd_mm_yyyy__HH_MM_SS');
gs.OBS_TYPE = 'C1';
gs.scenario_id = 1;
gs = snapshot_gps_input_dB(gs);

gs.sigmaCode    = 10e-9;

%%% assistance %%%
gs.ellBar = lla2ecef(gs.gt_lla)'; % ECEF

ellErrMag = 10; % 5e3
gs.ellBar = gs.ellBar + ellErrMag*randn(3,1); % naive 3d error
gs.bBar = 0;

end