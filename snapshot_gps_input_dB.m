function gs = snapshot_gps_input_dB(gs)

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

N = numel(scenario_n_files);
assert( isscalar(gs.scenario_id), 'id should be scalar');
assert( (1 <= gs.scenario_id) && (gs.scenario_id <= N), sprintf('invalid scenario id = %d', gs.scenario_id));

rinexe(scenario_n_files{gs.scenario_id},'eph.dat');
gs.ephemeris = get_eph('eph.dat');

ofile = scenario_o_files{gs.scenario_id};
gs.obs_reader = ObservationReader(ofile, gs);
gs.gt_lla = gt_per_file_lla(gs.scenario_id, :);
gs.gt = lla2ecef(gs.gt_lla)';

end

