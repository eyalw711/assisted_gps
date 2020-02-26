function snapshot_gps_main()
clear all
% rng(711);

gs = snapshot_gps_stats();
res = snapshot_gps_nav(gs);

res
end

