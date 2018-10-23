% MY EASY3
rng(711); % seed randomness

% Read RINEX ephemerides file and convert to
% internal Matlab format
rinexe('SITE247J.01N','eph.dat');
Eph = get_eph('eph.dat');

OBS_TYPE = 'C1';
light_ms = 299792458 * 0.001; % meters
TOW_MAX = 604800;

% We identify the observation file and open it
ofile1 = 'SITE247J.01O';
fid1 = fopen(ofile1,'rt');
[Obs_types1, rec_loc_assist, ~] = my_anheader(ofile1);
NoObs_types1 = size(Obs_types1,2)/2;
Pos = [];

% There are 20 epochs of data in ofile1: work only on 1st one

[time1, dt1, sats1, eof1] = fepoch_0(fid1);
NoSv1 = size(sats1,1);
obs1 = grabdata(fid1, NoSv1, NoObs_types1);
i = fobs_typ(Obs_types1, OBS_TYPE);

PRs = obs1(:,i); % Pseudoranges from the observations (meters)
PRs = PRs / light_ms; % Pseudoranges in milliseconds

% Now we're going to construct synthetic fractional PRs:
Ns = floor(PRs);
obs_frac_PRs = mod(PRs, 1);

% Generate assistance time:
two_min_worth_of_ms = 2*60*1000;
time_error_ms = randi([-two_min_worth_of_ms, two_min_worth_of_ms]);
TOW_assist = floor(time1*1000 + time_error_ms); % some integer millisecond

% Generate assistance position:
% someplace in Aalborg, 1.38 Km from the receiver: 
% 57.023997, 10.004533
% height 50m
rec_loc_assist = lla2ecef([57.023997 10.004533 50])';

pos = my_recPosCoarseTime_ls(obs_frac_PRs, sats1, Eph, TOW_assist, rec_loc_assist);

%     pos = recpo_ls(obs1(:,i),sats1,time1,Eph);
%     Pos = [Pos pos];
% me = mean(Pos,2);
