function [Obs_types,approx_loc, time_first_obs] = my_anheader(file)
%MY_ANHEADER Analyzes the header of a RINEX file and outputs
%	       the list of observation types and antenna offset.
%	       End of file is flagged 1, else 0. Likewise for the types.
%	       Typical call: anheader('pta.96o')

fid = fopen(file,'rt');
Obs_types = [];
approx_loc = zeros(3,1);
time_first_obs = zeros(2,1);

while 1			   % Gobbling the header
    line = fgetl(fid);
    answer = findstr(line,'END OF HEADER');
    if  ~isempty(answer)
        break;
    end
    if (line == -1)
        break;
    end
    answer = findstr(line,'# / TYPES OF OBSERV');
    if ~isempty(answer)
        [NObs, line] = strtok(line);
        NoObs = str2num(NObs);
        for k = 1:NoObs
            [ot, line] = strtok(line);
            Obs_types = [Obs_types ot];
        end
        ifound_types = 1;
    end
    answer = findstr(line,'APPROX POSITION XYZ');
    if ~isempty(answer)
        [approx_x, line] = strtok(line);
        approx_x = str2num(approx_x);
        [approx_y, line] = strtok(line);
        approx_y = str2num(approx_y);
        [approx_z, line] = strtok(line);
        approx_z = str2num(approx_z);
        approx_loc = [approx_x; approx_y; approx_z];
    end
    answer = findstr(line,'TIME OF FIRST OBS');
    if ~isempty(answer)
        [year, line] = strtok(line);
        year = str2num(year);
        [month, line] = strtok(line);
        month = str2num(month);
        [day, line] = strtok(line);
        day = str2num(day);
        [hour, line] = strtok(line);
        hour = str2num(hour);
        [minute, line] = strtok(line);
        minute = str2num(minute);
        [second, line] = strtok(line);
        second = str2num(second);
        jd = julday(year,month,day,hour+minute/60+second/3600);
        [week,sow] = gps_time(jd);
        time_first_obs = [week; sow];
    end
end
fclose(fid);
end
