function epochs = ubx_reader(file)
fid = fopen(file,'rt');

epochs = [];
in_epoch = false;
obs_inx = 1;
while 1 %epochs loop
    line = fgetl(fid);
    if (line == -1)
        break;
    end
    
    if ~in_epoch
        answer = findstr(line,'recvTOW');
        
        if  ~isempty(answer)
            curr_epoch = struct();
            [tok, rem] = strtok(line);
            [tok, rem] = strtok(rem);
            curr_epoch.recvTOW = str2double(strip(tok));
            curr_epoch.obs = struct();
            curr_epoch.obs.sv = [];
            curr_epoch.obs.pr = [];
            in_epoch = true;
            continue;
        end
        
    else
        
        answer = findstr(line,'UBX-RXM-RAWX');
        if  ~isempty(answer)
            epochs = [epochs; curr_epoch];
            in_epoch = false;
            obs_inx = 1;
            continue;
        end
        
        answer = findstr(line,'prMes');
        if  ~isempty(answer)
            [tok, rem] = strtok(line);
            [tok, rem] = strtok(rem);
            latest_prMes = str2double(strip(tok));
            continue;
        end
        
        answer = findstr(line,'gnssId');
        if  ~isempty(answer)
            [tok, rem] = strtok(line);
            [tok, rem] = strtok(rem);
            latest_gnssId = str2num(strip(tok));
            continue;
        end
        
        answer = findstr(line,'svId');
        if  ~isempty(answer)
            [tok, rem] = strtok(line);
            [tok, rem] = strtok(rem);
            latest_svId = str2num(strip(tok));
            continue;
        end
        
        answer = findstr(line,'sigId');
        if  ~isempty(answer)
            [tok, rem] = strtok(line);
            [tok, rem] = strtok(rem);
            latest_sigId = str2num(strip(tok));
            
            if latest_gnssId == 0  && latest_sigId == 0
                curr_epoch.obs.sv(obs_inx) = latest_svId;
                curr_epoch.obs.pr(obs_inx) = latest_prMes;
                obs_inx = obs_inx + 1;
            end
            
            continue;
        end
    end
end

fclose(fid);
end

