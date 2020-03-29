function [x,stats] = ilsSearch(R, b, u, options)

[m,n] = size(R);

stats.iterations = 0;
stats.candidates  = 0;
stats.outcome    = 'running';

if nargin < 3
    u = inf;
end

if nargin < 4
    options.maxIter = 100*n^3;
end
    
%x = round( R\b )

%u = (R*x-b)' * (R*x - b)

%xCandidate = zeros(n,1);
%lb      = zeros(n,1);
%ub      = zeros(n,1);
rounded = zeros(n,1);
offset  = zeros(n,1);
xLevel  = zeros(n,1);
dir     = zeros(n,1);

%[minR,minRind] = min(abs(diag(R)));
%fprintf('min abs(diag(R))=%e at %d\n',minR,minRind);

xReal = b(n) / R(n,n);
rounded(n) = round(xReal);
offset(n) = 0;
if rounded(n) >= xReal
    dir(n) = -1;
else
    dir(n) =  1;
end
stats.iterations = stats.iterations + 1;

level = n;

while true

    stats.iterations = stats.iterations + 1;
    if mod(stats.iterations,10000)==0
        fprintf('ilsSearch iteration %d bound=%e cadidates=%d xn=%d level=%d\n',stats.iterations,u,stats.candidates,xLevel(n),level);
    end
    
    if stats.iterations > options.maxIter
        outcome = 'too many iterations';
        break;
    end
    
    xLevel(level) = rounded(level) + offset(level);
    %fprintf('set x(%d) <- %f rounded=%d offset=%d dir=%d\n',level,xLevel(level),rounded(level),offset(level),dir(level));
    
    res = R(level:n,level:n)*xLevel(level:n)-b(level:n);
    
    if res' * res >= u
        %fprintf('popping from level %d, ssq=%e u=%e\n',level,res'*res,u);
        level = level + 1; % pop
        if level == n+1
            stats.outcome = 'done';
            break;
        end
        %pause
        continue;
    end
    
    % if we got here, xLevel(level) is admissible
    
    if level == 1 % we are at a leaf, found an integer candidate for x with a better residual
        if res' * res < u
            x = xLevel;
            u = res' * res;
            stats.candidates = stats.candidates + 1;
            %fprintf('its=%d lowering u to %e\n',iterations,u);
        else
            %fprintf('its=%d found a non-optimal assignment, ssq=%e\n',iterations,res'*res);
        end
        %pause
        continue;
    end
    
    % before doing to next level, set up offset for next iteration at this
    % level
    
    %if level==14; fprintf('level=%d dir=%d offset=%d\n',level,dir(level),offset(level));end;
    if dir(level)==-1
        if offset(level) >= 0
           offset(level) = -offset(level) - 1;
        else
           offset(level) = -offset(level);
        end
    else
        if offset(level) <= 0
           offset(level) = -offset(level) + 1;
        else
           offset(level) = -offset(level);
        end
    end    
    %if level==14; fprintf('next offset at level %d to %d\n',level,offset(level));end
    
    % now drop a level
    
    level = level-1;
    
    xReal = (b(level) - R(level,level+1:n)*xLevel(level+1:n)) / R(level,level);
    rounded(level) = round(xReal);
    
    offset(level) = 0;
    if rounded(level) >= xReal
      dir(level) = -1;
    else
      dir(level) =  1;
    end

    %if level==14; fprintf('going down to level %d rounded %d offset %d dir %d\n',level,rounded(level),offset(level),dir(level));end
    %pause
end



