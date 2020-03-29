function [Q,R,G,Gi,stats] = ilsReduce(Q,R,options)

[m,n] = size(R);
Gi = eye(n);
G  = eye(n);
%Q  = eye(n);

if nargin < 3
  options.algorithm = 'LLL';
  options.delta     = 1;
end

stats.iterations = 0;

if strcmp(options.algorithm,'LLL')==1
  i = 2;
  while i <= n
    if abs(R(i-1,i)) > abs(R(i-1,i-1))/2
      decrease(i-1,i);
      stats.iterations = stats.iterations + 1;
    end
    
    if R(i,i)^2 < (options.delta - (R(i-1,i)/R(i-1,i-1))^2) * R(i-1,i-1)^2
      swap(i);
      i = max(i-1,2);
      stats.iterations = stats.iterations + 1;
    else
      if true
      for k=i-2:-1:1
        if abs(R(k,k)) < 2*abs(R(k,i))
          decrease(k,i);
          stats.iterations = stats.iterations + 1;
        end
      end
      end
      i = i+1;
    end
  end
  
  return
end


for t=1:100

for j=n-1:-1:1
    for k=j+1:n
        multiplier = round( R(j,k) / R(j,j) );
        if multiplier ~= 0
            R (1:j,k) = R (1:j,k) - multiplier*R (1:j,j);
            G (:,k) = G (:,k) - multiplier*G (:,j);
            Gi(j,:) = Gi(j,:) + multiplier*Gi(k,:);
        end
    end
end

for j=n:-1:2
    if abs(R(j,j)) < abs(R(j-1,j-1))
        %fprintf('swap candidate, j=%d %e %e\n',j, R(j,j), R(j-1,j-1));
        sub = R(j-1:j,j-1:j);
        % check if it is worth it to swap
       [q,r] = qr( sub * [ 0 1 ; 1 0 ]);
       if abs(r(2,2)) > abs(r(1,1))
           %fprintf('swapping, j=%d %e %e\n',j, r(1,1), r(2,2));
           % swap columns in R
           R(:,[j-1 j]) = R(:,[j j-1]); 
           G(:,[j-1 j]) = G(:,[j j-1]);
           Gi([j-1 j],:) = Gi([j j-1],:);
           % transform back into upper triangular form
           R([j-1 j],:) = q'*R([j-1 j],:);
           Q(:,[j-1 j]) = Q(:,[j-1 j])*q;         
           %R
           %pause
       else
           %fprintf('not reducing\n');
       end
    end
end

end

function decrease(ii,jj)
  multiplier = round( R(ii,jj) / R(ii,ii) );
  if multiplier ~= 0
     R (:,jj) = R (:,jj) - multiplier*R (:,ii);
     G (:,jj) = G (:,jj) - multiplier*G (:,ii);
     Gi(ii,:) = Gi(ii,:) + multiplier*Gi(jj,:);
   end 
end

function swap(ii)
  sub = R(ii-1:ii,ii-1:ii);
  [q,r] = qr( sub * [ 0 1 ; 1 0 ]);
  R(:,[ii-1 ii]) = R(:,[ii ii-1]); 
  G(:,[ii-1 ii]) = G(:,[ii ii-1]);
  Gi([ii-1 ii],:) = Gi([ii ii-1],:);
  % transform back into upper triangular form
  R([ii-1 ii],:) = q'*R([ii-1 ii],:);
  Q(:,[ii-1 ii]) = Q(:,[ii-1 ii])*q;         
end



end