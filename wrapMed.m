function mm = wrapMed( ss )

[N, m] = size(ss);

ss = mod(ss,2*pi);
phaseVec = [-2*pi, 2*pi];
for ii = 1:N
    iter = 2;
    pp = zeros(7,1);
    pp(1) = (ss(ii,1)-ss(ii,2))^2 + (ss(ii,1)-ss(ii,3))^2 + (ss(ii,2)-ss(ii,3))^2;
    for jj = 1:m
        for kk = 1:2
            ssTemp = ss(ii,:);
            ssTemp(jj) = ssTemp(jj) + phaseVec(kk);
            pp(iter) = (ssTemp(1)-ssTemp(2))^2 + (ssTemp(1)-ssTemp(3))^2 + (ssTemp(2)-ssTemp(3))^2;
            iter = iter + 1;
        end
    end
    
    [~,ch] = min(pp);
    if ch > 1
        if ch == 2
            ss(ii,1) = ss(ii,1) + phaseVec(1);
        elseif ch == 3
            ss(ii,1) = ss(ii,1) + phaseVec(2);
        elseif ch == 4
            ss(ii,2) = ss(ii,2) + phaseVec(1);
        elseif ch == 5
            ss(ii,2) = ss(ii,2) + phaseVec(2);
        elseif ch == 6
            ss(ii,3) = ss(ii,3) + phaseVec(1);
        elseif ch == 7
            ss(ii,3) = ss(ii,3) + phaseVec(2);
        end        
    end
end

mm = ss;
            
