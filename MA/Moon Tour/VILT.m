function [DV,r_plus,r_minus] = VILT(rp_iso,N,M,iso_target,iso_dep,mu,type)

    % Check if the target and departure curves are equal --> in this case no
    % VILT is needed
    if iso_target == iso_dep
        
        if N == 0 && M == 0
            ra_plus = iso_target(1,end);
            rp_plus = iso_target(2,end);
            
        else
            Ra_target = interp(iso_target(1,:),1000);
            Rp_target = interp(iso_target(2,:),1000);
            [~,index1] = min(abs(rp_iso(Ra_target,N,M) - Rp_target));
            ra_plus = Ra_target(index1);
            rp_plus = Rp_target(index1);
        end
        r_plus = [ra_plus;rp_plus];
        r_minus = r_plus;
        DV = 0;
        
    else
        % Find the target orbit
        Ra_target = interp(iso_target(1,:),1000);
        Rp_target = interp(iso_target(2,:),1000);
        [~,index1] = min(abs(rp_iso(Ra_target,N,M) - Rp_target));
    
        ra_plus = Ra_target(index1);
        rp_plus = Rp_target(index1);
        r_plus = [ra_plus;rp_plus];
        
    if strcmp(type,'int')
        rp_minus = rp_plus;
        Ra_dep = interp(iso_dep(1,:),1000);
        Rp_dep = interp(iso_dep(2,:),1000);
        [~,index2] = min(abs(Rp_dep - rp_minus));
        ra_minus = Ra_dep(index2);
    else

        % Find the orbit before the VILT
        ra_minus = ra_plus;
        Ra_dep = interp(iso_dep(1,:),1000);
        Rp_dep = interp(iso_dep(2,:),1000);
        [~,index2] = min(abs(Ra_dep - ra_minus));
        rp_minus = Rp_dep(index2);

    end

    r_minus = [ra_minus;rp_minus];
    
        % Compute the VILT
    
       a_plus = (ra_plus + rp_plus)/2;
       a_minus = (ra_minus + rp_minus)/2;
        
       DV = sqrt(mu)*abs(sqrt(2/ra_minus - 1/a_minus) - sqrt(2/ra_plus - 1/a_plus)) ;
    end
end