function ra_v = plot_isoline(rp_iso,N,M,radius_orb,perc)
    
    if N >= M
        ra0 = fzero(@(ra) rp_iso(ra,N,M) - radius_orb, radius_orb);
    else
        ra0 = radius_orb;
    end

    ra_v = linspace(ra0,perc*ra0);
    

end