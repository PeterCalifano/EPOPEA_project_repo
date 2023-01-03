function ra_v = plot_isoline(rp_iso,N,M,radius_orb,perc)

    ra0 = fzero(@(ra) rp_iso(ra,N,M) - radius_orb, radius_orb );
    ra_v = linspace(ra0,perc*ra0);

end