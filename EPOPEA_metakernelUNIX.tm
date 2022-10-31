KPL/TM

\begintext
This is the default spice kernel pool for EPOPEA project. 
NB: this kernel was generated on Windows PC, file paths and line-ending 
    shall be changed on MacOS and linux.



\begindata

    PATH_VALUES = ('../spice_kernels')

    PATH_SYMBOLS = ('KERNELS')

    KERNELS_TO_LOAD = (
                   '$KERNELS/naif0012.tls',
                   '$KERNELS/gm_de431.tpc', 
                   '$KERNELS/pck00010_mod.tpc',
                   '$KERNELS/de440s.bsp',
                   '$KERNELS/earth_latest_high_prec.bpc', 
                   '$KERNELS/dss.bsp',
                   '$KERNELS/plu058.bsp', 
                   '$KERNELS/dss.tf',
                   '$KERNELS/sat441.bsp',
                   '$KERNELS/LATs.tf',
                   '$KERNELS/LATs.bsp',
)