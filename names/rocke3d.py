"""Model-specific dictionaries of variable names and coordinates."""
from .base import Model

__all__ = ("rocke3d", )


rocke3d = Model(
    # Coordinates
    t="time",
    # z="thlev_zsea_theta",
    lev="level",
    # s="thlev_C_theta",
    y="latitude",
    x="longitude",
    # lw_band="lw_bband_lev_pseudo",
    # Variables
    u="u",
    v="v",
    w="w",
    pres="p_3d",
    sh="q",
    # t_sfc="STASH_m01s00i024",
    toa_net_sw="srnf_toa",
    toa_isr="incsw_toa",
    toa_crf_lw="lwcrf_toa",
    # toa_olr="STASH_m01s02i205",
    toa_olr_cs="olrcs",
    # toa_osr="STASH_m01s01i208",
    toa_osr_cs="swup_toa_clrsky",
    # sfc_dn_sw="STASH_m01s01i235",
    # sfc_net_down_lw="STASH_m01s02i201",
    # sfc_net_down_sw="m01s01i201",
    # lw_up="STASH_m01s02i509",
    # ocean_frac="STASH_m01s03i395",
    temp="t",
    # rh="STASH_m01s30i113",
    # cld_ice_mf="STASH_m01s00i012",
    # cld_liq_mf="STASH_m01s00i254",
    # cld_ice_v="STASH_m01s00i268",
    # cld_liq_v="STASH_m01s00i267",
    # cld_v="STASH_m01s00i266",
    # caf="STASH_m01s09i217",
    # dt_sw="STASH_m01s01i232",
    # dt_lw="STASH_m01s02i232",
    # lwp="STASH_m01s30i405",
    # iwp="STASH_m01s30i406",
    # wvp="STASH_m01s30i461",
)
