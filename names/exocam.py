"""Model-specific dictionaries of variable names and coordinates."""
from .base import Model

__all__ = ("exocam", )


exocam = Model(
    # Coordinates
    t="time",
    # z="thlev_zsea_theta",
    lev="lev",
    # s="thlev_C_theta",
    y="latitude",
    x="longitude",
    # lw_band="lw_bband_lev_pseudo",
    # Variables
    u="U",
    v="V",
    w="OMEGA",
    # pres=None,
    sh="Q",
    # t_sfc="STASH_m01s00i024",
    toa_net_lw="FLNT",
    toa_net_sw="FSNT",
    toa_isr="FDSTOA",
    # toa_olr=None,
    # toa_olr_cs="STASH_m01s02i206",
    # toa_osr=None,
    # toa_osr_cs="STASH_m01s01i209",
    # sfc_dn_sw="STASH_m01s01i235",
    # sfc_net_down_lw="STASH_m01s02i201",
    # sfc_net_down_sw="m01s01i201",
    # lw_up="STASH_m01s02i509",
    # ocean_frac="STASH_m01s03i395",
    temp="T",
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
