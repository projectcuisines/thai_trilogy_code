# -*- coding: utf-8 -*-
"""Model-specific dictionaries of variable names and coordinates."""
from .base import Model

__all__ = ("rocke3d",)


rocke3d = Model(
    # Coordinates
    t="time",
    z="z",
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
    temp_v="virtual_temperature",
    sh="q",
    t_sfc="tgrnd",
    toa_net_sw="srnf_toa",
    toa_isr="incsw_toa",
    toa_crf_lw="lwcrf_toa",
    # toa_olr="STASH_m01s02i205",
    toa_olr_cs="olrcs",
    # toa_osr="STASH_m01s01i208",
    toa_osr_cs="swup_toa_clrsky",
    sfc_dn_sw="swds",
    sfc_up_lw="lwus",
    sfc_dn_lw="lwds",
    # sfc_net_down_lw="STASH_m01s02i201",
    # sfc_net_down_sw="m01s01i201",
    # lw_up="STASH_m01s02i509",
    ocean_frac="FOOPN",
    temp="t",
    rh="rh",
    cld_ice_mf="icmmr",
    cld_liq_mf="wcmmr",
    cld_ice_v="icf",
    cld_liq_v="wcf",
    # cld_v="STASH_m01s00i266",
    caf="totcld",
    dt_sw="swhr",
    dt_lw="lwhr",
    lwp="lwp",
    iwp="iwp",
    # wvp="STASH_m01s30i461",
)
