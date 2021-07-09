# -*- coding: utf-8 -*-
"""Model-specific dictionaries of variable names and coordinates."""
from .base import Model

__all__ = ("exocam",)


exocam = Model(
    # Coordinates
    t="time",
    z="z",
    lev="lev",
    # s="thlev_C_theta",
    y="latitude",
    x="longitude",
    # lw_band="lw_bband_lev_pseudo",
    # Variables
    u="U",
    v="V",
    w="OMEGA",
    pres="air_pressure",
    temp_v="virtual_temperature",
    sh="Q",
    t_sfc="TS",
    toa_net_lw="FLNT",
    toa_net_lw_cs="FLNTC",
    toa_net_sw="FSNT",
    toa_net_sw_cs="FSNTC",
    toa_isr="FDSTOA",
    # toa_olr=None,
    # toa_olr_cs="STASH_m01s02i206",
    # toa_osr=None,
    # toa_osr_cs="STASH_m01s01i209",
    sfc_dn_sw="FSDS",
    sfc_net_down_lw="FLNS",
    # sfc_net_down_sw="m01s01i201",
    # lw_up="STASH_m01s02i509",
    ocean_frac="ICEFRAC",
    temp="T",
    rh="RELHUM",
    cld_ice_mf="CLDICE",
    cld_liq_mf="CLDLIQ",
    # cld_ice_v="STASH_m01s00i268",
    # cld_liq_v="STASH_m01s00i267",
    cld_v="CLOUD",
    caf="CLDTOT",
    dt_sw="QRS",
    dt_lw="QRL",
    lwp="TGCLDLWP",
    iwp="TGCLDIWP",
    wvp="TMQ",
)
