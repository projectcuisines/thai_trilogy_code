# -*- coding: utf-8 -*-
"""Model-specific dictionaries of variable names and coordinates."""
from .base import Model

__all__ = ("lmdg",)


lmdg = Model(
    # Coordinates
    t="Time",
    z="altitude",
    lev="level_height",
    # s="thlev_C_theta",
    y="latitude",
    x="longitude",
    # lw_band="lw_bband_lev_pseudo",
    # Variables
    u="u",
    v="v",
    w="w",
    pres="p",
    temp_v="virtual_temperature",
    sh="h2o_vap",
    t_sfc="tsurf",
    toa_isr="ISR",
    toa_olr="OLR",
    toa_olr_cs="OLRcs",
    toa_net_sw="ASR",
    toa_net_sw_cs="ASRcs",
    # toa_osr_cs="STASH_m01s01i209",
    sfc_dn_sw="fluxsurfsw",
    sfc_net_down_lw="netfluxsurflw",
    # sfc_net_down_sw="m01s01i201",
    # lw_up="STASH_m01s02i509",
    # 1- h2o_ice_surf/1000 for the open ocean fraction
    ocean_frac="h2o_ice_surf",
    temp="temp",
    rh="RH",
    cld_ice_mf="h2o_ice",
    # cld_liq_mf="STASH_m01s00i254",
    # cld_ice_v="STASH_m01s00i268",
    # cld_liq_v="STASH_m01s00i267",
    cld_v="CLF",
    caf="CLFt",
    dt_sw="zdtsw",
    dt_lw="zdtlw",
    cwp="h2o_ice_col",
    # lwp="",
    # iwp="",
    wvp="h2o_vap_col",
    # Temperature threshold for cloud variables:
    #  if T>273K, h2o_ice is liquid,
    #  if T<258K h2o_ice is ice,
    #  if 258K<T<273K, linear interpolation from 100% ice,
    #  0 % liquid at 258K and 1% ice, 100 % liquid at 273K.
)
