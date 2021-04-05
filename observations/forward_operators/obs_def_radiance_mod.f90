! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! This module supports the observation types from the AIRS instruments.
! http://winds.jpl.nasa.gov/missions/quikscat/index.cfm

! BEGIN DART PREPROCESS TYPE DEFINITIONS
! AMSUA_METOP_A_CH1,            QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH2,            QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH3,            QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH4,            QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH5,            QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH6,            QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH7,            QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH8,            QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH9,            QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH10,           QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH11,           QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH12,           QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH13,           QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH14,           QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_METOP_A_CH15,           QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH1,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH2,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH3,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH4,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH5,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH6,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH7,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH8,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH9,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH10,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH11,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH12,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH13,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH14,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N15_CH15,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH1,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH2,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH3,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH4,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH5,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH6,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH7,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH8,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH9,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH10,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH11,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH12,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH13,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH14,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N18_CH15,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH1,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH2,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH3,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH4,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH5,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH6,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH7,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH8,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH9,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH10,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH11,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH12,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH13,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH14,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_N19_CH15,                QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH1,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH2,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH3,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH4,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH5,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH6,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH7,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH8,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH9,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH10,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH11,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH12,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH13,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH14,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! AMSUA_AQUA_CH15,               QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_METOP_A_CH1,              QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_METOP_A_CH2,              QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_METOP_A_CH3,              QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_METOP_A_CH4,              QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_METOP_A_CH5,              QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_N18_CH1,                  QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_N18_CH2,                  QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_N18_CH3,                  QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_N18_CH4,                  QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_N18_CH5,                  QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_N19_CH1,                  QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_N19_CH2,                  QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_N19_CH3,                  QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_N19_CH4,                  QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! MHS_N19_CH5,                  QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! ABI_G16_CH8,                  QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! ABI_G16_CH9,                  QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! ABI_G16_CH10,                 QTY_BRIGHTNESS_TEMPERATURE,        COMMON_CODE
! END DART PREPROCESS TYPE DEFINITIONS
