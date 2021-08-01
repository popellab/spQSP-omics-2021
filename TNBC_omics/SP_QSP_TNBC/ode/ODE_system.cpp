#include "ODE_system.h"
    
#define SPVAR(x) NV_DATA_S(y)[x]
#define NSPVAR(x) ptrOde->_nonspecies_var[x]
#define PARAM(x) _class_parameter[x]
#define PFILE(x) param.getVal(x)

namespace CancerVCT {

#define QSP_W ODE_system::_QSP_weight
bool ODE_system::use_steady_state = false;
bool ODE_system::use_resection = false;
double ODE_system::_QSP_weight = 1.0;


ODE_system::ODE_system()
:CVODEBase()
{
    setupVariables();
    setupEvents();
    setupCVODE();
    update_y_other();
}

ODE_system::ODE_system(const ODE_system& c)
{
    setupCVODE();
}

ODE_system::~ODE_system()
{
}

void ODE_system::initSolver(realtype t){

    restore_y();
    int flag;

    flag = CVodeInit(_cvode_mem, f, t, _y);
    check_flag(&flag, "CVodeInit", 1);

    /* Call CVodeRootInit to specify the root function g */
    flag = CVodeRootInit(_cvode_mem, _nroot, g);
    check_flag(&flag, "CVodeRootInit", 1);
    
    	/*Do not do this. Event only trigger when turn from false to true.
	  If this is reset before trigger evaluation at the beginning of simulation,
	  t=0 events might be missed.*/
    //updateTriggerComponentConditionsOnValue(t);
    //resetEventTriggers();

    return;
}    

state_type ODE_system::_class_parameter = state_type(614, 0);

void ODE_system::setup_class_parameters(Param& param){
    //V_C, mw9284a139_4905_49c2_b085_4f1b555b8ec0, index: 0
    //Unit: metre^(3)
    _class_parameter[0] = PFILE(5) * 0.0010000000000000002;
    //V_P, mw4980ee4d_44c5_41a1_98d0_12d822200107, index: 1
    //Unit: metre^(3)
    _class_parameter[1] = PFILE(6) * 0.0010000000000000002;
    //V_LN, mw437c0fc8_ca27_488b_89a0_ad8e4013c716, index: 2
    //Unit: metre^(3)
    _class_parameter[2] = PFILE(8) * 1.0000000000000013e-09;
    //V_e, mwb1b98857_f6af_41c2_8939_1562424b46a1, index: 3
    //Unit: metre^(3)
    _class_parameter[3] = PFILE(9) * 0.0010000000000000002;
    //A_e, mw974d1d7e_daee_49c2_ba07_5505b61dd47f, index: 4
    //Unit: metre^(2)
    _class_parameter[4] = PFILE(10) * 1e-12;
    //A_s, mw9f2b7328_ff68_419d_af22_30968c278da5, index: 5
    //Unit: metre^(2)
    _class_parameter[5] = PFILE(11) * 1e-12;
    //syn_T_C1, mwbcafcddb_1046_4c20_9043_62ecda4ff607, index: 6
    //Unit: metre^(2)
    _class_parameter[6] = PFILE(12) * 1e-12;
    //syn_T_APC, mwbb03a3bd_68ec_488a_8c10_2af82233b692, index: 7
    //Unit: metre^(2)
    _class_parameter[7] = PFILE(13) * 1e-12;
    //k_cell_clear, mwd71960f5_1662_4b22_9d1f_b6017b8122e6, index: 8
    //Unit: second^(-1)
    _class_parameter[8] = PFILE(310) * 1.15740740740741e-05;
    //cell, mwa1011d3f_410f_4dcc_b8d8_1da2cf3a9185, index: 9
    //Unit: mole^(1)
    _class_parameter[9] = PFILE(311) * 1.66053872801495e-24;
    //day, mw07d56341_4716_436a_b1b8_3dd3260407d0, index: 10
    //Unit: second^(1)
    _class_parameter[10] = PFILE(312) * 86400.0;
    //vol_cell, mw9ab5c38d_eef1_4cb4_9704_da6deb49a711, index: 11
    //Unit: metre^(3)mole^(-1)
    _class_parameter[11] = PFILE(313) * 602214.1989999996;
    //vol_Tcell, mw899f70ba_f953_482d_a1a4_cee6a4d1fad7, index: 12
    //Unit: metre^(3)mole^(-1)
    _class_parameter[12] = PFILE(314) * 602214.1989999996;
    //V_Tmin, mw6cd808bc_3eaf_4d53_aa3b_d669b391c1e4, index: 13
    //Unit: metre^(3)
    _class_parameter[13] = PFILE(315) * 1.0000000000000006e-06;
    //k_C1_growth, mw5182b76e_7e58_4e35_8018_484bbabfd2e2, index: 14
    //Unit: second^(-1)
    _class_parameter[14] = PFILE(329) * 1.15740740740741e-05;
    //C_max, mwdfff49d6_7102_4494_8e62_de0af49f62fa, index: 15
    //Unit: mole^(1)
    _class_parameter[15] = PFILE(330) * 1.66053872801495e-24;
    //k_C1_death, mw6d653759_b926_427b_8601_12ab9bab968f, index: 16
    //Unit: second^(-1)
    _class_parameter[16] = PFILE(331) * 1.15740740740741e-05;
    //k_C1_therapy, mw7a1fd9bb_6e0c_45c0_99cb_194675a848fe, index: 17
    //Unit: second^(-1)
    _class_parameter[17] = PFILE(332) * 1.15740740740741e-05;
    //initial_tumour_diameter, mw5bbec115_04d2_4b95_bd44_0653fed82101, index: 18
    //Unit: metre^(1)
    _class_parameter[18] = PFILE(333) * 0.01;
    //n_T0_clones, mw3bd06d41_4f31_4e13_bc32_238d45425303, index: 19
    //Unit: dimensionless^(1)
    _class_parameter[19] = PFILE(334) * 1.0;
    //Q_T0_in, mwb1326d57_7fa3_4123_81f5_15b3e9b94b3e, index: 20
    //Unit: mole^(1)second^(-1)
    _class_parameter[20] = PFILE(335) * 1.92191982409137e-29;
    //Q_T0_out, mwc62ade04_deb5_4c0b_9965_58977be414b8, index: 21
    //Unit: second^(-1)
    _class_parameter[21] = PFILE(336) * 1.15740740740741e-05;
    //k_T0_act, mw6e8d05a8_bf67_455e_a407_9ac7eebce61e, index: 22
    //Unit: second^(-1)
    _class_parameter[22] = PFILE(337) * 1.15740740740741e-05;
    //k_T0_pro, mw4f20c495_9e4f_433a_b03a_38623a32f56d, index: 23
    //Unit: second^(-1)
    _class_parameter[23] = PFILE(338) * 1.15740740740741e-05;
    //k_T0_death, mw8f33ae75_64a8_4367_8785_b91537aca567, index: 24
    //Unit: second^(-1)
    _class_parameter[24] = PFILE(339) * 1.15740740740741e-05;
    //q_T0_P_in, mw5b10455b_ca07_4a9a_9499_a731690d1613, index: 25
    //Unit: second^(-1)
    _class_parameter[25] = PFILE(340) * 0.0166666666666667;
    //q_T0_P_out, mw25fa377f_0942_44a2_8de2_21edec7015f3, index: 26
    //Unit: second^(-1)
    _class_parameter[26] = PFILE(341) * 1.15740740740741e-05;
    //q_T0_T_in, mwcf9b4fcb_2fb0_4965_aa49_53418a8c96d7, index: 27
    //Unit: metre^(-3)second^(-1)
    _class_parameter[27] = PFILE(342) * 16666.666666666693;
    //q_T0_LN_out, mwab075e91_aa6d_4f17_8576_4740e172a64e, index: 28
    //Unit: second^(-1)
    _class_parameter[28] = PFILE(343) * 1.15740740740741e-05;
    //k_IL2_deg, mw68172044_042c_4da3_b116_85263f229fd8, index: 29
    //Unit: second^(-1)
    _class_parameter[29] = PFILE(344) * 0.0166666666666667;
    //k_IL2_cons, mwb040e6ec_c05e_47c5_a3a6_ed7da6bc3d4b, index: 30
    //Unit: second^(-1)
    _class_parameter[30] = PFILE(345) * 167281721944.444;
    //k_IL2_sec, mw56b1cb91_4a41_4dbc_bcc5_1a67e8b39d4c, index: 31
    //Unit: second^(-1)
    _class_parameter[31] = PFILE(346) * 167281721944.444;
    //IL2_50, mw0ce33574_f53a_4802_b354_79b79046b1e6, index: 32
    //Unit: metre^(-3)mole^(1)
    _class_parameter[32] = PFILE(347) * 1.0000000000000008e-06;
    //IL2_50_Treg, mwd5f7ac2e_505e_4a94_89ce_75627b75fe5a, index: 33
    //Unit: metre^(-3)mole^(1)
    _class_parameter[33] = PFILE(348) * 1.0000000000000008e-06;
    //N0, mw785d231f_bf70_426d_a9a4_2bbb15c3e47b, index: 34
    //Unit: dimensionless^(1)
    _class_parameter[34] = PFILE(349) * 1.0;
    //N_costim, mwab1a352d_5370_4f5a_9b0b_fe7e9f731cc3, index: 35
    //Unit: dimensionless^(1)
    _class_parameter[35] = PFILE(350) * 1.0;
    //N_IL2, mw9528c0b4_d5d6_4072_aeb9_7126777e0adc, index: 36
    //Unit: dimensionless^(1)
    _class_parameter[36] = PFILE(351) * 1.0;
    //k_Treg, mw03e48b27_4075_4b4d_9cc2_d39428205f8c, index: 37
    //Unit: second^(-1)
    _class_parameter[37] = PFILE(352) * 1.15740740740741e-05;
    //n_T1_clones, mwea449865_2da3_4925_83c4_f6d69b8c5a03, index: 38
    //Unit: dimensionless^(1)
    _class_parameter[38] = PFILE(355) * 1.0;
    //Q_T1_in, mwa4a49dc4_599e_4f99_955e_3a07df7a7386, index: 39
    //Unit: mole^(1)second^(-1)
    _class_parameter[39] = PFILE(356) * 1.92191982409137e-29;
    //Q_T1_out, mw3076097a_236b_4222_be77_887b092d5030, index: 40
    //Unit: second^(-1)
    _class_parameter[40] = PFILE(357) * 1.15740740740741e-05;
    //k_T1_act, mw8952ed25_4946_4e38_924b_389da55a2936, index: 41
    //Unit: second^(-1)
    _class_parameter[41] = PFILE(358) * 1.15740740740741e-05;
    //k_T1_pro, mwd7e76bd5_da2c_46a1_a1aa_e649abfb8afb, index: 42
    //Unit: second^(-1)
    _class_parameter[42] = PFILE(359) * 1.15740740740741e-05;
    //k_T1_death, mw17c347b7_cb7c_4947_9a13_4be9ee85a566, index: 43
    //Unit: second^(-1)
    _class_parameter[43] = PFILE(360) * 1.15740740740741e-05;
    //q_T1_P_in, mw874fe11f_ca1e_4f76_932b_5d16e76ec39b, index: 44
    //Unit: second^(-1)
    _class_parameter[44] = PFILE(361) * 0.0166666666666667;
    //q_T1_P_out, mw70b6d4fc_b2e0_4f61_9f9f_63361c79871b, index: 45
    //Unit: second^(-1)
    _class_parameter[45] = PFILE(362) * 1.15740740740741e-05;
    //q_T1_T_in, mw3bed352d_53a7_4a47_98e9_df6c373ea691, index: 46
    //Unit: metre^(-3)second^(-1)
    _class_parameter[46] = PFILE(363) * 16666.666666666693;
    //q_T1_LN_out, mw5ca86745_ed29_4054_a3d4_2e6fcc058412, index: 47
    //Unit: second^(-1)
    _class_parameter[47] = PFILE(364) * 1.15740740740741e-05;
    //k_T1, mw9fa62699_b511_4589_bd54_805757663f11, index: 48
    //Unit: second^(-1)
    _class_parameter[48] = PFILE(365) * 1.15740740740741e-05;
    //k_C_T1, mwa78e360c_986a_4454_b2d8_25565450cb27, index: 49
    //Unit: second^(-1)
    _class_parameter[49] = PFILE(366) * 1.15740740740741e-05;
    //n_T2_clones, mwdec2211e_1121_4606_a483_1f4658d7d0b2, index: 50
    //Unit: dimensionless^(1)
    _class_parameter[50] = PFILE(368) * 1.0;
    //Q_T2_in, mwe8c595c6_0044_48b7_9c60_9c88fa67eb3b, index: 51
    //Unit: mole^(1)second^(-1)
    _class_parameter[51] = PFILE(369) * 1.92191982409137e-29;
    //Q_T2_out, mwd0348803_ce73_4f92_9866_6dfd466b28f4, index: 52
    //Unit: second^(-1)
    _class_parameter[52] = PFILE(370) * 1.15740740740741e-05;
    //k_T2_act, mw11d04e80_346c_4092_9df9_cb291791ec4c, index: 53
    //Unit: second^(-1)
    _class_parameter[53] = PFILE(371) * 1.15740740740741e-05;
    //k_T2_pro, mw591f89ff_d9a7_4a3d_957e_884d2d1fc53a, index: 54
    //Unit: second^(-1)
    _class_parameter[54] = PFILE(372) * 1.15740740740741e-05;
    //k_T2_death, mw724718ed_f49d_475e_a959_517bcb540793, index: 55
    //Unit: second^(-1)
    _class_parameter[55] = PFILE(373) * 1.15740740740741e-05;
    //q_T2_P_in, mw9f6b98c9_b2aa_4f7b_bb78_03e244586a44, index: 56
    //Unit: second^(-1)
    _class_parameter[56] = PFILE(374) * 0.0166666666666667;
    //q_T2_P_out, mwa0070000_67cd_44fc_9b96_7b35c6b7d5c8, index: 57
    //Unit: second^(-1)
    _class_parameter[57] = PFILE(375) * 1.15740740740741e-05;
    //q_T2_T_in, mw497eb8fb_c5de_4c56_b1ce_7a0bdb362dd3, index: 58
    //Unit: metre^(-3)second^(-1)
    _class_parameter[58] = PFILE(376) * 16666.666666666693;
    //q_T2_LN_out, mw25388fc9_1844_457d_abc5_1a12d7204f31, index: 59
    //Unit: second^(-1)
    _class_parameter[59] = PFILE(377) * 1.15740740740741e-05;
    //k_T2, mw37e03dda_4c88_4fef_bb6f_068968892f8a, index: 60
    //Unit: second^(-1)
    _class_parameter[60] = PFILE(378) * 1.15740740740741e-05;
    //k_C_T2, mwb0797b8e_e1ac_4392_9cd6_103e9404e3a2, index: 61
    //Unit: second^(-1)
    _class_parameter[61] = PFILE(379) * 1.15740740740741e-05;
    //n_T3_clones, mw328cf362_44ab_4edf_9a53_3b60948d0b18, index: 62
    //Unit: dimensionless^(1)
    _class_parameter[62] = PFILE(381) * 1.0;
    //Q_T3_in, mw24cb00d1_c65a_49c7_bd10_26e2467b8afd, index: 63
    //Unit: mole^(1)second^(-1)
    _class_parameter[63] = PFILE(382) * 1.92191982409137e-29;
    //Q_T3_out, mwab04704c_81b5_4e5d_965a_0e8b183342cc, index: 64
    //Unit: second^(-1)
    _class_parameter[64] = PFILE(383) * 1.15740740740741e-05;
    //k_T3_act, mw829d59f5_5fb7_40fe_8e82_37daf644b8d6, index: 65
    //Unit: second^(-1)
    _class_parameter[65] = PFILE(384) * 1.15740740740741e-05;
    //k_T3_pro, mwd2f45950_d091_416a_9765_f243c27dfce9, index: 66
    //Unit: second^(-1)
    _class_parameter[66] = PFILE(385) * 1.15740740740741e-05;
    //k_T3_death, mw3d0dab81_52cc_4440_b4d3_e8a235bf4750, index: 67
    //Unit: second^(-1)
    _class_parameter[67] = PFILE(386) * 1.15740740740741e-05;
    //q_T3_P_in, mwdcf65756_9d53_4516_b4cc_59e2036f819a, index: 68
    //Unit: second^(-1)
    _class_parameter[68] = PFILE(387) * 0.0166666666666667;
    //q_T3_P_out, mw4716591d_41a4_4aa2_a3d0_8caf1dbda0a7, index: 69
    //Unit: second^(-1)
    _class_parameter[69] = PFILE(388) * 1.15740740740741e-05;
    //q_T3_T_in, mwe27e69d3_eeae_430c_ab4f_9e47758ed3c3, index: 70
    //Unit: metre^(-3)second^(-1)
    _class_parameter[70] = PFILE(389) * 16666.666666666693;
    //q_T3_LN_out, mw901eb103_fa89_42ed_8c73_0f4c2a1ba7b9, index: 71
    //Unit: second^(-1)
    _class_parameter[71] = PFILE(390) * 1.15740740740741e-05;
    //k_T3, mw99ea9977_a05b_4fa5_97d1_d6429e00b25b, index: 72
    //Unit: second^(-1)
    _class_parameter[72] = PFILE(391) * 1.15740740740741e-05;
    //k_C_T3, mwef58b594_1a03_440c_b6b2_5ab2361aafda, index: 73
    //Unit: second^(-1)
    _class_parameter[73] = PFILE(392) * 1.15740740740741e-05;
    //n_T4_clones, mw5b48d7f4_4418_4802_ab85_b35a1adc5b33, index: 74
    //Unit: dimensionless^(1)
    _class_parameter[74] = PFILE(394) * 1.0;
    //Q_T4_in, mw1781a2f4_8a90_4db5_8059_52b3655828d4, index: 75
    //Unit: mole^(1)second^(-1)
    _class_parameter[75] = PFILE(395) * 1.92191982409137e-29;
    //Q_T4_out, mwc1f8319b_bd35_41d2_8e90_1e07be9436d4, index: 76
    //Unit: second^(-1)
    _class_parameter[76] = PFILE(396) * 1.15740740740741e-05;
    //k_T4_act, mwc7525c32_e52a_4080_b484_035c038370ba, index: 77
    //Unit: second^(-1)
    _class_parameter[77] = PFILE(397) * 1.15740740740741e-05;
    //k_T4_pro, mw93435ca1_ff97_440f_9dd6_dfe65ca248aa, index: 78
    //Unit: second^(-1)
    _class_parameter[78] = PFILE(398) * 1.15740740740741e-05;
    //k_T4_death, mwae49b6bc_a832_457e_8878_4b18409dba3f, index: 79
    //Unit: second^(-1)
    _class_parameter[79] = PFILE(399) * 1.15740740740741e-05;
    //q_T4_P_in, mw8519aa5c_9c54_4da7_9991_568e4c240ada, index: 80
    //Unit: second^(-1)
    _class_parameter[80] = PFILE(400) * 0.0166666666666667;
    //q_T4_P_out, mwdb444db1_29b9_416e_a728_caadd7b5517e, index: 81
    //Unit: second^(-1)
    _class_parameter[81] = PFILE(401) * 1.15740740740741e-05;
    //q_T4_T_in, mw93edfaa4_a81c_4d6e_980e_bd5880d867e5, index: 82
    //Unit: metre^(-3)second^(-1)
    _class_parameter[82] = PFILE(402) * 16666.666666666693;
    //q_T4_LN_out, mw9e899af6_6625_48a3_bd80_c60fae7006f2, index: 83
    //Unit: second^(-1)
    _class_parameter[83] = PFILE(403) * 1.15740740740741e-05;
    //k_T4, mwc4a3c106_4d87_441d_b721_144329ad5b90, index: 84
    //Unit: second^(-1)
    _class_parameter[84] = PFILE(404) * 1.15740740740741e-05;
    //k_C_T4, mw139617ac_064e_46be_9666_b097a3cf1e7e, index: 85
    //Unit: second^(-1)
    _class_parameter[85] = PFILE(405) * 1.15740740740741e-05;
    //n_T5_clones, mw68a57781_5054_4c81_989b_e468a58ddb85, index: 86
    //Unit: dimensionless^(1)
    _class_parameter[86] = PFILE(407) * 1.0;
    //Q_T5_in, mwcf2afb1e_1893_483d_a552_b2e7c2ec2251, index: 87
    //Unit: mole^(1)second^(-1)
    _class_parameter[87] = PFILE(408) * 1.92191982409137e-29;
    //Q_T5_out, mwe1c43569_3cd7_4b05_bba3_40b290b30470, index: 88
    //Unit: second^(-1)
    _class_parameter[88] = PFILE(409) * 1.15740740740741e-05;
    //k_T5_act, mwd0fab0b4_a437_43b4_9e17_dbb32d2c8e5f, index: 89
    //Unit: second^(-1)
    _class_parameter[89] = PFILE(410) * 1.15740740740741e-05;
    //k_T5_pro, mw2c5c2709_c77a_4f05_9804_512e0a96fb16, index: 90
    //Unit: second^(-1)
    _class_parameter[90] = PFILE(411) * 1.15740740740741e-05;
    //k_T5_death, mw6ee12fc9_cee9_4398_8cdb_f67241081808, index: 91
    //Unit: second^(-1)
    _class_parameter[91] = PFILE(412) * 1.15740740740741e-05;
    //q_T5_P_in, mwb526f03b_4934_4c05_bc00_2cbf1e258857, index: 92
    //Unit: second^(-1)
    _class_parameter[92] = PFILE(413) * 0.0166666666666667;
    //q_T5_P_out, mwa4bb061a_330b_407f_a3c0_61355884fc39, index: 93
    //Unit: second^(-1)
    _class_parameter[93] = PFILE(414) * 1.15740740740741e-05;
    //q_T5_T_in, mw2572c8e9_2512_4860_8f18_a32eac5d4a88, index: 94
    //Unit: metre^(-3)second^(-1)
    _class_parameter[94] = PFILE(415) * 16666.666666666693;
    //q_T5_LN_out, mw173f4ab8_ad0f_4f25_9206_c1f300408dcd, index: 95
    //Unit: second^(-1)
    _class_parameter[95] = PFILE(416) * 1.15740740740741e-05;
    //k_T5, mw28f266b2_4948_4458_8538_06e5674d9955, index: 96
    //Unit: second^(-1)
    _class_parameter[96] = PFILE(417) * 1.15740740740741e-05;
    //k_C_T5, mw08699a08_b6cd_47fb_93eb_91abc3be0d6f, index: 97
    //Unit: second^(-1)
    _class_parameter[97] = PFILE(418) * 1.15740740740741e-05;
    //n_T6_clones, mwd27c16f1_5622_40f4_a855_e613305066fd, index: 98
    //Unit: dimensionless^(1)
    _class_parameter[98] = PFILE(420) * 1.0;
    //Q_T6_in, mw003fb718_5e11_4b08_b018_fa150208df55, index: 99
    //Unit: mole^(1)second^(-1)
    _class_parameter[99] = PFILE(421) * 1.92191982409137e-29;
    //Q_T6_out, mwb7f0cb1e_8b3f_41b8_9689_1c78c6984d05, index: 100
    //Unit: second^(-1)
    _class_parameter[100] = PFILE(422) * 1.15740740740741e-05;
    //k_T6_act, mwb239df6b_aad1_4527_ae38_b2015aa8273b, index: 101
    //Unit: second^(-1)
    _class_parameter[101] = PFILE(423) * 1.15740740740741e-05;
    //k_T6_pro, mwc006c2b9_7030_4dbf_8eab_055cd863e20b, index: 102
    //Unit: second^(-1)
    _class_parameter[102] = PFILE(424) * 1.15740740740741e-05;
    //k_T6_death, mw5b5688d0_df5b_4eb4_a1fc_7d14045c41bd, index: 103
    //Unit: second^(-1)
    _class_parameter[103] = PFILE(425) * 1.15740740740741e-05;
    //q_T6_P_in, mw19c5bdbb_1482_4a2d_b327_9bbe1d48eb80, index: 104
    //Unit: second^(-1)
    _class_parameter[104] = PFILE(426) * 0.0166666666666667;
    //q_T6_P_out, mw217264ce_412f_4a45_b1d1_c313192dc367, index: 105
    //Unit: second^(-1)
    _class_parameter[105] = PFILE(427) * 1.15740740740741e-05;
    //q_T6_T_in, mwa1e4ec2a_aabf_4002_92ec_9b9761131b57, index: 106
    //Unit: metre^(-3)second^(-1)
    _class_parameter[106] = PFILE(428) * 16666.666666666693;
    //q_T6_LN_out, mw6f09385e_3531_4023_8c80_d699d4f89ffa, index: 107
    //Unit: second^(-1)
    _class_parameter[107] = PFILE(429) * 1.15740740740741e-05;
    //k_T6, mwa078cd66_0a89_4814_a8d4_66273facdc3e, index: 108
    //Unit: second^(-1)
    _class_parameter[108] = PFILE(430) * 1.15740740740741e-05;
    //k_C_T6, mw4176d329_3d74_425e_88b8_e50d9cd85f15, index: 109
    //Unit: second^(-1)
    _class_parameter[109] = PFILE(431) * 1.15740740740741e-05;
    //n_T7_clones, mw0ea51dcb_ec20_446b_afeb_6c58616e13ae, index: 110
    //Unit: dimensionless^(1)
    _class_parameter[110] = PFILE(433) * 1.0;
    //Q_T7_in, mw44c40498_5396_43f2_97d7_b9058a32a2fa, index: 111
    //Unit: mole^(1)second^(-1)
    _class_parameter[111] = PFILE(434) * 1.92191982409137e-29;
    //Q_T7_out, mw9609411c_ba38_4a10_bb42_c13a1932a3fd, index: 112
    //Unit: second^(-1)
    _class_parameter[112] = PFILE(435) * 1.15740740740741e-05;
    //k_T7_act, mw27ea2e2e_396b_4dba_88fa_b83ee5f6a89b, index: 113
    //Unit: second^(-1)
    _class_parameter[113] = PFILE(436) * 1.15740740740741e-05;
    //k_T7_pro, mw7e58c7ac_0a00_428d_be04_68ef4ccd8255, index: 114
    //Unit: second^(-1)
    _class_parameter[114] = PFILE(437) * 1.15740740740741e-05;
    //k_T7_death, mw6a9aa03a_9fb8_449f_b7bb_fd309e9efdc5, index: 115
    //Unit: second^(-1)
    _class_parameter[115] = PFILE(438) * 1.15740740740741e-05;
    //q_T7_P_in, mwa388e882_60a8_4970_9cee_f0fb3743b981, index: 116
    //Unit: second^(-1)
    _class_parameter[116] = PFILE(439) * 0.0166666666666667;
    //q_T7_P_out, mwd02fcf1d_1f53_4651_8db5_33615614ea2e, index: 117
    //Unit: second^(-1)
    _class_parameter[117] = PFILE(440) * 1.15740740740741e-05;
    //q_T7_T_in, mw6727fc88_b9c6_47fa_afb4_6fb748669e9a, index: 118
    //Unit: metre^(-3)second^(-1)
    _class_parameter[118] = PFILE(441) * 16666.666666666693;
    //q_T7_LN_out, mw613e724c_24c6_4aba_aa74_ec8a417e09ec, index: 119
    //Unit: second^(-1)
    _class_parameter[119] = PFILE(442) * 1.15740740740741e-05;
    //k_T7, mwa6998d68_0818_4c8c_8d72_d88029b8ed60, index: 120
    //Unit: second^(-1)
    _class_parameter[120] = PFILE(443) * 1.15740740740741e-05;
    //k_C_T7, mwa1050278_acb0_4c88_9d01_2ca00989912c, index: 121
    //Unit: second^(-1)
    _class_parameter[121] = PFILE(444) * 1.15740740740741e-05;
    //n_T8_clones, mwa440c9d5_7f52_4055_be9e_c39746b05d99, index: 122
    //Unit: dimensionless^(1)
    _class_parameter[122] = PFILE(446) * 1.0;
    //Q_T8_in, mwbdb36bbc_bcda_42e7_8d8c_4fec9c29ac86, index: 123
    //Unit: mole^(1)second^(-1)
    _class_parameter[123] = PFILE(447) * 1.92191982409137e-29;
    //Q_T8_out, mw576f531a_7994_41fd_8ea9_427766fa9708, index: 124
    //Unit: second^(-1)
    _class_parameter[124] = PFILE(448) * 1.15740740740741e-05;
    //k_T8_act, mw469bfbee_7e55_438e_a521_74324d710789, index: 125
    //Unit: second^(-1)
    _class_parameter[125] = PFILE(449) * 1.15740740740741e-05;
    //k_T8_pro, mwc0a5a842_c0b1_4f7d_b563_6017e0c920c3, index: 126
    //Unit: second^(-1)
    _class_parameter[126] = PFILE(450) * 1.15740740740741e-05;
    //k_T8_death, mw7ea4a09a_3add_4e02_98a7_21fac964dc28, index: 127
    //Unit: second^(-1)
    _class_parameter[127] = PFILE(451) * 1.15740740740741e-05;
    //q_T8_P_in, mw7decf784_b0f9_49a0_911e_a7b8c322a6a1, index: 128
    //Unit: second^(-1)
    _class_parameter[128] = PFILE(452) * 0.0166666666666667;
    //q_T8_P_out, mwd08fe5c9_9843_4cae_9e7b_24247890495e, index: 129
    //Unit: second^(-1)
    _class_parameter[129] = PFILE(453) * 1.15740740740741e-05;
    //q_T8_T_in, mw11a5da32_c798_4623_a0c4_ff086b0a4c68, index: 130
    //Unit: metre^(-3)second^(-1)
    _class_parameter[130] = PFILE(454) * 16666.666666666693;
    //q_T8_LN_out, mw76c823aa_f2de_494d_92be_8cedeab90117, index: 131
    //Unit: second^(-1)
    _class_parameter[131] = PFILE(455) * 1.15740740740741e-05;
    //k_T8, mw667185eb_db65_40f3_a562_26d5d20853c8, index: 132
    //Unit: second^(-1)
    _class_parameter[132] = PFILE(456) * 1.15740740740741e-05;
    //k_C_T8, mwc511855f_1cbd_468a_ba11_1b1db30f3932, index: 133
    //Unit: second^(-1)
    _class_parameter[133] = PFILE(457) * 1.15740740740741e-05;
    //n_T9_clones, mw538bd258_0bfb_4de9_abb8_388f8011e28c, index: 134
    //Unit: dimensionless^(1)
    _class_parameter[134] = PFILE(459) * 1.0;
    //Q_T9_in, mwb5430389_5faa_4cfd_aab6_64834abaf069, index: 135
    //Unit: mole^(1)second^(-1)
    _class_parameter[135] = PFILE(460) * 1.92191982409137e-29;
    //Q_T9_out, mw90965d2c_9460_480f_9d90_8d0d57806a81, index: 136
    //Unit: second^(-1)
    _class_parameter[136] = PFILE(461) * 1.15740740740741e-05;
    //k_T9_act, mw9d4b0306_fc33_46d4_a979_8374f5798dd0, index: 137
    //Unit: second^(-1)
    _class_parameter[137] = PFILE(462) * 1.15740740740741e-05;
    //k_T9_pro, mw45fbf76f_6748_42c5_a83a_105c97134aa4, index: 138
    //Unit: second^(-1)
    _class_parameter[138] = PFILE(463) * 1.15740740740741e-05;
    //k_T9_death, mw92650f60_b9fb_47a6_9a59_d34a6a93c53d, index: 139
    //Unit: second^(-1)
    _class_parameter[139] = PFILE(464) * 1.15740740740741e-05;
    //q_T9_P_in, mwb86ea2de_b20e_4b18_9681_f8eb27d0195a, index: 140
    //Unit: second^(-1)
    _class_parameter[140] = PFILE(465) * 0.0166666666666667;
    //q_T9_P_out, mw66ba46c6_73e5_4004_8d7c_9c159070cb35, index: 141
    //Unit: second^(-1)
    _class_parameter[141] = PFILE(466) * 1.15740740740741e-05;
    //q_T9_T_in, mw8298bcbb_d69f_4e74_af94_cf2e29dc6ab4, index: 142
    //Unit: metre^(-3)second^(-1)
    _class_parameter[142] = PFILE(467) * 16666.666666666693;
    //q_T9_LN_out, mw16065964_2b7f_431e_9ecd_09c61ef7240d, index: 143
    //Unit: second^(-1)
    _class_parameter[143] = PFILE(468) * 1.15740740740741e-05;
    //k_T9, mwd5f899e6_d58e_4ed4_97f1_8715bf01de07, index: 144
    //Unit: second^(-1)
    _class_parameter[144] = PFILE(469) * 1.15740740740741e-05;
    //k_C_T9, mw436d7335_57e3_42f8_a8bb_d7dbcfa42115, index: 145
    //Unit: second^(-1)
    _class_parameter[145] = PFILE(470) * 1.15740740740741e-05;
    //n_T10_clones, mwf6753623_bf88_40ea_aef4_a246ccbbb24b, index: 146
    //Unit: dimensionless^(1)
    _class_parameter[146] = PFILE(472) * 1.0;
    //Q_T10_in, mw697f6309_0ee7_4693_8edd_bc674cab85a3, index: 147
    //Unit: mole^(1)second^(-1)
    _class_parameter[147] = PFILE(473) * 1.92191982409137e-29;
    //Q_T10_out, mw72ab6360_efaf_46c3_b064_22ba49e84864, index: 148
    //Unit: second^(-1)
    _class_parameter[148] = PFILE(474) * 1.15740740740741e-05;
    //k_T10_act, mw1c7506c5_e0fd_449e_ae83_2f2f757bd879, index: 149
    //Unit: second^(-1)
    _class_parameter[149] = PFILE(475) * 1.15740740740741e-05;
    //k_T10_pro, mw35b2af30_edf1_4fd1_8059_f1d68819f5ca, index: 150
    //Unit: second^(-1)
    _class_parameter[150] = PFILE(476) * 1.15740740740741e-05;
    //k_T10_death, mwd52dd43e_5a20_40f5_88a1_3b9be87b082d, index: 151
    //Unit: second^(-1)
    _class_parameter[151] = PFILE(477) * 1.15740740740741e-05;
    //q_T10_P_in, mw919e546e_72ec_4706_bba8_d70b4188f266, index: 152
    //Unit: second^(-1)
    _class_parameter[152] = PFILE(478) * 0.0166666666666667;
    //q_T10_P_out, mw7e403267_0ff5_479c_8145_3da7e02a3142, index: 153
    //Unit: second^(-1)
    _class_parameter[153] = PFILE(479) * 1.15740740740741e-05;
    //q_T10_T_in, mwb5cb32e1_3b85_42ef_8d47_f1177f13e84b, index: 154
    //Unit: metre^(-3)second^(-1)
    _class_parameter[154] = PFILE(480) * 16666.666666666693;
    //q_T10_LN_out, mw5679e3b1_0bcd_4443_adf6_f58c5767c28f, index: 155
    //Unit: second^(-1)
    _class_parameter[155] = PFILE(481) * 1.15740740740741e-05;
    //k_T10, mw2384aae1_3939_45f0_93c7_4f35e51a5c77, index: 156
    //Unit: second^(-1)
    _class_parameter[156] = PFILE(482) * 1.15740740740741e-05;
    //k_C_T10, mw7dcdb3a3_6151_48d0_9283_cdf1496cd282, index: 157
    //Unit: second^(-1)
    _class_parameter[157] = PFILE(483) * 1.15740740740741e-05;
    //n_T11_clones, mw9afc5c57_aa70_46fd_8919_3e0754bc16ad, index: 158
    //Unit: dimensionless^(1)
    _class_parameter[158] = PFILE(485) * 1.0;
    //Q_T11_in, mwfbf60e7c_84b5_40b3_9392_fbc560a0e0d4, index: 159
    //Unit: mole^(1)second^(-1)
    _class_parameter[159] = PFILE(486) * 1.92191982409137e-29;
    //Q_T11_out, mw1f23fcbc_a169_48e6_b646_c9c156e35edc, index: 160
    //Unit: second^(-1)
    _class_parameter[160] = PFILE(487) * 1.15740740740741e-05;
    //k_T11_act, mw2df65934_eefd_4581_b658_a448cb809c0e, index: 161
    //Unit: second^(-1)
    _class_parameter[161] = PFILE(488) * 1.15740740740741e-05;
    //k_T11_pro, mw634c4f24_63ae_4eca_9bf2_0434cd8cd858, index: 162
    //Unit: second^(-1)
    _class_parameter[162] = PFILE(489) * 1.15740740740741e-05;
    //k_T11_death, mwdb99f9d3_0cb5_4a3d_9814_63ada55dd498, index: 163
    //Unit: second^(-1)
    _class_parameter[163] = PFILE(490) * 1.15740740740741e-05;
    //q_T11_P_in, mw7d2222d8_74c2_48e8_befe_0c9d3500038d, index: 164
    //Unit: second^(-1)
    _class_parameter[164] = PFILE(491) * 0.0166666666666667;
    //q_T11_P_out, mwb01f525e_b9bc_49cc_bf3e_acc8eaa15b96, index: 165
    //Unit: second^(-1)
    _class_parameter[165] = PFILE(492) * 1.15740740740741e-05;
    //q_T11_T_in, mw680a3013_5ad3_4b6d_951d_989c91407b35, index: 166
    //Unit: metre^(-3)second^(-1)
    _class_parameter[166] = PFILE(493) * 16666.666666666693;
    //q_T11_LN_out, mw75fda0c3_9a76_41ec_8d54_7f5052c3a4d4, index: 167
    //Unit: second^(-1)
    _class_parameter[167] = PFILE(494) * 1.15740740740741e-05;
    //k_T11, mwc8e82180_017e_4605_a04f_5b2df14a56d9, index: 168
    //Unit: second^(-1)
    _class_parameter[168] = PFILE(495) * 1.15740740740741e-05;
    //k_C_T11, mw3cf0090a_81ad_4c2a_9098_24c407b94566, index: 169
    //Unit: second^(-1)
    _class_parameter[169] = PFILE(496) * 1.15740740740741e-05;
    //n_T12_clones, mw363f8525_0b3f_4872_bd3b_7c20e46551f1, index: 170
    //Unit: dimensionless^(1)
    _class_parameter[170] = PFILE(498) * 1.0;
    //Q_T12_in, mw2be26f51_35eb_4917_8e97_f154f58c0fc7, index: 171
    //Unit: mole^(1)second^(-1)
    _class_parameter[171] = PFILE(499) * 1.92191982409137e-29;
    //Q_T12_out, mw7f2982a7_c3b4_4073_a62c_2d25736ad485, index: 172
    //Unit: second^(-1)
    _class_parameter[172] = PFILE(500) * 1.15740740740741e-05;
    //k_T12_act, mwbfeb6c0c_861e_4c1a_83b9_2897c0a3eda3, index: 173
    //Unit: second^(-1)
    _class_parameter[173] = PFILE(501) * 1.15740740740741e-05;
    //k_T12_pro, mw7d1fca68_5920_4185_9033_a251f2402a14, index: 174
    //Unit: second^(-1)
    _class_parameter[174] = PFILE(502) * 1.15740740740741e-05;
    //k_T12_death, mw12e2fd89_5b38_459f_8be7_f7f559e637e9, index: 175
    //Unit: second^(-1)
    _class_parameter[175] = PFILE(503) * 1.15740740740741e-05;
    //q_T12_P_in, mwa2b0f349_cd93_4298_b403_25e116ed925e, index: 176
    //Unit: second^(-1)
    _class_parameter[176] = PFILE(504) * 0.0166666666666667;
    //q_T12_P_out, mw211ba81f_fc4c_47d5_bac8_fd8e9f13bf18, index: 177
    //Unit: second^(-1)
    _class_parameter[177] = PFILE(505) * 1.15740740740741e-05;
    //q_T12_T_in, mw32800965_54fd_4a26_87b9_302ef38b8b48, index: 178
    //Unit: metre^(-3)second^(-1)
    _class_parameter[178] = PFILE(506) * 16666.666666666693;
    //q_T12_LN_out, mwa7dd638f_a94b_445f_8ba0_093f31856a02, index: 179
    //Unit: second^(-1)
    _class_parameter[179] = PFILE(507) * 1.15740740740741e-05;
    //k_T12, mwc2e325d5_2550_439f_bf05_48ba705ddd76, index: 180
    //Unit: second^(-1)
    _class_parameter[180] = PFILE(508) * 1.15740740740741e-05;
    //k_C_T12, mw9abb78bf_b991_42b3_8e1e_dd2a436151d3, index: 181
    //Unit: second^(-1)
    _class_parameter[181] = PFILE(509) * 1.15740740740741e-05;
    //n_T13_clones, mw2a333031_e6c9_4967_b091_d00078a7d762, index: 182
    //Unit: dimensionless^(1)
    _class_parameter[182] = PFILE(511) * 1.0;
    //Q_T13_in, mwb15e3b5a_ee55_4e9f_a24f_cf792cb15c70, index: 183
    //Unit: mole^(1)second^(-1)
    _class_parameter[183] = PFILE(512) * 1.92191982409137e-29;
    //Q_T13_out, mwe97bd858_53f1_431a_b643_1b1dc33c22ba, index: 184
    //Unit: second^(-1)
    _class_parameter[184] = PFILE(513) * 1.15740740740741e-05;
    //k_T13_act, mw358b7d0e_a517_4455_a589_114d80bfd5bc, index: 185
    //Unit: second^(-1)
    _class_parameter[185] = PFILE(514) * 1.15740740740741e-05;
    //k_T13_pro, mwc6a4df6d_3471_4d31_8b02_43b47d2e1f3c, index: 186
    //Unit: second^(-1)
    _class_parameter[186] = PFILE(515) * 1.15740740740741e-05;
    //k_T13_death, mw4cf68282_7f87_4f31_9fe6_c19993aac16f, index: 187
    //Unit: second^(-1)
    _class_parameter[187] = PFILE(516) * 1.15740740740741e-05;
    //q_T13_P_in, mw8d5c8ca5_dcce_4f42_9f17_1eed5fdfbec5, index: 188
    //Unit: second^(-1)
    _class_parameter[188] = PFILE(517) * 0.0166666666666667;
    //q_T13_P_out, mw671ea70b_e040_4d2b_bc4f_0358a0204d9e, index: 189
    //Unit: second^(-1)
    _class_parameter[189] = PFILE(518) * 1.15740740740741e-05;
    //q_T13_T_in, mwdbae01a0_8219_434b_adee_8258c4b8d713, index: 190
    //Unit: metre^(-3)second^(-1)
    _class_parameter[190] = PFILE(519) * 16666.666666666693;
    //q_T13_LN_out, mwd561fb03_ded3_4cc3_bd4a_35c3223e4694, index: 191
    //Unit: second^(-1)
    _class_parameter[191] = PFILE(520) * 1.15740740740741e-05;
    //k_T13, mw6feb3f34_23d9_48ba_b2c5_3606348bd301, index: 192
    //Unit: second^(-1)
    _class_parameter[192] = PFILE(521) * 1.15740740740741e-05;
    //k_C_T13, mw57679c7c_4a71_44db_b8d2_813cf2e45738, index: 193
    //Unit: second^(-1)
    _class_parameter[193] = PFILE(522) * 1.15740740740741e-05;
    //n_T14_clones, mwa294c4b3_4a80_4a4d_a375_c44646087e46, index: 194
    //Unit: dimensionless^(1)
    _class_parameter[194] = PFILE(524) * 1.0;
    //Q_T14_in, mwc9465b41_7e0d_4906_a1c7_ebf56c3be6dc, index: 195
    //Unit: mole^(1)second^(-1)
    _class_parameter[195] = PFILE(525) * 1.92191982409137e-29;
    //Q_T14_out, mwc8d46719_121b_494e_b091_51f689d7a797, index: 196
    //Unit: second^(-1)
    _class_parameter[196] = PFILE(526) * 1.15740740740741e-05;
    //k_T14_act, mwef0cd5d7_db97_4f35_970c_15861614fe7a, index: 197
    //Unit: second^(-1)
    _class_parameter[197] = PFILE(527) * 1.15740740740741e-05;
    //k_T14_pro, mw39380552_c721_45ac_9f44_c993da344e9d, index: 198
    //Unit: second^(-1)
    _class_parameter[198] = PFILE(528) * 1.15740740740741e-05;
    //k_T14_death, mw9d1bacde_b3d3_40ce_8f62_caa0af93c422, index: 199
    //Unit: second^(-1)
    _class_parameter[199] = PFILE(529) * 1.15740740740741e-05;
    //q_T14_P_in, mwd4c42467_caf1_4fad_96d0_7bd2692782c5, index: 200
    //Unit: second^(-1)
    _class_parameter[200] = PFILE(530) * 0.0166666666666667;
    //q_T14_P_out, mwa12b61d9_ad3e_4f5b_9d15_d7aa94b8a6ea, index: 201
    //Unit: second^(-1)
    _class_parameter[201] = PFILE(531) * 1.15740740740741e-05;
    //q_T14_T_in, mweec242f2_fa66_4695_8551_b905a1bb2c27, index: 202
    //Unit: metre^(-3)second^(-1)
    _class_parameter[202] = PFILE(532) * 16666.666666666693;
    //q_T14_LN_out, mwa9b037e4_9da7_414e_8bc1_3ffbebc01aaf, index: 203
    //Unit: second^(-1)
    _class_parameter[203] = PFILE(533) * 1.15740740740741e-05;
    //k_T14, mw8bbfa42e_56b0_4cc3_8ded_57d0a2556f02, index: 204
    //Unit: second^(-1)
    _class_parameter[204] = PFILE(534) * 1.15740740740741e-05;
    //k_C_T14, mw578e9acd_43c1_4ea5_9374_ddccda94790d, index: 205
    //Unit: second^(-1)
    _class_parameter[205] = PFILE(535) * 1.15740740740741e-05;
    //n_T15_clones, mwfa5be138_6f5b_4647_915b_3f6df8c8fd4d, index: 206
    //Unit: dimensionless^(1)
    _class_parameter[206] = PFILE(537) * 1.0;
    //Q_T15_in, mw349afc58_09c1_4dad_a2af_d067551fe65e, index: 207
    //Unit: mole^(1)second^(-1)
    _class_parameter[207] = PFILE(538) * 1.92191982409137e-29;
    //Q_T15_out, mwc561b1f7_55c9_454c_b0d9_b6327b5a4fab, index: 208
    //Unit: second^(-1)
    _class_parameter[208] = PFILE(539) * 1.15740740740741e-05;
    //k_T15_act, mw77eb17f6_ec25_47cc_8c50_b897e81423aa, index: 209
    //Unit: second^(-1)
    _class_parameter[209] = PFILE(540) * 1.15740740740741e-05;
    //k_T15_pro, mw053149fb_bce8_4801_bf0c_3bae857c960f, index: 210
    //Unit: second^(-1)
    _class_parameter[210] = PFILE(541) * 1.15740740740741e-05;
    //k_T15_death, mw7d9ef403_9d01_4c08_94f5_31b738d6783d, index: 211
    //Unit: second^(-1)
    _class_parameter[211] = PFILE(542) * 1.15740740740741e-05;
    //q_T15_P_in, mw4a649314_49ca_48f6_9926_3e149269f184, index: 212
    //Unit: second^(-1)
    _class_parameter[212] = PFILE(543) * 0.0166666666666667;
    //q_T15_P_out, mw78c280bd_690d_4292_91a1_8e9ca3e09acf, index: 213
    //Unit: second^(-1)
    _class_parameter[213] = PFILE(544) * 1.15740740740741e-05;
    //q_T15_T_in, mwff6408bc_ca7b_49d2_8d58_f9169899bbfd, index: 214
    //Unit: metre^(-3)second^(-1)
    _class_parameter[214] = PFILE(545) * 16666.666666666693;
    //q_T15_LN_out, mw44a423f6_3dc2_4a1c_837b_318986a02052, index: 215
    //Unit: second^(-1)
    _class_parameter[215] = PFILE(546) * 1.15740740740741e-05;
    //k_T15, mwe2bff87c_fb5c_49c7_b4f8_b0416cfc2d15, index: 216
    //Unit: second^(-1)
    _class_parameter[216] = PFILE(547) * 1.15740740740741e-05;
    //k_C_T15, mw7bf7e6d8_6cf2_4824_8b73_938d73b5329f, index: 217
    //Unit: second^(-1)
    _class_parameter[217] = PFILE(548) * 1.15740740740741e-05;
    //n_T16_clones, mw6c7840f4_d2c7_4ac0_ab86_4604b82897d3, index: 218
    //Unit: dimensionless^(1)
    _class_parameter[218] = PFILE(550) * 1.0;
    //Q_T16_in, mwacee25d1_ec3e_46d8_92d9_7476b42f8431, index: 219
    //Unit: mole^(1)second^(-1)
    _class_parameter[219] = PFILE(551) * 1.92191982409137e-29;
    //Q_T16_out, mwff175bc6_3c4f_46e3_bc9c_ae0b91f3748e, index: 220
    //Unit: second^(-1)
    _class_parameter[220] = PFILE(552) * 1.15740740740741e-05;
    //k_T16_act, mw82b8b440_fde9_4f8f_bfdb_5ab59efdeaac, index: 221
    //Unit: second^(-1)
    _class_parameter[221] = PFILE(553) * 1.15740740740741e-05;
    //k_T16_pro, mw3633fa87_7b5b_46bd_8944_c2dd7ba0249e, index: 222
    //Unit: second^(-1)
    _class_parameter[222] = PFILE(554) * 1.15740740740741e-05;
    //k_T16_death, mw1be2673e_06eb_4cb9_adbc_69f0869165fb, index: 223
    //Unit: second^(-1)
    _class_parameter[223] = PFILE(555) * 1.15740740740741e-05;
    //q_T16_P_in, mw360aba1d_bc30_4390_9a73_4a83acdeedfe, index: 224
    //Unit: second^(-1)
    _class_parameter[224] = PFILE(556) * 0.0166666666666667;
    //q_T16_P_out, mw1ffadcc7_1660_4fc9_ae6a_64801866a37f, index: 225
    //Unit: second^(-1)
    _class_parameter[225] = PFILE(557) * 1.15740740740741e-05;
    //q_T16_T_in, mwdd6098e8_034b_44e6_b112_3cb5bc279ec4, index: 226
    //Unit: metre^(-3)second^(-1)
    _class_parameter[226] = PFILE(558) * 16666.666666666693;
    //q_T16_LN_out, mwf114ddc8_f537_418c_a626_87fe9b33e6b2, index: 227
    //Unit: second^(-1)
    _class_parameter[227] = PFILE(559) * 1.15740740740741e-05;
    //k_T16, mw90663403_b58d_415c_9651_9deb749473a9, index: 228
    //Unit: second^(-1)
    _class_parameter[228] = PFILE(560) * 1.15740740740741e-05;
    //k_C_T16, mw211a9c47_a936_4636_99b5_400e6f69e0e3, index: 229
    //Unit: second^(-1)
    _class_parameter[229] = PFILE(561) * 1.15740740740741e-05;
    //n_T17_clones, mw63e32f5a_54b0_45ec_9cd2_de315a96408f, index: 230
    //Unit: dimensionless^(1)
    _class_parameter[230] = PFILE(563) * 1.0;
    //Q_T17_in, mw74a36d9f_cb1a_4d3c_815a_46fb826bcd9e, index: 231
    //Unit: mole^(1)second^(-1)
    _class_parameter[231] = PFILE(564) * 1.92191982409137e-29;
    //Q_T17_out, mw7bc2b879_82ef_4e4b_88cc_ece07a0af63e, index: 232
    //Unit: second^(-1)
    _class_parameter[232] = PFILE(565) * 1.15740740740741e-05;
    //k_T17_act, mw4c52d501_2d4d_4d10_93dc_e2dd721115e2, index: 233
    //Unit: second^(-1)
    _class_parameter[233] = PFILE(566) * 1.15740740740741e-05;
    //k_T17_pro, mw1fc64367_e0d5_4d31_9a13_75fe0fa9b1b8, index: 234
    //Unit: second^(-1)
    _class_parameter[234] = PFILE(567) * 1.15740740740741e-05;
    //k_T17_death, mwb93fb2a1_07d3_441b_8c19_322da6af9ed9, index: 235
    //Unit: second^(-1)
    _class_parameter[235] = PFILE(568) * 1.15740740740741e-05;
    //q_T17_P_in, mwbd39a596_492b_4706_886e_701f9e05428d, index: 236
    //Unit: second^(-1)
    _class_parameter[236] = PFILE(569) * 0.0166666666666667;
    //q_T17_P_out, mw04311bc4_4bd2_4f4d_89f2_e23293c03d21, index: 237
    //Unit: second^(-1)
    _class_parameter[237] = PFILE(570) * 1.15740740740741e-05;
    //q_T17_T_in, mw170e27a3_99fd_4673_b9df_845827a22240, index: 238
    //Unit: metre^(-3)second^(-1)
    _class_parameter[238] = PFILE(571) * 16666.666666666693;
    //q_T17_LN_out, mw90512b6a_e42c_4aa9_8516_77a1b795b7a0, index: 239
    //Unit: second^(-1)
    _class_parameter[239] = PFILE(572) * 1.15740740740741e-05;
    //k_T17, mw64d79f1e_93e3_4dae_a5a0_c3ecfe6b34a2, index: 240
    //Unit: second^(-1)
    _class_parameter[240] = PFILE(573) * 1.15740740740741e-05;
    //k_C_T17, mwdb120b5f_b278_4f09_9ac6_39f78945cebd, index: 241
    //Unit: second^(-1)
    _class_parameter[241] = PFILE(574) * 1.15740740740741e-05;
    //k_APC_mat, mwb0ef13c3_e61d_46cb_bec5_4ad32b6bd9b8, index: 242
    //Unit: second^(-1)
    _class_parameter[242] = PFILE(576) * 1.15740740740741e-05;
    //k_APC_mig, mw4a059608_ae48_406a_84f0_846b2b097892, index: 243
    //Unit: second^(-1)
    _class_parameter[243] = PFILE(577) * 1.15740740740741e-05;
    //k_APC_death, mwe5f5eeb7_0771_4c3e_90ba_c5c7c4b84bc0, index: 244
    //Unit: second^(-1)
    _class_parameter[244] = PFILE(578) * 1.15740740740741e-05;
    //k_mAPC_death, mwfda375e0_7055_4fcf_96a4_5c7171e17b03, index: 245
    //Unit: second^(-1)
    _class_parameter[245] = PFILE(579) * 1.15740740740741e-05;
    //APC0_T, mwafedb5f6_990b_4a44_bb23_953492c5e2e8, index: 246
    //Unit: metre^(-3)mole^(1)
    _class_parameter[246] = PFILE(580) * 1.6605387280149534e-18;
    //APC0_LN, mwc9036db5_7b09_4edc_a9b3_f320e1264675, index: 247
    //Unit: metre^(-3)mole^(1)
    _class_parameter[247] = PFILE(581) * 1.6605387280149534e-18;
    //k_c, mw7d19df24_aefa_42b8_a7f3_60d35f358ee6, index: 248
    //Unit: second^(-1)
    _class_parameter[248] = PFILE(582) * 1.15740740740741e-05;
    //c0, mwdc1faa36_b099_44c8_90b7_c0d5b6623265, index: 249
    //Unit: metre^(-3)mole^(1)
    _class_parameter[249] = PFILE(583) * 999.9999999999994;
    //c50, mw6da1a69d_a37e_4a39_abfc_c5d64be20169, index: 250
    //Unit: metre^(-3)mole^(1)
    _class_parameter[250] = PFILE(584) * 999.9999999999994;
    //DAMPs, mw875f757c_ec44_44cb_b6d9_eae94b4d81ce, index: 251
    //Unit: dimensionless^(1)
    _class_parameter[251] = PFILE(585) * 6.02214199e+23;
    //n_sites_APC, mw5ee2be99_6ae0_45ae_a10a_0b2cfc2d7cba, index: 252
    //Unit: dimensionless^(1)
    _class_parameter[252] = PFILE(586) * 1.0;
    //kin, mwd308ad80_f4c3_4c8c_8c27_4be35a5860e9, index: 253
    //Unit: second^(-1)
    _class_parameter[253] = PFILE(587) * 1.15740740740741e-05;
    //kout, mw90dc0fd6_7889_4ada_836c_af0f01b5e7e8, index: 254
    //Unit: second^(-1)
    _class_parameter[254] = PFILE(588) * 1.15740740740741e-05;
    //k_P0_up, mw017c86ce_440a_4c08_8deb_6463e4aec578, index: 255
    //Unit: mole^(-1)second^(-1)
    _class_parameter[255] = PFILE(589) * 6.97007174768519e+18;
    //k_xP0_deg, mw831db726_1d9d_43b9_88d0_7e56bb71d6ca, index: 256
    //Unit: second^(-1)
    _class_parameter[256] = PFILE(590) * 1.15740740740741e-05;
    //k_P0_deg, mwb5008dbc_2048_439f_b7a7_8db964582fba, index: 257
    //Unit: second^(-1)
    _class_parameter[257] = PFILE(591) * 1.15740740740741e-05;
    //k_p0_deg, mwcab8fc74_2387_4dac_905d_b43c78b67cb2, index: 258
    //Unit: second^(-1)
    _class_parameter[258] = PFILE(592) * 1.15740740740741e-05;
    //k_P0_on, mw9af118d4_793e_4d9b_be46_07e752e3f44f, index: 259
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[259] = PFILE(593) * 1.1574074074074112e-08;
    //k_P0_d1, mweb470b9d_a32d_485f_bbc8_245462866134, index: 260
    //Unit: metre^(-3)mole^(1)
    _class_parameter[260] = PFILE(594) * 999.9999999999994;
    //p0_50, mw1a1aca5c_e807_4805_b975_e14f09ce6bf9, index: 261
    //Unit: metre^(-2)mole^(1)
    _class_parameter[261] = PFILE(595) * 1.66053872801495e-12;
    //P0_C1, mwf81473a3_937c_4cca_9a08_843b453145ee, index: 262
    //Unit: metre^(-3)
    _class_parameter[262] = PFILE(596) * 6.022141989999979e+26;
    //A_syn, mw417eda6d_4943_42cf_b7f9_e0861c537885, index: 263
    //Unit: metre^(2)
    _class_parameter[263] = PFILE(597) * 1e-12;
    //A_Tcell, mw911bfe5d_26e2_468d_927e_831c1ae164eb, index: 264
    //Unit: metre^(2)
    _class_parameter[264] = PFILE(598) * 1e-12;
    //A_cell, mw417534ec_c150_4279_b9a3_397d8ce508cc, index: 265
    //Unit: metre^(2)
    _class_parameter[265] = PFILE(599) * 1e-12;
    //A_APC, mwa1356868_423a_4af0_bccd_4d9e625f5c4b, index: 266
    //Unit: metre^(2)
    _class_parameter[266] = PFILE(600) * 1e-12;
    //k_M1p0_TCR_on, mwea30906d_bd38_410c_a9f6_7f495fa50176, index: 267
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[267] = PFILE(601) * 602214199000.0;
    //k_M1p0_TCR_off, mwfd6fad6d_79c4_4a33_b789_5672fdf11de9, index: 268
    //Unit: second^(-1)
    _class_parameter[268] = PFILE(602) * 1.0;
    //TCR_p0_tot, mw79144076_94b5_410c_bcba_82d671115654, index: 269
    //Unit: metre^(-2)mole^(1)
    _class_parameter[269] = PFILE(603) * 1.66053872801495e-12;
    //k_P1_up, mw2ed5881b_e4f1_408c_b8f6_84cb25b3e55b, index: 270
    //Unit: mole^(-1)second^(-1)
    _class_parameter[270] = PFILE(605) * 6.97007174768519e+18;
    //k_xP1_deg, mw503e1832_5e4e_4629_9f4e_da8358b89296, index: 271
    //Unit: second^(-1)
    _class_parameter[271] = PFILE(606) * 1.15740740740741e-05;
    //k_P1_deg, mwd11fc5fd_e271_461d_9940_33c96bcc8347, index: 272
    //Unit: second^(-1)
    _class_parameter[272] = PFILE(607) * 1.15740740740741e-05;
    //k_p1_deg, mwfe5395cb_6a3b_4fc7_a829_6fc0022620ef, index: 273
    //Unit: second^(-1)
    _class_parameter[273] = PFILE(608) * 1.15740740740741e-05;
    //k_P1_on, mw1c38d665_81ed_4418_8676_9149d0349679, index: 274
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[274] = PFILE(609) * 1.1574074074074112e-08;
    //k_P1_d1, mw34873426_3072_4e70_b877_4ce58e05e564, index: 275
    //Unit: metre^(-3)mole^(1)
    _class_parameter[275] = PFILE(610) * 999.9999999999994;
    //p1_50, mw1e956a72_b095_42b9_95ea_fe3e6dde685e, index: 276
    //Unit: metre^(-2)mole^(1)
    _class_parameter[276] = PFILE(611) * 1.66053872801495e-12;
    //P1_C1, mw977cd85b_da2e_45b0_bb08_a7bb514700cc, index: 277
    //Unit: metre^(-3)
    _class_parameter[277] = PFILE(612) * 6.022141989999979e+26;
    //k_M1p1_TCR_on, mw0ff8aa52_f5fd_44f5_ba4f_cb389d2a31a5, index: 278
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[278] = PFILE(613) * 602214199000.0;
    //k_M1p1_TCR_off, mw9cdfb535_e69a_4c4c_b632_367ef93eb438, index: 279
    //Unit: second^(-1)
    _class_parameter[279] = PFILE(614) * 1.0;
    //k_M1p1_TCR_p, mwee24520f_bdcf_42f0_9514_c14f392e022a, index: 280
    //Unit: second^(-1)
    _class_parameter[280] = PFILE(615) * 1.0;
    //phi_M1p1_TCR, mwd4157901_9add_4c4d_a1b4_747916795ee3, index: 281
    //Unit: second^(-1)
    _class_parameter[281] = PFILE(616) * 1.0;
    //N_M1p1_TCR, mw99b14fb5_2a18_44a2_9980_8c96af2322c4, index: 282
    //Unit: dimensionless^(1)
    _class_parameter[282] = PFILE(617) * 1.0;
    //TCR_p1_tot, mw0445a51d_b6a6_416b_9a4b_6d7e750428f0, index: 283
    //Unit: metre^(-2)mole^(1)
    _class_parameter[283] = PFILE(618) * 1.66053872801495e-12;
    //k_P2_up, mwba74d1c6_7b9f_458b_ba3c_18a7d3313998, index: 284
    //Unit: mole^(-1)second^(-1)
    _class_parameter[284] = PFILE(620) * 6.97007174768519e+18;
    //k_xP2_deg, mw7126c072_cca6_48a7_8aff_2d41953618aa, index: 285
    //Unit: second^(-1)
    _class_parameter[285] = PFILE(621) * 1.15740740740741e-05;
    //k_P2_deg, mw1d19b685_c1e2_4ad2_a384_189b1aab676d, index: 286
    //Unit: second^(-1)
    _class_parameter[286] = PFILE(622) * 1.15740740740741e-05;
    //k_p2_deg, mwea9dc036_6163_437c_97ba_0b37642e8b47, index: 287
    //Unit: second^(-1)
    _class_parameter[287] = PFILE(623) * 1.15740740740741e-05;
    //k_P2_on, mw3f9d2793_1d2f_4a99_a5a9_2ecddbf615b1, index: 288
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[288] = PFILE(624) * 1.1574074074074112e-08;
    //k_P2_d1, mw86562714_ab1b_4231_b109_9250714a849c, index: 289
    //Unit: metre^(-3)mole^(1)
    _class_parameter[289] = PFILE(625) * 999.9999999999994;
    //p2_50, mw4b621d3a_f4c4_446c_9864_97da5057df12, index: 290
    //Unit: metre^(-2)mole^(1)
    _class_parameter[290] = PFILE(626) * 1.66053872801495e-12;
    //P2_C1, mwf7e66649_b90d_4acd_8150_aee9f2e2c29a, index: 291
    //Unit: metre^(-3)
    _class_parameter[291] = PFILE(627) * 6.022141989999979e+26;
    //k_M1p2_TCR_on, mwe1b377a3_5c29_4364_b8c6_aaa352af7b6f, index: 292
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[292] = PFILE(628) * 602214199000.0;
    //k_M1p2_TCR_off, mw78893b0d_751c_40fb_9ccd_3d39003ee68c, index: 293
    //Unit: second^(-1)
    _class_parameter[293] = PFILE(629) * 1.0;
    //k_M1p2_TCR_p, mwcc1988ee_36a6_42ff_b600_ef4993757051, index: 294
    //Unit: second^(-1)
    _class_parameter[294] = PFILE(630) * 1.0;
    //phi_M1p2_TCR, mw2b94cba8_f028_44ec_a748_065116b30d6a, index: 295
    //Unit: second^(-1)
    _class_parameter[295] = PFILE(631) * 1.0;
    //N_M1p2_TCR, mw0bc2a080_1e44_4cee_881b_0f3bb57724d8, index: 296
    //Unit: dimensionless^(1)
    _class_parameter[296] = PFILE(632) * 1.0;
    //TCR_p2_tot, mw8419dca7_d2cc_40f6_a71d_a843304d50ab, index: 297
    //Unit: metre^(-2)mole^(1)
    _class_parameter[297] = PFILE(633) * 1.66053872801495e-12;
    //k_P3_up, mw941bf168_bd41_45f4_85d6_a52850f9c2b0, index: 298
    //Unit: mole^(-1)second^(-1)
    _class_parameter[298] = PFILE(635) * 6.97007174768519e+18;
    //k_xP3_deg, mw1468a3a2_5a9b_405f_a42f_782a5c093f07, index: 299
    //Unit: second^(-1)
    _class_parameter[299] = PFILE(636) * 1.15740740740741e-05;
    //k_P3_deg, mw7f6e8ebd_5f54_4b1c_afe6_07ed83c95b8a, index: 300
    //Unit: second^(-1)
    _class_parameter[300] = PFILE(637) * 1.15740740740741e-05;
    //k_p3_deg, mw636ccce5_bf7f_4034_9e64_599dc3906668, index: 301
    //Unit: second^(-1)
    _class_parameter[301] = PFILE(638) * 1.15740740740741e-05;
    //k_P3_on, mwa5b36634_387c_4b4e_8658_5b1e14cd4929, index: 302
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[302] = PFILE(639) * 1.1574074074074112e-08;
    //k_P3_d1, mw9e5f5df7_a666_4292_9be1_9afe836fcf7a, index: 303
    //Unit: metre^(-3)mole^(1)
    _class_parameter[303] = PFILE(640) * 999.9999999999994;
    //p3_50, mw2458ceeb_89b7_4930_a364_7dc8d844ae1e, index: 304
    //Unit: metre^(-2)mole^(1)
    _class_parameter[304] = PFILE(641) * 1.66053872801495e-12;
    //P3_C1, mw2dd8ab85_ba2c_49b5_a686_0fec349821c8, index: 305
    //Unit: metre^(-3)
    _class_parameter[305] = PFILE(642) * 6.022141989999979e+26;
    //k_M1p3_TCR_on, mw77259946_df93_4a70_881d_c1a537e4f041, index: 306
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[306] = PFILE(643) * 602214199000.0;
    //k_M1p3_TCR_off, mw385d0cc8_81f2_4c55_b96d_211ead454bd4, index: 307
    //Unit: second^(-1)
    _class_parameter[307] = PFILE(644) * 1.0;
    //k_M1p3_TCR_p, mwe3867384_a411_40cf_be14_a8c34417ef65, index: 308
    //Unit: second^(-1)
    _class_parameter[308] = PFILE(645) * 1.0;
    //phi_M1p3_TCR, mwab848304_eb4b_4946_bf4c_537eb40cd2d9, index: 309
    //Unit: second^(-1)
    _class_parameter[309] = PFILE(646) * 1.0;
    //N_M1p3_TCR, mwb6007dbf_9205_482f_8277_f6216b1cc390, index: 310
    //Unit: dimensionless^(1)
    _class_parameter[310] = PFILE(647) * 1.0;
    //TCR_p3_tot, mw54d97cfa_909c_4855_81bd_9db6b157fe1c, index: 311
    //Unit: metre^(-2)mole^(1)
    _class_parameter[311] = PFILE(648) * 1.66053872801495e-12;
    //k_P4_up, mw9b033f03_c809_4d58_94bf_bbf59e416eb4, index: 312
    //Unit: mole^(-1)second^(-1)
    _class_parameter[312] = PFILE(650) * 6.97007174768519e+18;
    //k_xP4_deg, mwbcf6bfdc_d36d_4899_a1e4_551d8af81cd2, index: 313
    //Unit: second^(-1)
    _class_parameter[313] = PFILE(651) * 1.15740740740741e-05;
    //k_P4_deg, mw669dd490_37b1_46a6_88b9_01c0980f85f6, index: 314
    //Unit: second^(-1)
    _class_parameter[314] = PFILE(652) * 1.15740740740741e-05;
    //k_p4_deg, mwfe1d10a6_35e8_44f1_9927_2c717b25a35d, index: 315
    //Unit: second^(-1)
    _class_parameter[315] = PFILE(653) * 1.15740740740741e-05;
    //k_P4_on, mwf274edd5_e72c_4c36_8537_f3959d648c14, index: 316
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[316] = PFILE(654) * 1.1574074074074112e-08;
    //k_P4_d1, mwcab0dc31_db71_45b3_bd5f_06002a1c6d62, index: 317
    //Unit: metre^(-3)mole^(1)
    _class_parameter[317] = PFILE(655) * 999.9999999999994;
    //p4_50, mw4bfdf969_c5f4_4e9b_82b2_6bf5ca4b16bb, index: 318
    //Unit: metre^(-2)mole^(1)
    _class_parameter[318] = PFILE(656) * 1.66053872801495e-12;
    //P4_C1, mwef22207f_c785_4cf5_a733_46aedc35fe31, index: 319
    //Unit: metre^(-3)
    _class_parameter[319] = PFILE(657) * 6.022141989999979e+26;
    //k_M1p4_TCR_on, mwbc1ed9ba_7ccb_49a0_be71_17593b9eca92, index: 320
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[320] = PFILE(658) * 602214199000.0;
    //k_M1p4_TCR_off, mwac650f0e_bf25_4cb4_99f4_3c8be9a72587, index: 321
    //Unit: second^(-1)
    _class_parameter[321] = PFILE(659) * 1.0;
    //k_M1p4_TCR_p, mwe9ddad2e_5e44_4269_9721_6d907613ba2f, index: 322
    //Unit: second^(-1)
    _class_parameter[322] = PFILE(660) * 1.0;
    //phi_M1p4_TCR, mwe4ee7ea2_baac_4c30_9b06_c94e56813feb, index: 323
    //Unit: second^(-1)
    _class_parameter[323] = PFILE(661) * 1.0;
    //N_M1p4_TCR, mwc132c073_3350_49e1_966e_ef1fa066980f, index: 324
    //Unit: dimensionless^(1)
    _class_parameter[324] = PFILE(662) * 1.0;
    //TCR_p4_tot, mw41fd5f0d_0e18_4335_b552_5656d98abfa9, index: 325
    //Unit: metre^(-2)mole^(1)
    _class_parameter[325] = PFILE(663) * 1.66053872801495e-12;
    //k_P5_up, mw3db7a56c_c25a_4119_bf83_a1e073360b68, index: 326
    //Unit: mole^(-1)second^(-1)
    _class_parameter[326] = PFILE(665) * 6.97007174768519e+18;
    //k_xP5_deg, mw33eec3ce_665a_4c9e_86a4_20fce877ce0e, index: 327
    //Unit: second^(-1)
    _class_parameter[327] = PFILE(666) * 1.15740740740741e-05;
    //k_P5_deg, mw5d03e414_12bc_403c_a7c7_1abfedf1bf21, index: 328
    //Unit: second^(-1)
    _class_parameter[328] = PFILE(667) * 1.15740740740741e-05;
    //k_p5_deg, mw45637eed_9cfa_49e6_9388_eef5cec0df90, index: 329
    //Unit: second^(-1)
    _class_parameter[329] = PFILE(668) * 1.15740740740741e-05;
    //k_P5_on, mwce17b5c7_4065_4087_91fb_0412908afd04, index: 330
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[330] = PFILE(669) * 1.1574074074074112e-08;
    //k_P5_d1, mw0711f760_dc47_478f_9609_011772888050, index: 331
    //Unit: metre^(-3)mole^(1)
    _class_parameter[331] = PFILE(670) * 999.9999999999994;
    //p5_50, mwc7d3be93_be01_4f16_8349_723a553f86ef, index: 332
    //Unit: metre^(-2)mole^(1)
    _class_parameter[332] = PFILE(671) * 1.66053872801495e-12;
    //P5_C1, mwda381890_572f_4212_8e35_25acd84e7650, index: 333
    //Unit: metre^(-3)
    _class_parameter[333] = PFILE(672) * 6.022141989999979e+26;
    //k_M1p5_TCR_on, mwd0c783fc_47ce_4295_a6dc_30e6dbea01c9, index: 334
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[334] = PFILE(673) * 602214199000.0;
    //k_M1p5_TCR_off, mw7c4f24c5_9ecb_4abe_8fa0_70935bc884c8, index: 335
    //Unit: second^(-1)
    _class_parameter[335] = PFILE(674) * 1.0;
    //k_M1p5_TCR_p, mwc20bf4df_c85a_41ed_9652_fdddf5a66b28, index: 336
    //Unit: second^(-1)
    _class_parameter[336] = PFILE(675) * 1.0;
    //phi_M1p5_TCR, mwabcf4515_2711_4920_9c4b_c37d43245b67, index: 337
    //Unit: second^(-1)
    _class_parameter[337] = PFILE(676) * 1.0;
    //N_M1p5_TCR, mwf02dddd0_5897_438f_8194_362f90cd1e7d, index: 338
    //Unit: dimensionless^(1)
    _class_parameter[338] = PFILE(677) * 1.0;
    //TCR_p5_tot, mw325a40f1_5f05_42c7_9cbd_d8cde93fccbc, index: 339
    //Unit: metre^(-2)mole^(1)
    _class_parameter[339] = PFILE(678) * 1.66053872801495e-12;
    //k_P6_up, mwc1f411ac_c726_4ddc_8980_8008527cc5c0, index: 340
    //Unit: mole^(-1)second^(-1)
    _class_parameter[340] = PFILE(680) * 6.97007174768519e+18;
    //k_xP6_deg, mw58604e6c_ff21_4766_be30_7c564036a3e7, index: 341
    //Unit: second^(-1)
    _class_parameter[341] = PFILE(681) * 1.15740740740741e-05;
    //k_P6_deg, mw596c7162_8d24_4bc0_a7d9_efa140276f26, index: 342
    //Unit: second^(-1)
    _class_parameter[342] = PFILE(682) * 1.15740740740741e-05;
    //k_p6_deg, mwa8e0e539_aa5f_4c0f_bcbd_956121c6dbae, index: 343
    //Unit: second^(-1)
    _class_parameter[343] = PFILE(683) * 1.15740740740741e-05;
    //k_P6_on, mw981981ee_b67b_43d4_a8b3_8a27ca739880, index: 344
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[344] = PFILE(684) * 1.1574074074074112e-08;
    //k_P6_d1, mw13075872_edb3_4a4e_8a87_821ff146a162, index: 345
    //Unit: metre^(-3)mole^(1)
    _class_parameter[345] = PFILE(685) * 999.9999999999994;
    //p6_50, mwc996b24b_b918_40c8_8c81_bf2a4aa27208, index: 346
    //Unit: metre^(-2)mole^(1)
    _class_parameter[346] = PFILE(686) * 1.66053872801495e-12;
    //P6_C1, mwa6fa47c1_875f_4504_81bc_db51bd4b04d4, index: 347
    //Unit: metre^(-3)
    _class_parameter[347] = PFILE(687) * 6.022141989999979e+26;
    //k_M1p6_TCR_on, mwb494f703_24da_42ab_8941_75e6ff4e90ab, index: 348
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[348] = PFILE(688) * 602214199000.0;
    //k_M1p6_TCR_off, mwbd89a406_aa02_4bc0_adea_c2fa7a0f7dca, index: 349
    //Unit: second^(-1)
    _class_parameter[349] = PFILE(689) * 1.0;
    //k_M1p6_TCR_p, mw712dd34c_cb50_44cb_9550_952458c8860b, index: 350
    //Unit: second^(-1)
    _class_parameter[350] = PFILE(690) * 1.0;
    //phi_M1p6_TCR, mw97477bfc_cfb2_44bd_960d_1436e5963f07, index: 351
    //Unit: second^(-1)
    _class_parameter[351] = PFILE(691) * 1.0;
    //N_M1p6_TCR, mw7c33514d_9e5d_4938_9249_065a78159f5a, index: 352
    //Unit: dimensionless^(1)
    _class_parameter[352] = PFILE(692) * 1.0;
    //TCR_p6_tot, mwe9f9123a_8898_4d40_af61_162e04edfa1f, index: 353
    //Unit: metre^(-2)mole^(1)
    _class_parameter[353] = PFILE(693) * 1.66053872801495e-12;
    //k_P7_up, mwc48e155f_77cc_43ca_9c6a_cbf2c77440f7, index: 354
    //Unit: mole^(-1)second^(-1)
    _class_parameter[354] = PFILE(695) * 6.97007174768519e+18;
    //k_xP7_deg, mwe3d53c16_5709_4f6d_b85d_3b2cd458e160, index: 355
    //Unit: second^(-1)
    _class_parameter[355] = PFILE(696) * 1.15740740740741e-05;
    //k_P7_deg, mw999d1230_1bfe_4b4e_a246_c753846a7e72, index: 356
    //Unit: second^(-1)
    _class_parameter[356] = PFILE(697) * 1.15740740740741e-05;
    //k_p7_deg, mw61a6eca3_d1d7_48cb_9ab2_43050a1767bd, index: 357
    //Unit: second^(-1)
    _class_parameter[357] = PFILE(698) * 1.15740740740741e-05;
    //k_P7_on, mw7f05f9fa_5066_44d3_917b_da543cfa2907, index: 358
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[358] = PFILE(699) * 1.1574074074074112e-08;
    //k_P7_d1, mwd4a515b5_a0e5_4446_a8cd_a80ec0177d9d, index: 359
    //Unit: metre^(-3)mole^(1)
    _class_parameter[359] = PFILE(700) * 999.9999999999994;
    //p7_50, mw6670c5f6_2c22_4610_b9e7_a8aac870a014, index: 360
    //Unit: metre^(-2)mole^(1)
    _class_parameter[360] = PFILE(701) * 1.66053872801495e-12;
    //P7_C1, mwde2d7506_4ec2_4d79_ab8d_fce5e1ba4200, index: 361
    //Unit: metre^(-3)
    _class_parameter[361] = PFILE(702) * 6.022141989999979e+26;
    //k_M1p7_TCR_on, mw6aa1f115_d8de_4ffa_a0a7_c488f1dde033, index: 362
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[362] = PFILE(703) * 602214199000.0;
    //k_M1p7_TCR_off, mw0a68bf84_c6d1_408c_8cdf_b5b4b56b9fc6, index: 363
    //Unit: second^(-1)
    _class_parameter[363] = PFILE(704) * 1.0;
    //k_M1p7_TCR_p, mw2b5f9f54_fe20_4863_a29c_7b470cca0b0e, index: 364
    //Unit: second^(-1)
    _class_parameter[364] = PFILE(705) * 1.0;
    //phi_M1p7_TCR, mw4155caa8_c2fb_4af6_aa87_d87c1c10bec9, index: 365
    //Unit: second^(-1)
    _class_parameter[365] = PFILE(706) * 1.0;
    //N_M1p7_TCR, mwa775a1d8_6698_4c8a_88d7_9bcabba9e316, index: 366
    //Unit: dimensionless^(1)
    _class_parameter[366] = PFILE(707) * 1.0;
    //TCR_p7_tot, mw67561213_67f3_466e_84d1_a40d2db6b607, index: 367
    //Unit: metre^(-2)mole^(1)
    _class_parameter[367] = PFILE(708) * 1.66053872801495e-12;
    //k_P8_up, mwa6bd2520_b379_4e19_957a_189999a6a231, index: 368
    //Unit: mole^(-1)second^(-1)
    _class_parameter[368] = PFILE(710) * 6.97007174768519e+18;
    //k_xP8_deg, mw9d84abd2_6118_4d15_94f8_1c84879e094d, index: 369
    //Unit: second^(-1)
    _class_parameter[369] = PFILE(711) * 1.15740740740741e-05;
    //k_P8_deg, mwd2510216_587a_4896_a24a_afe3f2c7bfe7, index: 370
    //Unit: second^(-1)
    _class_parameter[370] = PFILE(712) * 1.15740740740741e-05;
    //k_p8_deg, mw17cc4ea3_0707_4e95_b3bb_49fb4146d5d0, index: 371
    //Unit: second^(-1)
    _class_parameter[371] = PFILE(713) * 1.15740740740741e-05;
    //k_P8_on, mw2a975c7d_b424_4b85_a08e_aa9100b111c5, index: 372
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[372] = PFILE(714) * 1.1574074074074112e-08;
    //k_P8_d1, mweb61ee85_2afa_46b5_b246_fdb763ffc0fa, index: 373
    //Unit: metre^(-3)mole^(1)
    _class_parameter[373] = PFILE(715) * 999.9999999999994;
    //p8_50, mw68a4486d_f0ed_4762_94da_ab802d5b8fe1, index: 374
    //Unit: metre^(-2)mole^(1)
    _class_parameter[374] = PFILE(716) * 1.66053872801495e-12;
    //P8_C1, mwdd75344a_c492_4482_9ebd_418a3eb52b24, index: 375
    //Unit: metre^(-3)
    _class_parameter[375] = PFILE(717) * 6.022141989999979e+26;
    //k_M1p8_TCR_on, mw807c962f_6ae0_4bff_a909_67b72e9b764e, index: 376
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[376] = PFILE(718) * 602214199000.0;
    //k_M1p8_TCR_off, mw0f7f77d5_16a8_43fd_88cd_54a93ea58222, index: 377
    //Unit: second^(-1)
    _class_parameter[377] = PFILE(719) * 1.0;
    //k_M1p8_TCR_p, mw4237ab67_d030_4d9f_a6c0_46acb0aacdc4, index: 378
    //Unit: second^(-1)
    _class_parameter[378] = PFILE(720) * 1.0;
    //phi_M1p8_TCR, mwcbd41acd_a4ca_4aa5_b55f_b380404100a7, index: 379
    //Unit: second^(-1)
    _class_parameter[379] = PFILE(721) * 1.0;
    //N_M1p8_TCR, mw025e6c89_8f42_4495_8d97_de9e515eeed9, index: 380
    //Unit: dimensionless^(1)
    _class_parameter[380] = PFILE(722) * 1.0;
    //TCR_p8_tot, mwf5428332_3373_49b6_adf5_823dc20616c2, index: 381
    //Unit: metre^(-2)mole^(1)
    _class_parameter[381] = PFILE(723) * 1.66053872801495e-12;
    //k_P9_up, mw0e016d98_90f4_4865_9a7c_f7b1a24fdf81, index: 382
    //Unit: mole^(-1)second^(-1)
    _class_parameter[382] = PFILE(725) * 6.97007174768519e+18;
    //k_xP9_deg, mw7973b47a_ff98_4b4a_be16_70d633b6f2b3, index: 383
    //Unit: second^(-1)
    _class_parameter[383] = PFILE(726) * 1.15740740740741e-05;
    //k_P9_deg, mw37343031_0c72_41a8_9114_9839cf9dd4fc, index: 384
    //Unit: second^(-1)
    _class_parameter[384] = PFILE(727) * 1.15740740740741e-05;
    //k_p9_deg, mwf6cb0f21_c35e_4bc4_a738_7174d55cad04, index: 385
    //Unit: second^(-1)
    _class_parameter[385] = PFILE(728) * 1.15740740740741e-05;
    //k_P9_on, mwe4abff26_9f18_485d_b2b3_6bbcc32873da, index: 386
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[386] = PFILE(729) * 1.1574074074074112e-08;
    //k_P9_d1, mw24a30218_d677_4821_a937_8f5327b3cf41, index: 387
    //Unit: metre^(-3)mole^(1)
    _class_parameter[387] = PFILE(730) * 999.9999999999994;
    //p9_50, mwa3c11297_dbfc_4531_a93f_6161b4df9364, index: 388
    //Unit: metre^(-2)mole^(1)
    _class_parameter[388] = PFILE(731) * 1.66053872801495e-12;
    //P9_C1, mwf327f365_6457_4c88_a8df_f590b9be787c, index: 389
    //Unit: metre^(-3)
    _class_parameter[389] = PFILE(732) * 6.022141989999979e+26;
    //k_M1p9_TCR_on, mw07106ffe_efbe_4fe3_80dd_efade7e75d26, index: 390
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[390] = PFILE(733) * 602214199000.0;
    //k_M1p9_TCR_off, mw768f68fe_1078_4010_925e_3b846a8433d1, index: 391
    //Unit: second^(-1)
    _class_parameter[391] = PFILE(734) * 1.0;
    //k_M1p9_TCR_p, mw20cdf02b_4da6_413d_b01e_8bb661c8a2c5, index: 392
    //Unit: second^(-1)
    _class_parameter[392] = PFILE(735) * 1.0;
    //phi_M1p9_TCR, mw3d793d5d_f464_4ff0_a49a_254a287d20d5, index: 393
    //Unit: second^(-1)
    _class_parameter[393] = PFILE(736) * 1.0;
    //N_M1p9_TCR, mw2765cd48_2c61_4fee_85aa_c50bc4be0dd4, index: 394
    //Unit: dimensionless^(1)
    _class_parameter[394] = PFILE(737) * 1.0;
    //TCR_p9_tot, mw9a7d72e1_8e2d_4054_90b6_b24b55d5f4fb, index: 395
    //Unit: metre^(-2)mole^(1)
    _class_parameter[395] = PFILE(738) * 1.66053872801495e-12;
    //k_P10_up, mw20153460_1cb0_4b9b_b817_191c77fd5c13, index: 396
    //Unit: mole^(-1)second^(-1)
    _class_parameter[396] = PFILE(740) * 6.97007174768519e+18;
    //k_xP10_deg, mw4394fe95_6a19_4e88_8f1c_7f9c2adbf745, index: 397
    //Unit: second^(-1)
    _class_parameter[397] = PFILE(741) * 1.15740740740741e-05;
    //k_P10_deg, mw5df99daf_b3de_4b02_b54b_9e8d6fcfc2b9, index: 398
    //Unit: second^(-1)
    _class_parameter[398] = PFILE(742) * 1.15740740740741e-05;
    //k_p10_deg, mw1e5efbaf_9baf_474d_a7bf_fda23274848f, index: 399
    //Unit: second^(-1)
    _class_parameter[399] = PFILE(743) * 1.15740740740741e-05;
    //k_P10_on, mw0111a3c5_441e_490e_b9d7_515f59092ba1, index: 400
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[400] = PFILE(744) * 1.1574074074074112e-08;
    //k_P10_d1, mwfe19a6a5_829c_42e5_a2d5_ee4e25d73108, index: 401
    //Unit: metre^(-3)mole^(1)
    _class_parameter[401] = PFILE(745) * 999.9999999999994;
    //p10_50, mweab692ff_689d_41d4_b553_1054f6aa1294, index: 402
    //Unit: metre^(-2)mole^(1)
    _class_parameter[402] = PFILE(746) * 1.66053872801495e-12;
    //P10_C1, mw26a1100a_e71e_4d6a_95fb_fe7e281c646c, index: 403
    //Unit: metre^(-3)
    _class_parameter[403] = PFILE(747) * 6.022141989999979e+26;
    //k_M1p10_TCR_on, mw7eb42c73_43d5_4092_bc77_2b15ffa8d253, index: 404
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[404] = PFILE(748) * 602214199000.0;
    //k_M1p10_TCR_off, mw789f6fd0_8e53_4693_a82d_adcd304a9228, index: 405
    //Unit: second^(-1)
    _class_parameter[405] = PFILE(749) * 1.0;
    //k_M1p10_TCR_p, mw34181542_8f36_466a_9bdc_3650165fccf7, index: 406
    //Unit: second^(-1)
    _class_parameter[406] = PFILE(750) * 1.0;
    //phi_M1p10_TCR, mw5a797c6b_ed8d_46df_87ad_d5217d2ed6ce, index: 407
    //Unit: second^(-1)
    _class_parameter[407] = PFILE(751) * 1.0;
    //N_M1p10_TCR, mw7f4b4c02_face_4e6e_8469_ca878fa0fcc1, index: 408
    //Unit: dimensionless^(1)
    _class_parameter[408] = PFILE(752) * 1.0;
    //TCR_p10_tot, mw79b9f003_db4b_4ab4_910d_0fb0b813050c, index: 409
    //Unit: metre^(-2)mole^(1)
    _class_parameter[409] = PFILE(753) * 1.66053872801495e-12;
    //k_P11_up, mw0fe7e0b6_fc63_4221_b750_df23f5c8f286, index: 410
    //Unit: mole^(-1)second^(-1)
    _class_parameter[410] = PFILE(755) * 6.97007174768519e+18;
    //k_xP11_deg, mw83cdd5ea_cc5e_4888_8808_a4472fc6986f, index: 411
    //Unit: second^(-1)
    _class_parameter[411] = PFILE(756) * 1.15740740740741e-05;
    //k_P11_deg, mwddb6533a_88ed_4003_ac58_20d3fee6e2c1, index: 412
    //Unit: second^(-1)
    _class_parameter[412] = PFILE(757) * 1.15740740740741e-05;
    //k_p11_deg, mw503439a7_19e6_46c5_8172_1c39d07ff524, index: 413
    //Unit: second^(-1)
    _class_parameter[413] = PFILE(758) * 1.15740740740741e-05;
    //k_P11_on, mwc15d6cdc_71d6_4eae_abf3_b24367e6c4b4, index: 414
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[414] = PFILE(759) * 1.1574074074074112e-08;
    //k_P11_d1, mw97534dda_dfc5_4793_b791_a591b6a8f181, index: 415
    //Unit: metre^(-3)mole^(1)
    _class_parameter[415] = PFILE(760) * 999.9999999999994;
    //p11_50, mw705343cb_489c_42f0_96d7_f6aaa762a63e, index: 416
    //Unit: metre^(-2)mole^(1)
    _class_parameter[416] = PFILE(761) * 1.66053872801495e-12;
    //P11_C1, mwca81dc85_9901_41a0_a1cd_05971f1008d9, index: 417
    //Unit: metre^(-3)
    _class_parameter[417] = PFILE(762) * 6.022141989999979e+26;
    //k_M1p11_TCR_on, mwca556706_7cd8_47fb_8ea5_52073d58fcb2, index: 418
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[418] = PFILE(763) * 602214199000.0;
    //k_M1p11_TCR_off, mw8f596c91_fa9d_4bd9_99be_370d466e7641, index: 419
    //Unit: second^(-1)
    _class_parameter[419] = PFILE(764) * 1.0;
    //k_M1p11_TCR_p, mw14ecd973_0200_48e3_bc2a_5d0da6b555f0, index: 420
    //Unit: second^(-1)
    _class_parameter[420] = PFILE(765) * 1.0;
    //phi_M1p11_TCR, mwefe51424_f9f4_41a8_b76d_2f7a686d5422, index: 421
    //Unit: second^(-1)
    _class_parameter[421] = PFILE(766) * 1.0;
    //N_M1p11_TCR, mw319b9334_d8b6_4008_8a8c_7ff1f4e0814e, index: 422
    //Unit: dimensionless^(1)
    _class_parameter[422] = PFILE(767) * 1.0;
    //TCR_p11_tot, mwc7c3d658_251c_4a14_a252_eff031d1f1e1, index: 423
    //Unit: metre^(-2)mole^(1)
    _class_parameter[423] = PFILE(768) * 1.66053872801495e-12;
    //k_P12_up, mwf02f6e42_e700_4f92_bbb0_1ef0d1ca51a7, index: 424
    //Unit: mole^(-1)second^(-1)
    _class_parameter[424] = PFILE(770) * 6.97007174768519e+18;
    //k_xP12_deg, mwff3f65bc_ba7b_462c_b0ba_cd8b00566a88, index: 425
    //Unit: second^(-1)
    _class_parameter[425] = PFILE(771) * 1.15740740740741e-05;
    //k_P12_deg, mwe0c6e2c0_4b5e_4b12_9608_6b639af9770b, index: 426
    //Unit: second^(-1)
    _class_parameter[426] = PFILE(772) * 1.15740740740741e-05;
    //k_p12_deg, mw28667e35_bedc_447d_8fbc_38a4ae4f17dd, index: 427
    //Unit: second^(-1)
    _class_parameter[427] = PFILE(773) * 1.15740740740741e-05;
    //k_P12_on, mwca551c56_1f67_41a3_ba58_ba736fff4fce, index: 428
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[428] = PFILE(774) * 1.1574074074074112e-08;
    //k_P12_d1, mw64a2c3e6_5602_40fb_b786_a771a40549dc, index: 429
    //Unit: metre^(-3)mole^(1)
    _class_parameter[429] = PFILE(775) * 999.9999999999994;
    //p12_50, mw1ee3487d_83d3_4d2d_babb_61cd163fd927, index: 430
    //Unit: metre^(-2)mole^(1)
    _class_parameter[430] = PFILE(776) * 1.66053872801495e-12;
    //P12_C1, mw5bbfeb72_bd79_401c_8720_9487bfd56a23, index: 431
    //Unit: metre^(-3)
    _class_parameter[431] = PFILE(777) * 6.022141989999979e+26;
    //k_M1p12_TCR_on, mw75cef1a6_95ef_4ec5_b640_51fdeed860f0, index: 432
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[432] = PFILE(778) * 602214199000.0;
    //k_M1p12_TCR_off, mw4a87490f_b111_498b_8839_eb44076f76ed, index: 433
    //Unit: second^(-1)
    _class_parameter[433] = PFILE(779) * 1.0;
    //k_M1p12_TCR_p, mw0f8c4eb5_56ea_4840_8b0b_3e0ee065140c, index: 434
    //Unit: second^(-1)
    _class_parameter[434] = PFILE(780) * 1.0;
    //phi_M1p12_TCR, mw8814ed67_2164_48e6_87ed_1072a6124768, index: 435
    //Unit: second^(-1)
    _class_parameter[435] = PFILE(781) * 1.0;
    //N_M1p12_TCR, mw5fa761a3_c9db_4007_a66d_f95dfd0cabd2, index: 436
    //Unit: dimensionless^(1)
    _class_parameter[436] = PFILE(782) * 1.0;
    //TCR_p12_tot, mw46f1dbb3_ed5b_4a6e_99e2_336f4fb790e9, index: 437
    //Unit: metre^(-2)mole^(1)
    _class_parameter[437] = PFILE(783) * 1.66053872801495e-12;
    //k_P13_up, mw94c90f39_e05f_4c57_9a0c_6939751ac042, index: 438
    //Unit: mole^(-1)second^(-1)
    _class_parameter[438] = PFILE(785) * 6.97007174768519e+18;
    //k_xP13_deg, mwb396992e_a8a9_4035_8f1d_527801535ca8, index: 439
    //Unit: second^(-1)
    _class_parameter[439] = PFILE(786) * 1.15740740740741e-05;
    //k_P13_deg, mw7d190164_9290_4e7e_af3e_501bb9c343aa, index: 440
    //Unit: second^(-1)
    _class_parameter[440] = PFILE(787) * 1.15740740740741e-05;
    //k_p13_deg, mw48694d83_5bfc_4532_b83c_6dc35df7e929, index: 441
    //Unit: second^(-1)
    _class_parameter[441] = PFILE(788) * 1.15740740740741e-05;
    //k_P13_on, mw93646a12_d5af_448b_bd24_4cce0d7bfda8, index: 442
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[442] = PFILE(789) * 1.1574074074074112e-08;
    //k_P13_d1, mw6adec3d9_5a6f_476e_bee0_41e46b2ee16a, index: 443
    //Unit: metre^(-3)mole^(1)
    _class_parameter[443] = PFILE(790) * 999.9999999999994;
    //p13_50, mw6657664b_c642_4dc5_8fb5_c2f820ac87fd, index: 444
    //Unit: metre^(-2)mole^(1)
    _class_parameter[444] = PFILE(791) * 1.66053872801495e-12;
    //P13_C1, mwaadec916_5b95_4a99_8e13_c74890cb59b9, index: 445
    //Unit: metre^(-3)
    _class_parameter[445] = PFILE(792) * 6.022141989999979e+26;
    //k_M1p13_TCR_on, mwd505f7ec_19ba_44f3_9fd6_c430bc57cc41, index: 446
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[446] = PFILE(793) * 602214199000.0;
    //k_M1p13_TCR_off, mw5567727b_546a_498a_ad80_595c5ae93527, index: 447
    //Unit: second^(-1)
    _class_parameter[447] = PFILE(794) * 1.0;
    //k_M1p13_TCR_p, mwef84df07_d490_4170_b199_65164a8cdc3c, index: 448
    //Unit: second^(-1)
    _class_parameter[448] = PFILE(795) * 1.0;
    //phi_M1p13_TCR, mwf397ee2f_cb79_4056_91ec_7e5284cac082, index: 449
    //Unit: second^(-1)
    _class_parameter[449] = PFILE(796) * 1.0;
    //N_M1p13_TCR, mwceaa6435_ceb6_4c25_99b0_4338aebad9d9, index: 450
    //Unit: dimensionless^(1)
    _class_parameter[450] = PFILE(797) * 1.0;
    //TCR_p13_tot, mwed7b661c_be02_4d0a_8249_6646d33f1ec8, index: 451
    //Unit: metre^(-2)mole^(1)
    _class_parameter[451] = PFILE(798) * 1.66053872801495e-12;
    //k_P14_up, mwc8a722dd_0b59_414a_8c63_18f932145f46, index: 452
    //Unit: mole^(-1)second^(-1)
    _class_parameter[452] = PFILE(800) * 6.97007174768519e+18;
    //k_xP14_deg, mwc6fcbf4a_b2a8_4e60_8959_d479483f125a, index: 453
    //Unit: second^(-1)
    _class_parameter[453] = PFILE(801) * 1.15740740740741e-05;
    //k_P14_deg, mw64345816_5fac_453c_b0bc_6311330e4c29, index: 454
    //Unit: second^(-1)
    _class_parameter[454] = PFILE(802) * 1.15740740740741e-05;
    //k_p14_deg, mwbe2b83d3_3fee_4798_b4f4_2c5eb8886256, index: 455
    //Unit: second^(-1)
    _class_parameter[455] = PFILE(803) * 1.15740740740741e-05;
    //k_P14_on, mw83420379_e776_4b7c_ae12_dd2a590dd123, index: 456
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[456] = PFILE(804) * 1.1574074074074112e-08;
    //k_P14_d1, mw4ccc5be3_c69d_4f03_a28b_77fac808d766, index: 457
    //Unit: metre^(-3)mole^(1)
    _class_parameter[457] = PFILE(805) * 999.9999999999994;
    //p14_50, mw5009d0e0_65a3_4997_b511_3b26d7f6b778, index: 458
    //Unit: metre^(-2)mole^(1)
    _class_parameter[458] = PFILE(806) * 1.66053872801495e-12;
    //P14_C1, mwf95a4c8b_1e6c_4616_9f0d_60fbb0c9c7f5, index: 459
    //Unit: metre^(-3)
    _class_parameter[459] = PFILE(807) * 6.022141989999979e+26;
    //k_M1p14_TCR_on, mw7aa296f9_c5b0_418e_a00d_db49748733cb, index: 460
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[460] = PFILE(808) * 602214199000.0;
    //k_M1p14_TCR_off, mw3581efb3_de82_4f57_a719_89ae2fd4e2c0, index: 461
    //Unit: second^(-1)
    _class_parameter[461] = PFILE(809) * 1.0;
    //k_M1p14_TCR_p, mw48b06332_c671_4ab0_9bab_71cb8924d167, index: 462
    //Unit: second^(-1)
    _class_parameter[462] = PFILE(810) * 1.0;
    //phi_M1p14_TCR, mwf3d28581_a3bc_44e2_9f1e_69fed5203d5a, index: 463
    //Unit: second^(-1)
    _class_parameter[463] = PFILE(811) * 1.0;
    //N_M1p14_TCR, mw595e7bed_e3d4_467d_8c2f_487efe096f92, index: 464
    //Unit: dimensionless^(1)
    _class_parameter[464] = PFILE(812) * 1.0;
    //TCR_p14_tot, mw67223784_f514_4d05_8f1d_7ad3f2410318, index: 465
    //Unit: metre^(-2)mole^(1)
    _class_parameter[465] = PFILE(813) * 1.66053872801495e-12;
    //k_P15_up, mweb82057d_139c_485a_9aab_db5f9c516252, index: 466
    //Unit: mole^(-1)second^(-1)
    _class_parameter[466] = PFILE(815) * 6.97007174768519e+18;
    //k_xP15_deg, mw801df6d9_b5ac_4bfc_b157_83a61575c3d1, index: 467
    //Unit: second^(-1)
    _class_parameter[467] = PFILE(816) * 1.15740740740741e-05;
    //k_P15_deg, mwe04220d8_93ed_4b11_b944_e144e191cbfc, index: 468
    //Unit: second^(-1)
    _class_parameter[468] = PFILE(817) * 1.15740740740741e-05;
    //k_p15_deg, mwfcfc5444_3574_46d6_97ae_a24ad16e05fc, index: 469
    //Unit: second^(-1)
    _class_parameter[469] = PFILE(818) * 1.15740740740741e-05;
    //k_P15_on, mw3d983db9_2ffd_43c4_8359_c21427604436, index: 470
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[470] = PFILE(819) * 1.1574074074074112e-08;
    //k_P15_d1, mwb6240d35_5621_4d0c_9004_da287e692f3e, index: 471
    //Unit: metre^(-3)mole^(1)
    _class_parameter[471] = PFILE(820) * 999.9999999999994;
    //p15_50, mw8a24a5c3_185b_402d_9538_8e589ecaf2a2, index: 472
    //Unit: metre^(-2)mole^(1)
    _class_parameter[472] = PFILE(821) * 1.66053872801495e-12;
    //P15_C1, mweaf7f189_26ff_46db_a106_b4f43ae3275a, index: 473
    //Unit: metre^(-3)
    _class_parameter[473] = PFILE(822) * 6.022141989999979e+26;
    //k_M1p15_TCR_on, mw929f00e8_72c8_46af_a540_8ddc34bee5d2, index: 474
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[474] = PFILE(823) * 602214199000.0;
    //k_M1p15_TCR_off, mwcdd6b8d0_f3ad_446f_8a02_64d2c22c15ea, index: 475
    //Unit: second^(-1)
    _class_parameter[475] = PFILE(824) * 1.0;
    //k_M1p15_TCR_p, mwe8c7fcf7_cf0e_4d34_ab7e_75dd98f0c021, index: 476
    //Unit: second^(-1)
    _class_parameter[476] = PFILE(825) * 1.0;
    //phi_M1p15_TCR, mwc390b9d7_d69d_4ace_a1f1_04a1c24fdc63, index: 477
    //Unit: second^(-1)
    _class_parameter[477] = PFILE(826) * 1.0;
    //N_M1p15_TCR, mwb06b1967_aa9a_425a_9dcb_a5bfa946d626, index: 478
    //Unit: dimensionless^(1)
    _class_parameter[478] = PFILE(827) * 1.0;
    //TCR_p15_tot, mw3f5b3e03_c295_4525_b8df_310f334bdbc1, index: 479
    //Unit: metre^(-2)mole^(1)
    _class_parameter[479] = PFILE(828) * 1.66053872801495e-12;
    //k_P16_up, mw785c215d_650e_4de0_98e3_7d78eddc5381, index: 480
    //Unit: mole^(-1)second^(-1)
    _class_parameter[480] = PFILE(830) * 6.97007174768519e+18;
    //k_xP16_deg, mw44e2a0c8_0354_4267_9ad6_357d622fd383, index: 481
    //Unit: second^(-1)
    _class_parameter[481] = PFILE(831) * 1.15740740740741e-05;
    //k_P16_deg, mwde5bc21e_e781_4149_ab97_0900de467d94, index: 482
    //Unit: second^(-1)
    _class_parameter[482] = PFILE(832) * 1.15740740740741e-05;
    //k_p16_deg, mwe12d4631_1675_4f53_b485_1087766a1b13, index: 483
    //Unit: second^(-1)
    _class_parameter[483] = PFILE(833) * 1.15740740740741e-05;
    //k_P16_on, mwf73894d9_765a_4ba0_856d_abc0917f2208, index: 484
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[484] = PFILE(834) * 1.1574074074074112e-08;
    //k_P16_d1, mwb59295c3_5d56_4ab7_bc28_a2568f6af94a, index: 485
    //Unit: metre^(-3)mole^(1)
    _class_parameter[485] = PFILE(835) * 999.9999999999994;
    //p16_50, mwfabe78ef_43c1_4628_9516_3486fe336bb2, index: 486
    //Unit: metre^(-2)mole^(1)
    _class_parameter[486] = PFILE(836) * 1.66053872801495e-12;
    //P16_C1, mwf920f490_cfb6_4e45_84d9_ec752261f585, index: 487
    //Unit: metre^(-3)
    _class_parameter[487] = PFILE(837) * 6.022141989999979e+26;
    //k_M1p16_TCR_on, mwbad1f034_c8a7_42f2_92cc_e3588bbe0d53, index: 488
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[488] = PFILE(838) * 602214199000.0;
    //k_M1p16_TCR_off, mwe6934db8_b459_41c2_a709_dc9676de7e9e, index: 489
    //Unit: second^(-1)
    _class_parameter[489] = PFILE(839) * 1.0;
    //k_M1p16_TCR_p, mw562fb8c9_7bce_4da7_b16b_61a588ac6fea, index: 490
    //Unit: second^(-1)
    _class_parameter[490] = PFILE(840) * 1.0;
    //phi_M1p16_TCR, mwc05f4a52_c1e2_4b64_a8e5_ba8cbe669760, index: 491
    //Unit: second^(-1)
    _class_parameter[491] = PFILE(841) * 1.0;
    //N_M1p16_TCR, mw2d2866b8_ca7b_4be7_a339_7ca6c9c610c6, index: 492
    //Unit: dimensionless^(1)
    _class_parameter[492] = PFILE(842) * 1.0;
    //TCR_p16_tot, mwc48d3b4c_6eba_4f77_ac79_20d5c5f3324f, index: 493
    //Unit: metre^(-2)mole^(1)
    _class_parameter[493] = PFILE(843) * 1.66053872801495e-12;
    //k_P17_up, mw00917373_dfb8_475a_9b82_a8d78d81bb40, index: 494
    //Unit: mole^(-1)second^(-1)
    _class_parameter[494] = PFILE(845) * 6.97007174768519e+18;
    //k_xP17_deg, mw71cb3f73_73d8_49d2_9928_5254785b6519, index: 495
    //Unit: second^(-1)
    _class_parameter[495] = PFILE(846) * 1.15740740740741e-05;
    //k_P17_deg, mwab138d35_a73c_462d_b97e_59f98c547f58, index: 496
    //Unit: second^(-1)
    _class_parameter[496] = PFILE(847) * 1.15740740740741e-05;
    //k_p17_deg, mw2d739c36_7a5f_47d5_bd2e_f8e514f2c739, index: 497
    //Unit: second^(-1)
    _class_parameter[497] = PFILE(848) * 1.15740740740741e-05;
    //k_P17_on, mwcc420fe1_259b_4d95_966f_91fc42d698ca, index: 498
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[498] = PFILE(849) * 1.1574074074074112e-08;
    //k_P17_d1, mwdba663a7_86d3_402d_8cf6_d7ce67c56409, index: 499
    //Unit: metre^(-3)mole^(1)
    _class_parameter[499] = PFILE(850) * 999.9999999999994;
    //p17_50, mw55b6b5fd_5084_42bf_835e_fdc52404527b, index: 500
    //Unit: metre^(-2)mole^(1)
    _class_parameter[500] = PFILE(851) * 1.66053872801495e-12;
    //P17_C1, mwd5b8eaff_5e27_4e30_9ed4_38206d77fab4, index: 501
    //Unit: metre^(-3)
    _class_parameter[501] = PFILE(852) * 6.022141989999979e+26;
    //k_M1p17_TCR_on, mw09c0a184_1436_4ce6_bc43_a0ab61f9e654, index: 502
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[502] = PFILE(853) * 602214199000.0;
    //k_M1p17_TCR_off, mw042e564d_3e8e_44de_a67d_823e364784fb, index: 503
    //Unit: second^(-1)
    _class_parameter[503] = PFILE(854) * 1.0;
    //k_M1p17_TCR_p, mw9aa08f40_b15c_481d_b11b_0503a65eab4a, index: 504
    //Unit: second^(-1)
    _class_parameter[504] = PFILE(855) * 1.0;
    //phi_M1p17_TCR, mw75eac6f9_1075_49ed_bd53_334d0a84b3d2, index: 505
    //Unit: second^(-1)
    _class_parameter[505] = PFILE(856) * 1.0;
    //N_M1p17_TCR, mw7140ebf3_7fa3_4965_af73_3f8648b15dd1, index: 506
    //Unit: dimensionless^(1)
    _class_parameter[506] = PFILE(857) * 1.0;
    //TCR_p17_tot, mw00a6ca38_7765_4fc8_a1d2_8c14f1b2b619, index: 507
    //Unit: metre^(-2)mole^(1)
    _class_parameter[507] = PFILE(858) * 1.66053872801495e-12;
    //kon_PD1_PDL1, mwb8181402_42e4_4f2e_b1b2_4e0ba8c1d855, index: 508
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[508] = PFILE(860) * 1000000000000.0;
    //q_P_nivo, mw1118de81_ad78_4042_abf7_d58c98041ce1, index: 509
    //Unit: metre^(3)second^(-1)
    _class_parameter[509] = PFILE(861) * 0.0010000000000000007;
    //q_T_nivo, mw9c6c8627_997b_41bf_a722_2e6d6a585034, index: 510
    //Unit: metre^(3)second^(-1)
    _class_parameter[510] = PFILE(862) * 1.0000000000000006e-06;
    //q_LN_nivo, mwdd2499c1_8a0a_425a_ac1b_ace7ffcd0b2f, index: 511
    //Unit: metre^(3)second^(-1)
    _class_parameter[511] = PFILE(863) * 1.0000000000000006e-06;
    //q_LD_nivo, mwfc68aa9c_9309_4fc1_8e69_2577d858282a, index: 512
    //Unit: second^(-1)
    _class_parameter[512] = PFILE(864) * 0.0166666666666667;
    //k_cl_nivo, mw51049dbe_a8d9_4107_ac95_85cd28efbf64, index: 513
    //Unit: metre^(3)second^(-1)
    _class_parameter[513] = PFILE(865) * 1.1574074074074112e-08;
    //gamma_C_nivo, mw2ea82991_2fd6_4143_96e2_8335857a821f, index: 514
    //Unit: dimensionless^(1)
    _class_parameter[514] = PFILE(866) * 1.0;
    //gamma_P_nivo, mwcd6558c0_5df7_476e_9164_0e75c5728c8c, index: 515
    //Unit: dimensionless^(1)
    _class_parameter[515] = PFILE(867) * 1.0;
    //gamma_T_nivo, mwfac548f7_da13_41c1_9b6b_359a2b905440, index: 516
    //Unit: dimensionless^(1)
    _class_parameter[516] = PFILE(868) * 1.0;
    //gamma_LN_nivo, mw9da87375_4b0f_428f_ab78_76aa221d8ef2, index: 517
    //Unit: dimensionless^(1)
    _class_parameter[517] = PFILE(869) * 1.0;
    //q_P_durv, mwbad92b96_48d9_4297_aa06_04da85f192ab, index: 518
    //Unit: metre^(3)second^(-1)
    _class_parameter[518] = PFILE(870) * 0.0010000000000000007;
    //q_T_durv, mwe18adf59_40b7_417b_b506_c050b792c558, index: 519
    //Unit: metre^(3)second^(-1)
    _class_parameter[519] = PFILE(871) * 1.0000000000000006e-06;
    //q_LN_durv, mwff4d3221_c556_42b9_ad84_98d3cad02a37, index: 520
    //Unit: metre^(3)second^(-1)
    _class_parameter[520] = PFILE(872) * 1.0000000000000006e-06;
    //q_LD_durv, mwa3734600_66f6_4f18_b6ce_ff2354a9ec0c, index: 521
    //Unit: second^(-1)
    _class_parameter[521] = PFILE(873) * 0.0166666666666667;
    //k_cl_durv, mwb2777d3a_7a08_473c_81c7_757b7e0d7cf6, index: 522
    //Unit: metre^(3)second^(-1)
    _class_parameter[522] = PFILE(874) * 1.1574074074074112e-08;
    //gamma_C_durv, mwc4e956e0_2e0a_4168_9c79_a3662657fe25, index: 523
    //Unit: dimensionless^(1)
    _class_parameter[523] = PFILE(875) * 1.0;
    //gamma_P_durv, mwa58f74d2_ba43_4248_ac61_35414aaf93bb, index: 524
    //Unit: dimensionless^(1)
    _class_parameter[524] = PFILE(876) * 1.0;
    //gamma_T_durv, mw99231c1c_3d32_47da_a2d3_2be3cba92e24, index: 525
    //Unit: dimensionless^(1)
    _class_parameter[525] = PFILE(877) * 1.0;
    //gamma_LN_durv, mw4f713a5b_69d9_4d59_91d1_9edaca75d00c, index: 526
    //Unit: dimensionless^(1)
    _class_parameter[526] = PFILE(878) * 1.0;
    //q_P_ipi, mw651f7937_756e_4fd1_92f8_de3831981346, index: 527
    //Unit: metre^(3)second^(-1)
    _class_parameter[527] = PFILE(879) * 0.0010000000000000007;
    //q_T_ipi, mw66e1f287_db96_4ee3_b628_e447e86d182d, index: 528
    //Unit: metre^(3)second^(-1)
    _class_parameter[528] = PFILE(880) * 1.0000000000000006e-06;
    //q_LN_ipi, mwc249bcdb_8833_40af_9438_1fef554cbf84, index: 529
    //Unit: metre^(3)second^(-1)
    _class_parameter[529] = PFILE(881) * 1.0000000000000006e-06;
    //q_LD_ipi, mwac7a3ffb_1b1c_45c2_b8a2_68bcdf74a5c5, index: 530
    //Unit: second^(-1)
    _class_parameter[530] = PFILE(882) * 0.0166666666666667;
    //k_cl_ipi, mw3bf3fc49_1cf2_4590_911a_176b17cc813d, index: 531
    //Unit: metre^(3)second^(-1)
    _class_parameter[531] = PFILE(883) * 1.1574074074074112e-08;
    //gamma_C_ipi, mw2061ce87_5261_4512_bd15_5e33a644a0e5, index: 532
    //Unit: dimensionless^(1)
    _class_parameter[532] = PFILE(884) * 1.0;
    //gamma_P_ipi, mwd40ea5b1_0921_4225_9114_d3848c8aff4d, index: 533
    //Unit: dimensionless^(1)
    _class_parameter[533] = PFILE(885) * 1.0;
    //gamma_T_ipi, mwf16081dd_6c66_4557_b6bc_df226c5d5f79, index: 534
    //Unit: dimensionless^(1)
    _class_parameter[534] = PFILE(886) * 1.0;
    //gamma_LN_ipi, mwe2323f7f_0329_4450_9c47_0f97b85623d5, index: 535
    //Unit: dimensionless^(1)
    _class_parameter[535] = PFILE(887) * 1.0;
    //kon_PD1_PDL2, mw622bfd23_dbe3_4106_b26c_90050f8cb5dc, index: 536
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[536] = PFILE(888) * 1000000000000.0;
    //kon_PD1_nivo, mw1d0c74e2_ceac_48b7_a000_df4c73abeb0d, index: 537
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[537] = PFILE(889) * 0.0010000000000000007;
    //kon_PDL1_durv, mw08edba64_9ba2_42c4_98f5_38582e096cc2, index: 538
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[538] = PFILE(890) * 0.0010000000000000007;
    //kon_CD28_CD80, mwb30bb81a_8438_400e_b6b9_fd6403bdc348, index: 539
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[539] = PFILE(891) * 1000000000000.0;
    //kon_CD28_CD86, mw0802c7ea_5528_4b71_8ee3_d4f7ce1f29cc, index: 540
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[540] = PFILE(892) * 1000000000000.0;
    //kon_CTLA4_CD80, mw30c478fe_b827_4e66_93f3_374d12945edd, index: 541
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[541] = PFILE(893) * 1000000000000.0;
    //kon_CTLA4_CD86, mw6349e81a_eb06_4700_8a9c_2c87c34af8e5, index: 542
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[542] = PFILE(894) * 1000000000000.0;
    //kon_CD80_PDL1, mw52e60d52_5d60_4654_8cdf_f5edaac28d6c, index: 543
    //Unit: metre^(2)mole^(-1)second^(-1)
    _class_parameter[543] = PFILE(895) * 1000000000000.0;
    //kon_CTLA4_ipi, mw0dbf25a3_b47a_483c_891d_01536da74974, index: 544
    //Unit: metre^(3)mole^(-1)second^(-1)
    _class_parameter[544] = PFILE(896) * 0.0010000000000000007;
    //koff_PD1_PDL1, mwe34c46f4_b4d7_435f_b8b7_2273b12dd86b, index: 545
    //Unit: second^(-1)
    _class_parameter[545] = PFILE(897) * 1.0;
    //koff_PD1_PDL2, mwe0533580_3920_481f_98e0_d2e1a28c288e, index: 546
    //Unit: second^(-1)
    _class_parameter[546] = PFILE(898) * 1.0;
    //koff_PD1_nivo, mw45a25844_601f_4d85_a0cc_bf65dc8510f3, index: 547
    //Unit: second^(-1)
    _class_parameter[547] = PFILE(899) * 1.0;
    //koff_PDL1_durv, mwb75a39b2_9ac4_4dfe_a2ab_1bbb98206892, index: 548
    //Unit: second^(-1)
    _class_parameter[548] = PFILE(900) * 1.0;
    //koff_CD28_CD80, mwc204f3d5_2f78_43f0_bbdb_e8036dabf823, index: 549
    //Unit: second^(-1)
    _class_parameter[549] = PFILE(901) * 1.0;
    //koff_CD28_CD86, mw1489846f_cca6_49aa_a08f_ef9da867d083, index: 550
    //Unit: second^(-1)
    _class_parameter[550] = PFILE(902) * 1.0;
    //koff_CTLA4_CD80, mwce6af464_d45e_4e55_ad02_31232d7563be, index: 551
    //Unit: second^(-1)
    _class_parameter[551] = PFILE(903) * 1.0;
    //koff_CTLA4_CD86, mw6ef1cb98_a5c3_40e4_8854_39352f4c3238, index: 552
    //Unit: second^(-1)
    _class_parameter[552] = PFILE(904) * 1.0;
    //koff_CD80_PDL1, mwc734bfd3_9ebb_4281_9bdc_289720c5b85a, index: 553
    //Unit: second^(-1)
    _class_parameter[553] = PFILE(905) * 1.0;
    //koff_CTLA4_ipi, mwee92919b_6672_4ccf_9da4_f53c08b94cdb, index: 554
    //Unit: second^(-1)
    _class_parameter[554] = PFILE(906) * 1.0;
    //Chi_PD1_nivo, mw635f93b7_1afc_4140_b2b5_375e9c1b3a85, index: 555
    //Unit: metre^(-1)
    _class_parameter[555] = PFILE(907) * 999999999.9999999;
    //Chi_PDL1_durv, mw2060d5b0_badc_4503_b9c2_92e85678ed39, index: 556
    //Unit: metre^(-1)
    _class_parameter[556] = PFILE(908) * 999999999.9999999;
    //Chi_CTLA4_ipi, mw14e1fdfd_df86_4d44_b4c8_fb0c75796675, index: 557
    //Unit: metre^(-1)
    _class_parameter[557] = PFILE(909) * 999999999.9999999;
    //PD1_50, mw70b7fe16_006d_4d57_b9d4_297b69e9bea5, index: 558
    //Unit: metre^(-2)mole^(1)
    _class_parameter[558] = PFILE(910) * 1.66053872801495e-12;
    //n_PD1, mwcf57eb97_dd1c_4fcb_b7bd_7fa469d44f7f, index: 559
    //Unit: dimensionless^(1)
    _class_parameter[559] = PFILE(911) * 1.0;
    //CD28_CD8X_50, mwc0074702_53a3_4182_96e8_4cde8df5719f, index: 560
    //Unit: metre^(-2)mole^(1)
    _class_parameter[560] = PFILE(912) * 1.66053872801495e-12;
    //n_CD28_CD8X, mw61ac6690_7201_4ba4_83e9_a47f27571bef, index: 561
    //Unit: dimensionless^(1)
    _class_parameter[561] = PFILE(913) * 1.0;
    //T_PD1_total, mwd71a1f52_01a1_474d_93ec_797c8eb06faa, index: 562
    //Unit: mole^(1)
    _class_parameter[562] = PFILE(914) * 1.66053872801495e-24;
    //T_CD28_total, mw9cf7847c_375b_4b8c_92bf_c4cbfaec77b5, index: 563
    //Unit: mole^(1)
    _class_parameter[563] = PFILE(915) * 1.66053872801495e-24;
    //T_CTLA4_syn, mw664b1465_1717_4523_84d5_d62579e3392b, index: 564
    //Unit: mole^(1)
    _class_parameter[564] = PFILE(916) * 1.66053872801495e-24;
    //T_PDL1_total, mw3112e554_66c8_40a4_8d8a_92c09b3702c1, index: 565
    //Unit: mole^(1)
    _class_parameter[565] = PFILE(917) * 1.66053872801495e-24;
    //C1_PDL1_total, mw37dd9b9c_aa2b_4723_bbd9_e312895a27f8, index: 566
    //Unit: mole^(1)
    _class_parameter[566] = PFILE(918) * 1.66053872801495e-24;
    //C1_PDL2_total, mwd5171c86_7e30_44d6_9d07_134fab827402, index: 567
    //Unit: mole^(1)
    _class_parameter[567] = PFILE(919) * 1.66053872801495e-24;
    //C1_CD80_total, mwb17cfb4f_5730_4c3f_8635_f63cc43d68df, index: 568
    //Unit: mole^(1)
    _class_parameter[568] = PFILE(920) * 1.66053872801495e-24;
    //C1_CD86_total, mwdc4cfa19_6f32_447b_84d8_9baf0e7ef894, index: 569
    //Unit: mole^(1)
    _class_parameter[569] = PFILE(921) * 1.66053872801495e-24;
    //APC_PDL1_total, mw9a660315_7235_41e2_921f_71a63e456a56, index: 570
    //Unit: mole^(1)
    _class_parameter[570] = PFILE(922) * 1.66053872801495e-24;
    //APC_PDL2_total, mwdf170272_5222_4a9e_a39f_95b6fe574bb1, index: 571
    //Unit: mole^(1)
    _class_parameter[571] = PFILE(923) * 1.66053872801495e-24;
    //APC_CD80_total, mwbf1eee9f_0e5b_485c_8576_d62b6ae5f799, index: 572
    //Unit: mole^(1)
    _class_parameter[572] = PFILE(924) * 1.66053872801495e-24;
    //APC_CD86_total, mw1fb6111c_b4ac_48c7_a88f_37febd247658, index: 573
    //Unit: mole^(1)
    _class_parameter[573] = PFILE(925) * 1.66053872801495e-24;
    //Treg_CTLA4_tot, mwa99e9620_9835_4936_8dac_fbeb1f90e5b9, index: 574
    //Unit: mole^(1)
    _class_parameter[574] = PFILE(926) * 1.66053872801495e-24;
    //Treg_CTLA4_50, mw9a40a96a_9c56_44a8_96fe_8d4f0e930d87, index: 575
    //Unit: mole^(1)
    _class_parameter[575] = PFILE(927) * 1.66053872801495e-24;
    //n_Treg_CTLA4, mwf65307d8_2656_47f4_9d4b_8c3d5252d18b, index: 576
    //Unit: dimensionless^(1)
    _class_parameter[576] = PFILE(928) * 1.0;
    //k_CTLA4_ADCC, mwed062cf1_6892_42d4_a374_a5e1c506fc3e, index: 577
    //Unit: second^(-1)
    _class_parameter[577] = PFILE(929) * 1.15740740740741e-05;
    //k_rec_MDSC, mwdf185fcf_61e7_4aac_9a2a_525eab2ebe5e, index: 578
    //Unit: second^(-1)
    _class_parameter[578] = PFILE(932) * 1.15740740740741e-05;
    //kd_MDSC, mw2f45c6d2_a45a_4db9_8e59_0d5c44b2b7cb, index: 579
    //Unit: second^(-1)
    _class_parameter[579] = PFILE(933) * 1.15740740740741e-05;
    //IC50_ENT_C, mw192a0047_8b54_4c25_9aa0_8f250985a801, index: 580
    //Unit: metre^(-3)mole^(1)
    _class_parameter[580] = PFILE(934) * 999.9999999999994;
    //k_deg_CCL2, mwd2dbd533_a0fe_4230_9d6b_a568191b00fc, index: 581
    //Unit: second^(-1)
    _class_parameter[581] = PFILE(935) * 0.000277777777777778;
    //k_deg_NO, mw4a8a9125_38c7_4e14_9b8a_4cf5bb48413c, index: 582
    //Unit: second^(-1)
    _class_parameter[582] = PFILE(936) * 1.15740740740741e-05;
    //k_deg_ArgI, mwb668761b_199e_4a89_8558_2366c09f4740, index: 583
    //Unit: second^(-1)
    _class_parameter[583] = PFILE(937) * 1.15740740740741e-05;
    //k_sec_CCL2, mw0343546d_b1d5_4901_9f49_da54fd8d5a83, index: 584
    //Unit: second^(-1)
    _class_parameter[584] = PFILE(938) * 6970071747.68519;
    //k_sec_NO, mw7587e62a_3cab_4665_942e_a1e288bcca44, index: 585
    //Unit: second^(-1)
    _class_parameter[585] = PFILE(939) * 6970071747.68519;
    //k_sec_ArgI, mwa6ab5b6f_4866_4760_9af5_dc3b1ab42197, index: 586
    //Unit: second^(-1)
    _class_parameter[586] = PFILE(940) * 6970071747685.19;
    //IC50_ENT_NO, mw99d7bc72_df64_4e61_bc2e_d5f19e945ab0, index: 587
    //Unit: metre^(-3)mole^(1)
    _class_parameter[587] = PFILE(941) * 999.9999999999994;
    //ki_Treg, mw871ff489_c07e_4b9f_864a_1565a80c29d1, index: 588
    //Unit: second^(-1)
    _class_parameter[588] = PFILE(942) * 1.15740740740741e-05;
    //IC50_ArgI_CTL, mw9d51a770_02b9_4d85_b8f4_c587c606d7d6, index: 589
    //Unit: metre^(-3)mole^(1)
    _class_parameter[589] = PFILE(943) * 999.9999999999994;
    //IC50_NO_CTL, mw998bcba9_c10d_4c20_8eda_a09b6dbae42b, index: 590
    //Unit: metre^(-3)mole^(1)
    _class_parameter[590] = PFILE(944) * 999.9999999999994;
    //EC50_CCL2_rec, mw15c5925c_0d40_44e8_b79c_fc27f36f76b2, index: 591
    //Unit: metre^(-3)mole^(1)
    _class_parameter[591] = PFILE(945) * 999.9999999999994;
    //EC50_ArgI_Treg, mweb363f06_8bfe_4ca0_a61e_42bd6185398f, index: 592
    //Unit: metre^(-3)mole^(1)
    _class_parameter[592] = PFILE(946) * 999.9999999999994;
    //MDSC_max, mw008468c2_2b9f_4397_9358_7752867acaa1, index: 593
    //Unit: metre^(-3)mole^(1)
    _class_parameter[593] = PFILE(947) * 1.6605387280149534e-18;
    //Treg_max, mwbb374fb9_1d4f_422f_b242_3302378b483c, index: 594
    //Unit: metre^(-3)mole^(1)
    _class_parameter[594] = PFILE(948) * 1.6605387280149534e-18;
    //IC50_ENT_CCL2, mw6e4ce77b_fd74_4bf6_976b_b98624c23b37, index: 595
    //Unit: metre^(-3)mole^(1)
    _class_parameter[595] = PFILE(949) * 999.9999999999994;
    //k_brec_MDSC, mwc26d2a88_0c6c_43be_904d_063832113fea, index: 596
    //Unit: second^(-1)
    _class_parameter[596] = PFILE(950) * 1.15740740740741e-05;
    //IC50_ENT_ArgI, mw3ed52353_7bd1_4f33_bed7_680d261e6446, index: 597
    //Unit: metre^(-3)mole^(1)
    _class_parameter[597] = PFILE(951) * 999.9999999999994;
    //k_a1_ENT, mw7b740318_7519_425d_9cac_d1e09443c510, index: 598
    //Unit: second^(-1)
    _class_parameter[598] = PFILE(952) * 0.000277777777777778;
    //k_a2_ENT, mwd33e4ffc_b590_44bd_a41a_46983c988592, index: 599
    //Unit: second^(-1)
    _class_parameter[599] = PFILE(953) * 0.000277777777777778;
    //k_cln_ENT, mw06cf230e_b449_4bdd_81d1_8cbcf4c9dcd8, index: 600
    //Unit: metre^(-3)mole^(1)second^(-1)
    _class_parameter[600] = PFILE(954) * 0.277777777777778;
    //Kc_ENT, mwa27e250f_3e51_46d6_b877_4705296e7108, index: 601
    //Unit: metre^(-3)mole^(1)
    _class_parameter[601] = PFILE(955) * 999.9999999999994;
    //lagP, mw17dbc3e7_70fa_4b90_9e78_c96fc1cad353, index: 602
    //Unit: second^(1)
    _class_parameter[602] = PFILE(956) * 3600.0;
    //durP, mw3ab724c5_1ba9_456c_bb30_38858b58a2f3, index: 603
    //Unit: second^(1)
    _class_parameter[603] = PFILE(957) * 3600.0;
    //k_dose2, mwcf83a00a_4a3f_41f4_9aef_221c30932dad, index: 604
    //Unit: second^(-1)
    _class_parameter[604] = PFILE(958) * 0.000277777777777778;
    //q_P_ENT, mw64ea5aae_9a89_4e9a_a42a_96f0120e0a53, index: 605
    //Unit: metre^(3)second^(-1)
    _class_parameter[605] = PFILE(959) * 1.0000000000000006e-06;
    //q_T_ENT, mw32374ece_a173_42c6_ae78_331f4f81a4aa, index: 606
    //Unit: metre^(3)second^(-1)
    _class_parameter[606] = PFILE(960) * 1.0000000000000006e-06;
    //q_LN_ENT, mw1bf8a6b5_80b1_4a32_9147_42be4bde05e6, index: 607
    //Unit: metre^(3)second^(-1)
    _class_parameter[607] = PFILE(961) * 1.0000000000000006e-06;
    //q_LD_ENT, mw9058340c_58e3_4d28_b100_8394c3844ea3, index: 608
    //Unit: second^(-1)
    _class_parameter[608] = PFILE(962) * 0.0166666666666667;
    //k_cl_ENT, mwd697706c_8c97_4d47_9a93_cf3786e14e81, index: 609
    //Unit: second^(-1)
    _class_parameter[609] = PFILE(963) * 0.000277777777777778;
    //gamma_C_ENT, mw175d86f5_e28e_4d39_9085_f13c4049e7fc, index: 610
    //Unit: dimensionless^(1)
    _class_parameter[610] = PFILE(964) * 1.0;
    //gamma_P_ENT, mw9913da12_68ae_4fd2_b717_bb2f0e6835a1, index: 611
    //Unit: dimensionless^(1)
    _class_parameter[611] = PFILE(965) * 1.0;
    //gamma_T_ENT, mw5eac4dec_f7b0_4c0e_b994_79dfd46aa1f7, index: 612
    //Unit: dimensionless^(1)
    _class_parameter[612] = PFILE(966) * 1.0;
    //gamma_LN_ENT, mw7c8dbea6_db4f_4620_b5b2_3be2c8fac8a9, index: 613
    //Unit: dimensionless^(1)
    _class_parameter[613] = PFILE(967) * 1.0;
}

void ODE_system::setupVariables(void){

    _species_var = std::vector<realtype>(296, 0);
    _nonspecies_var = std::vector<realtype>(0, 0);
    //species not part of ode left-hand side
    _species_other =  std::vector<realtype>(0, 0);
    
    return;
}


void ODE_system::setup_instance_varaibles(Param& param){

    //V_C.T0, mwaaaef7c6_e588_4d2c_bed7_47f96f061a2f, index: 0
    //Unit: mole^(1)
    _species_var[0] = PFILE(14) * 1.66053872801495e-24;
    //V_C.T1, mwcd3af0b0_d9fd_4a3a_af63_846a7422bb48, index: 1
    //Unit: mole^(1)
    _species_var[1] = PFILE(15) * 1.66053872801495e-24;
    //V_C.T2, mw8b3d85d3_0e27_471e_ba7e_daf2974a38b7, index: 2
    //Unit: mole^(1)
    _species_var[2] = PFILE(16) * 1.66053872801495e-24;
    //V_C.T3, mwe295babc_975c_451b_bb87_d516c6820a61, index: 3
    //Unit: mole^(1)
    _species_var[3] = PFILE(17) * 1.66053872801495e-24;
    //V_C.T4, mw0c2514b6_aea0_4f56_a825_409fe8443f64, index: 4
    //Unit: mole^(1)
    _species_var[4] = PFILE(18) * 1.66053872801495e-24;
    //V_C.T5, mwd595ad47_2b68_41eb_9a7b_0b62d916d9b5, index: 5
    //Unit: mole^(1)
    _species_var[5] = PFILE(19) * 1.66053872801495e-24;
    //V_C.T6, mw03e2e10d_ecb2_4fe0_a52e_2012023f0d55, index: 6
    //Unit: mole^(1)
    _species_var[6] = PFILE(20) * 1.66053872801495e-24;
    //V_C.T7, mw884d370b_81c2_4c28_9c08_71d44fa2da6c, index: 7
    //Unit: mole^(1)
    _species_var[7] = PFILE(21) * 1.66053872801495e-24;
    //V_C.T8, mw153c4d18_5905_4fe9_a533_e9a715802bf1, index: 8
    //Unit: mole^(1)
    _species_var[8] = PFILE(22) * 1.66053872801495e-24;
    //V_C.T9, mw3f217db0_e47b_4a40_a129_c52f7c173e1e, index: 9
    //Unit: mole^(1)
    _species_var[9] = PFILE(23) * 1.66053872801495e-24;
    //V_C.T10, mw21f6907b_0c36_402d_a4ac_d8daf5fe2f23, index: 10
    //Unit: mole^(1)
    _species_var[10] = PFILE(24) * 1.66053872801495e-24;
    //V_C.T11, mw3120f99d_839d_4b8b_996e_7f809ecdf57b, index: 11
    //Unit: mole^(1)
    _species_var[11] = PFILE(25) * 1.66053872801495e-24;
    //V_C.T12, mw8f628ccd_d582_4f5a_b5ff_606099edf059, index: 12
    //Unit: mole^(1)
    _species_var[12] = PFILE(26) * 1.66053872801495e-24;
    //V_C.T13, mw3dbd936e_d9d4_4a5c_b571_aea42c85e08e, index: 13
    //Unit: mole^(1)
    _species_var[13] = PFILE(27) * 1.66053872801495e-24;
    //V_C.T14, mw4752f771_499c_44ca_87e8_8ee91e080756, index: 14
    //Unit: mole^(1)
    _species_var[14] = PFILE(28) * 1.66053872801495e-24;
    //V_C.T15, mw663e3d1d_16cd_49ea_b983_43b1c44ea91c, index: 15
    //Unit: mole^(1)
    _species_var[15] = PFILE(29) * 1.66053872801495e-24;
    //V_C.T16, mw8e92e19f_abb1_4a34_bbba_9ffce8c42472, index: 16
    //Unit: mole^(1)
    _species_var[16] = PFILE(30) * 1.66053872801495e-24;
    //V_C.T17, mw22348f08_4711_4643_801e_ad0acfb0f87e, index: 17
    //Unit: mole^(1)
    _species_var[17] = PFILE(31) * 1.66053872801495e-24;
    //V_C.nivo, mw63c02d58_863e_4361_a077_e21b8d485948, index: 18
    //Unit: mole^(1)metre^(-3)
    _species_var[18] = PFILE(32) * 999.9999999999999;
    //V_C.durv, mw028533fe_f91e_499b_9ce0_ae17c7567d57, index: 19
    //Unit: mole^(1)metre^(-3)
    _species_var[19] = PFILE(33) * 999.9999999999999;
    //V_C.ipi, mw94268282_5fa1_42d0_9e89_331093861419, index: 20
    //Unit: mole^(1)metre^(-3)
    _species_var[20] = PFILE(34) * 999.9999999999999;
    //V_C.ENT, mw5b552ea7_7f80_427f_a541_6bbcc3fdd140, index: 21
    //Unit: mole^(1)metre^(-3)
    _species_var[21] = PFILE(35) * 999.9999999999999;
    //V_C.ENT_Buccal, mwe7dd5c6e_cf32_45e7_9245_847793d4c6df, index: 22
    //Unit: mole^(1)metre^(-3)
    _species_var[22] = PFILE(36) * 999.9999999999999;
    //V_C.ENT_GI, mw252a9877_9789_415c_84d5_3e41ee75f280, index: 23
    //Unit: mole^(1)metre^(-3)
    _species_var[23] = PFILE(37) * 999.9999999999999;
    //V_C.Dose2, mw7571b801_fbb9_4058_bffe_434c1b5d573d, index: 24
    //Unit: mole^(1)metre^(-3)
    _species_var[24] = PFILE(38) * 999.9999999999999;
    //V_P.T0, mw0d02324e_6766_4a9b_ada5_1fe0a1fd7935, index: 25
    //Unit: mole^(1)
    _species_var[25] = PFILE(39) * 1.66053872801495e-24;
    //V_P.T1, mwfd7bc971_e62b_4005_812d_863661b5197c, index: 26
    //Unit: mole^(1)
    _species_var[26] = PFILE(40) * 1.66053872801495e-24;
    //V_P.T2, mw6b15ff60_351c_45fe_a766_fb9f8174b4aa, index: 27
    //Unit: mole^(1)
    _species_var[27] = PFILE(41) * 1.66053872801495e-24;
    //V_P.T3, mwb12f8a41_3d54_4e7c_84f7_8cf19d702f8b, index: 28
    //Unit: mole^(1)
    _species_var[28] = PFILE(42) * 1.66053872801495e-24;
    //V_P.T4, mw853e4c5c_693c_42de_a86d_738188c350d1, index: 29
    //Unit: mole^(1)
    _species_var[29] = PFILE(43) * 1.66053872801495e-24;
    //V_P.T5, mwf6e3139c_a3af_4731_ac63_8828dc0b70c1, index: 30
    //Unit: mole^(1)
    _species_var[30] = PFILE(44) * 1.66053872801495e-24;
    //V_P.T6, mwf81bdb37_89a5_425d_9004_907db14585ff, index: 31
    //Unit: mole^(1)
    _species_var[31] = PFILE(45) * 1.66053872801495e-24;
    //V_P.T7, mwaa20148a_c4d6_4a36_aab3_76e704148dc2, index: 32
    //Unit: mole^(1)
    _species_var[32] = PFILE(46) * 1.66053872801495e-24;
    //V_P.T8, mw4a8b497c_4f9a_4ce3_94ad_2825c1136163, index: 33
    //Unit: mole^(1)
    _species_var[33] = PFILE(47) * 1.66053872801495e-24;
    //V_P.T9, mw131f4b3a_421c_4850_ba3c_2e1960bda421, index: 34
    //Unit: mole^(1)
    _species_var[34] = PFILE(48) * 1.66053872801495e-24;
    //V_P.T10, mw886a9450_7bb8_4d70_b235_29e97b20775a, index: 35
    //Unit: mole^(1)
    _species_var[35] = PFILE(49) * 1.66053872801495e-24;
    //V_P.T11, mwe94f3b7c_0d65_4f7e_a6a9_646709cfd788, index: 36
    //Unit: mole^(1)
    _species_var[36] = PFILE(50) * 1.66053872801495e-24;
    //V_P.T12, mwb55944ec_c742_4ad9_9321_041c85d05bc3, index: 37
    //Unit: mole^(1)
    _species_var[37] = PFILE(51) * 1.66053872801495e-24;
    //V_P.T13, mwa9493a01_bc9b_440d_b71d_20687bcc9044, index: 38
    //Unit: mole^(1)
    _species_var[38] = PFILE(52) * 1.66053872801495e-24;
    //V_P.T14, mwc4737530_b26b_4035_8c8d_97b59bd47417, index: 39
    //Unit: mole^(1)
    _species_var[39] = PFILE(53) * 1.66053872801495e-24;
    //V_P.T15, mw27bb1cad_7782_46a3_be09_a5a7be812d18, index: 40
    //Unit: mole^(1)
    _species_var[40] = PFILE(54) * 1.66053872801495e-24;
    //V_P.T16, mwb4be94e9_b04d_41ae_8164_5dbbfabb04f9, index: 41
    //Unit: mole^(1)
    _species_var[41] = PFILE(55) * 1.66053872801495e-24;
    //V_P.T17, mw89636a3c_b7fc_429b_80ca_435a399c7b74, index: 42
    //Unit: mole^(1)
    _species_var[42] = PFILE(56) * 1.66053872801495e-24;
    //V_P.nivo, mw3f8bdae7_8238_4539_a7a5_5eb96f79ef5e, index: 43
    //Unit: mole^(1)metre^(-3)
    _species_var[43] = PFILE(57) * 999.9999999999999;
    //V_P.durv, mwd757dd80_9c20_46dc_95cf_80c2b7be1f42, index: 44
    //Unit: mole^(1)metre^(-3)
    _species_var[44] = PFILE(58) * 999.9999999999999;
    //V_P.ipi, mw4929199a_deff_4344_9065_064c38a55ed2, index: 45
    //Unit: mole^(1)metre^(-3)
    _species_var[45] = PFILE(59) * 999.9999999999999;
    //V_P.Treg_CTLA4, mw3221af27_b92f_4e3c_bd21_1f6f1d2258a9, index: 46
    //Unit: mole^(1)
    _species_var[46] = PFILE(60) * 1.66053872801495e-24;
    //V_P.Treg_CTLA4_ipi, mwa1528c57_1975_4924_9bd3_149ccd6f0c91, index: 47
    //Unit: mole^(1)
    _species_var[47] = PFILE(61) * 1.66053872801495e-24;
    //V_P.Treg_CTLA4_ipi_CTLA4, mw88f8ddae_19a1_4fff_9a25_a8639e40e10a, index: 48
    //Unit: mole^(1)
    _species_var[48] = PFILE(62) * 1.66053872801495e-24;
    //V_P.ENT, mw45079498_4103_498e_876a_5c694739a519, index: 49
    //Unit: mole^(1)metre^(-3)
    _species_var[49] = PFILE(63) * 999.9999999999999;
    //V_T.C_x, mw075cd2c7_3213_44ca_a039_1568ab74cfea, index: 50
    //Unit: mole^(1)
    _species_var[50] = PFILE(64) * 1.66053872801495e-24;
    //V_T.T_exh, mwf7d843c3_7b5a_4ee2_94af_4fef2aacf5f5, index: 51
    //Unit: mole^(1)
    _species_var[51] = PFILE(65) * 1.66053872801495e-24;
    //V_T.C1, mw166e0979_99b2_41df_895a_ea3c00a45838, index: 52
    //Unit: mole^(1)
    _species_var[52] = PFILE(66) * 1.66053872801495e-24;
    //V_T.T0, mw37ae46fd_f396_4927_9ed8_b3baf4f98151, index: 53
    //Unit: mole^(1)
    _species_var[53] = PFILE(67) * 1.66053872801495e-24;
    //V_T.T1, mw6f683497_36e9_4229_ae99_3efdbc50f159, index: 54
    //Unit: mole^(1)
    _species_var[54] = PFILE(68) * 1.66053872801495e-24;
    //V_T.T2, mw453cbed5_335c_4afe_b8ce_e3bfba7f95af, index: 55
    //Unit: mole^(1)
    _species_var[55] = PFILE(69) * 1.66053872801495e-24;
    //V_T.T3, mwc25d2f4f_0887_49c4_9974_54a9fffc6be4, index: 56
    //Unit: mole^(1)
    _species_var[56] = PFILE(70) * 1.66053872801495e-24;
    //V_T.T4, mw0fb1b0f2_650b_4217_a01f_7043c38d5e5d, index: 57
    //Unit: mole^(1)
    _species_var[57] = PFILE(71) * 1.66053872801495e-24;
    //V_T.T5, mwc374815a_7c05_4fb6_9112_d9c9eaa81458, index: 58
    //Unit: mole^(1)
    _species_var[58] = PFILE(72) * 1.66053872801495e-24;
    //V_T.T6, mw89f3907e_c8fe_4d38_bbc6_ac2c8f30e41e, index: 59
    //Unit: mole^(1)
    _species_var[59] = PFILE(73) * 1.66053872801495e-24;
    //V_T.T7, mwb4c32cbb_688c_43c2_ba8c_3d71cdc3704a, index: 60
    //Unit: mole^(1)
    _species_var[60] = PFILE(74) * 1.66053872801495e-24;
    //V_T.T8, mw73ef4e84_fc33_4fc1_8806_7f6b62fe028d, index: 61
    //Unit: mole^(1)
    _species_var[61] = PFILE(75) * 1.66053872801495e-24;
    //V_T.T9, mw2bc6166d_243f_4727_8ae9_ad2e8a4bd855, index: 62
    //Unit: mole^(1)
    _species_var[62] = PFILE(76) * 1.66053872801495e-24;
    //V_T.T10, mw0123a5c0_65c1_4fe9_8f33_72fbed24e02f, index: 63
    //Unit: mole^(1)
    _species_var[63] = PFILE(77) * 1.66053872801495e-24;
    //V_T.T11, mwfc3724f9_fbb6_41b6_9561_a05e9aa4dc38, index: 64
    //Unit: mole^(1)
    _species_var[64] = PFILE(78) * 1.66053872801495e-24;
    //V_T.T12, mw282f3315_a779_48d1_a29f_9034f7ae49cf, index: 65
    //Unit: mole^(1)
    _species_var[65] = PFILE(79) * 1.66053872801495e-24;
    //V_T.T13, mwa6915107_c043_48b6_b88c_84be459564ae, index: 66
    //Unit: mole^(1)
    _species_var[66] = PFILE(80) * 1.66053872801495e-24;
    //V_T.T14, mw28528359_591e_495f_b22e_a90369f5bc14, index: 67
    //Unit: mole^(1)
    _species_var[67] = PFILE(81) * 1.66053872801495e-24;
    //V_T.T15, mwf74d972f_de8b_459b_b687_87e83582a23c, index: 68
    //Unit: mole^(1)
    _species_var[68] = PFILE(82) * 1.66053872801495e-24;
    //V_T.T16, mw54792167_c6de_4a77_b01a_29254b1c790a, index: 69
    //Unit: mole^(1)
    _species_var[69] = PFILE(83) * 1.66053872801495e-24;
    //V_T.T17, mwf65e0441_0453_4898_81df_a3052c62f0da, index: 70
    //Unit: mole^(1)
    _species_var[70] = PFILE(84) * 1.66053872801495e-24;
    //V_T.APC, mw4423bdbf_a379_45f4_ae29_52201258722e, index: 71
    //Unit: mole^(1)
    _species_var[71] = PFILE(85) * 1.66053872801495e-24;
    //V_T.mAPC, mwdb34e6e8_16be_4062_b94d_309cf631dbc3, index: 72
    //Unit: mole^(1)
    _species_var[72] = PFILE(86) * 1.66053872801495e-24;
    //V_T.c, mw39588528_a3d0_4dbf_bd5e_46e7a7154d01, index: 73
    //Unit: mole^(1)metre^(-3)
    _species_var[73] = PFILE(87) * 1000000.0;
    //V_T.nivo, mw9ceb78e9_ecef_48ef_a199_50ebe79d1258, index: 74
    //Unit: mole^(1)metre^(-3)
    _species_var[74] = PFILE(88) * 1000000.0;
    //V_T.durv, mw9409587e_c189_4940_9d15_76653a62c851, index: 75
    //Unit: mole^(1)metre^(-3)
    _species_var[75] = PFILE(89) * 1000000.0;
    //V_T.ipi, mw2b266ab9_41ce_4c0c_9cb7_179052988e9a, index: 76
    //Unit: mole^(1)metre^(-3)
    _species_var[76] = PFILE(90) * 1000000.0;
    //V_T.Treg_CTLA4, mw71d524c2_e2fe_4148_a010_2bbea05850f9, index: 77
    //Unit: mole^(1)
    _species_var[77] = PFILE(91) * 1.66053872801495e-24;
    //V_T.Treg_CTLA4_ipi, mw27f0c60d_791d_4ca7_8f72_54fd6178ddb9, index: 78
    //Unit: mole^(1)
    _species_var[78] = PFILE(92) * 1.66053872801495e-24;
    //V_T.Treg_CTLA4_ipi_CTLA4, mw41ac9e9c_959c_4d45_bd66_9cb31d776f8c, index: 79
    //Unit: mole^(1)
    _species_var[79] = PFILE(93) * 1.66053872801495e-24;
    //V_T.MDSC, mwe7a12d73_7de9_4eea_bc39_c6ef99bc1fc1, index: 80
    //Unit: mole^(1)
    _species_var[80] = PFILE(94) * 1.66053872801495e-24;
    //V_T.CCL2, mwba56dc93_d4c2_48cf_b475_d1be1a6ba0a5, index: 81
    //Unit: mole^(1)metre^(-3)
    _species_var[81] = PFILE(95) * 1000000.0;
    //V_T.NO, mwf1e69623_307f_4a65_84d3_e1f080541e58, index: 82
    //Unit: mole^(1)metre^(-3)
    _species_var[82] = PFILE(96) * 1000000.0;
    //V_T.ArgI, mwe34fffe5_4646_44be_8e97_bf12c8539ea0, index: 83
    //Unit: mole^(1)metre^(-3)
    _species_var[83] = PFILE(97) * 1000000.0;
    //V_T.ENT, mwd8d8ca9c_7e30_4537_81be_1b79064bdd01, index: 84
    //Unit: mole^(1)metre^(-3)
    _species_var[84] = PFILE(98) * 1000000.0;
    //V_LN.nT0, mw920e67b9_62f6_48b3_8a9c_6544547bb235, index: 85
    //Unit: mole^(1)
    _species_var[85] = PFILE(99) * 1.66053872801495e-24;
    //V_LN.aT0, mw89a9ef97_3923_45b5_85c6_4cd4be9feefd, index: 86
    //Unit: mole^(1)
    _species_var[86] = PFILE(100) * 1.66053872801495e-24;
    //V_LN.T0, mwd2830144_a864_49fb_ab28_ec735c74c11a, index: 87
    //Unit: mole^(1)
    _species_var[87] = PFILE(101) * 1.66053872801495e-24;
    //V_LN.IL2, mwbf918460_0454_43de_a912_30843acc673c, index: 88
    //Unit: mole^(1)metre^(-3)
    _species_var[88] = PFILE(102) * 999999999.9999999;
    //V_LN.nT1, mw0ab82d1e_c61b_4761_8e9a_01362e11da39, index: 89
    //Unit: mole^(1)
    _species_var[89] = PFILE(103) * 1.66053872801495e-24;
    //V_LN.aT1, mwd8d25183_7dac_4707_918a_550db06f821d, index: 90
    //Unit: mole^(1)
    _species_var[90] = PFILE(104) * 1.66053872801495e-24;
    //V_LN.T1, mw69bb9ecc_1a53_4f24_9997_328b8bb33bef, index: 91
    //Unit: mole^(1)
    _species_var[91] = PFILE(105) * 1.66053872801495e-24;
    //V_LN.nT2, mw4881c045_03fb_4763_a6d3_d628d353ad25, index: 92
    //Unit: mole^(1)
    _species_var[92] = PFILE(106) * 1.66053872801495e-24;
    //V_LN.aT2, mw9dbf720a_83fa_4830_b7cf_8708aebc36e9, index: 93
    //Unit: mole^(1)
    _species_var[93] = PFILE(107) * 1.66053872801495e-24;
    //V_LN.T2, mw992e0289_f0db_436b_82f4_4c81bcb5c420, index: 94
    //Unit: mole^(1)
    _species_var[94] = PFILE(108) * 1.66053872801495e-24;
    //V_LN.nT3, mw51402d86_bd03_48e7_b097_4650354cc845, index: 95
    //Unit: mole^(1)
    _species_var[95] = PFILE(109) * 1.66053872801495e-24;
    //V_LN.aT3, mw0ec88a8d_394d_4dfe_8378_d96cd597e96d, index: 96
    //Unit: mole^(1)
    _species_var[96] = PFILE(110) * 1.66053872801495e-24;
    //V_LN.T3, mwa1f944cc_d860_451f_89f2_46176538c30d, index: 97
    //Unit: mole^(1)
    _species_var[97] = PFILE(111) * 1.66053872801495e-24;
    //V_LN.nT4, mw2d422ced_6bbd_4f0d_91cf_9eace75d58dc, index: 98
    //Unit: mole^(1)
    _species_var[98] = PFILE(112) * 1.66053872801495e-24;
    //V_LN.aT4, mwb77411c8_bfb8_42d6_854a_7ac6f4c13dc2, index: 99
    //Unit: mole^(1)
    _species_var[99] = PFILE(113) * 1.66053872801495e-24;
    //V_LN.T4, mwd11181a0_579e_45a3_ae3f_c00a7cf991d1, index: 100
    //Unit: mole^(1)
    _species_var[100] = PFILE(114) * 1.66053872801495e-24;
    //V_LN.nT5, mwb641a2a2_89c2_440d_ba14_fe94ca207b5c, index: 101
    //Unit: mole^(1)
    _species_var[101] = PFILE(115) * 1.66053872801495e-24;
    //V_LN.aT5, mw89a4059a_4b54_41de_a651_8def0ee81fe4, index: 102
    //Unit: mole^(1)
    _species_var[102] = PFILE(116) * 1.66053872801495e-24;
    //V_LN.T5, mwf079d340_fada_447e_bbbb_9f7801c27884, index: 103
    //Unit: mole^(1)
    _species_var[103] = PFILE(117) * 1.66053872801495e-24;
    //V_LN.nT6, mwf8bc1b17_9698_4a69_a1d5_caba680db09a, index: 104
    //Unit: mole^(1)
    _species_var[104] = PFILE(118) * 1.66053872801495e-24;
    //V_LN.aT6, mw143fe149_a878_439e_bd2c_cc4917a96651, index: 105
    //Unit: mole^(1)
    _species_var[105] = PFILE(119) * 1.66053872801495e-24;
    //V_LN.T6, mwde91eaa1_a0c7_4e62_8331_b419fead899d, index: 106
    //Unit: mole^(1)
    _species_var[106] = PFILE(120) * 1.66053872801495e-24;
    //V_LN.nT7, mw18ba6320_5fff_4d53_8bc4_7c0a87aa1022, index: 107
    //Unit: mole^(1)
    _species_var[107] = PFILE(121) * 1.66053872801495e-24;
    //V_LN.aT7, mw25588f00_3c01_4091_9817_b001f15216b7, index: 108
    //Unit: mole^(1)
    _species_var[108] = PFILE(122) * 1.66053872801495e-24;
    //V_LN.T7, mw2f5cd66f_5d5b_4a98_a2ca_a13a5635cb4f, index: 109
    //Unit: mole^(1)
    _species_var[109] = PFILE(123) * 1.66053872801495e-24;
    //V_LN.nT8, mw83577716_314a_4a7e_b4e0_4068eafc0f0f, index: 110
    //Unit: mole^(1)
    _species_var[110] = PFILE(124) * 1.66053872801495e-24;
    //V_LN.aT8, mw48a54685_49b2_4adc_a114_12ab366dc608, index: 111
    //Unit: mole^(1)
    _species_var[111] = PFILE(125) * 1.66053872801495e-24;
    //V_LN.T8, mwe2415429_d3b3_48b1_bcc6_661d15fb5c66, index: 112
    //Unit: mole^(1)
    _species_var[112] = PFILE(126) * 1.66053872801495e-24;
    //V_LN.nT9, mwb876fb85_e75a_430c_9bf9_1229186b6ced, index: 113
    //Unit: mole^(1)
    _species_var[113] = PFILE(127) * 1.66053872801495e-24;
    //V_LN.aT9, mw4164da68_0182_484f_994f_9593af1de84b, index: 114
    //Unit: mole^(1)
    _species_var[114] = PFILE(128) * 1.66053872801495e-24;
    //V_LN.T9, mw0051777d_f601_4a17_84d6_668d2141870c, index: 115
    //Unit: mole^(1)
    _species_var[115] = PFILE(129) * 1.66053872801495e-24;
    //V_LN.nT10, mwbb0a6ad0_c9ee_43cf_8557_5f061762a7ee, index: 116
    //Unit: mole^(1)
    _species_var[116] = PFILE(130) * 1.66053872801495e-24;
    //V_LN.aT10, mw107ff18e_9183_45b0_b003_0006ee0a5ee6, index: 117
    //Unit: mole^(1)
    _species_var[117] = PFILE(131) * 1.66053872801495e-24;
    //V_LN.T10, mwcf606105_63cc_41e4_939c_d42d06e2c6b8, index: 118
    //Unit: mole^(1)
    _species_var[118] = PFILE(132) * 1.66053872801495e-24;
    //V_LN.nT11, mw218bdaff_63e2_4ed7_b981_96bf29e1a32d, index: 119
    //Unit: mole^(1)
    _species_var[119] = PFILE(133) * 1.66053872801495e-24;
    //V_LN.aT11, mwb59891a3_6ce4_46b8_af32_85c3c9d87771, index: 120
    //Unit: mole^(1)
    _species_var[120] = PFILE(134) * 1.66053872801495e-24;
    //V_LN.T11, mw3b16f0c0_0bdf_4bfd_a2fd_b62eac257d89, index: 121
    //Unit: mole^(1)
    _species_var[121] = PFILE(135) * 1.66053872801495e-24;
    //V_LN.nT12, mwf6a0a1e8_93c3_40b0_bb81_c31ebd0b0979, index: 122
    //Unit: mole^(1)
    _species_var[122] = PFILE(136) * 1.66053872801495e-24;
    //V_LN.aT12, mw3f8b0429_b9c8_4f9b_b5cd_a52fc6b3f705, index: 123
    //Unit: mole^(1)
    _species_var[123] = PFILE(137) * 1.66053872801495e-24;
    //V_LN.T12, mw11f6f713_0766_415a_874f_96ae4be29790, index: 124
    //Unit: mole^(1)
    _species_var[124] = PFILE(138) * 1.66053872801495e-24;
    //V_LN.nT13, mw7c8eb662_5ae0_4345_9348_36d868b141e3, index: 125
    //Unit: mole^(1)
    _species_var[125] = PFILE(139) * 1.66053872801495e-24;
    //V_LN.aT13, mw99f44159_9cdf_4d3e_aa00_5b7c3d25f462, index: 126
    //Unit: mole^(1)
    _species_var[126] = PFILE(140) * 1.66053872801495e-24;
    //V_LN.T13, mwc7b3db6c_fd11_4079_b9b9_b47d4165d023, index: 127
    //Unit: mole^(1)
    _species_var[127] = PFILE(141) * 1.66053872801495e-24;
    //V_LN.nT14, mwe4f646d8_1c3d_42ce_a9a8_80c82f2128f3, index: 128
    //Unit: mole^(1)
    _species_var[128] = PFILE(142) * 1.66053872801495e-24;
    //V_LN.aT14, mwf418ccd9_65bc_4219_895d_825f79e8ff48, index: 129
    //Unit: mole^(1)
    _species_var[129] = PFILE(143) * 1.66053872801495e-24;
    //V_LN.T14, mwa61d1489_c17f_47ac_945c_bfefad9e470b, index: 130
    //Unit: mole^(1)
    _species_var[130] = PFILE(144) * 1.66053872801495e-24;
    //V_LN.nT15, mw1a9ea930_30c0_4539_83e9_29dbeddf8e8e, index: 131
    //Unit: mole^(1)
    _species_var[131] = PFILE(145) * 1.66053872801495e-24;
    //V_LN.aT15, mw79053dd3_8b11_47d5_83ef_c42244f23284, index: 132
    //Unit: mole^(1)
    _species_var[132] = PFILE(146) * 1.66053872801495e-24;
    //V_LN.T15, mw5ff766e9_af79_4701_921f_e413cdcb7a92, index: 133
    //Unit: mole^(1)
    _species_var[133] = PFILE(147) * 1.66053872801495e-24;
    //V_LN.nT16, mw899dbc66_7bf5_4c6a_a832_23f4b222c1dc, index: 134
    //Unit: mole^(1)
    _species_var[134] = PFILE(148) * 1.66053872801495e-24;
    //V_LN.aT16, mw1152f21f_3b27_41ab_a58a_eb82b6985a44, index: 135
    //Unit: mole^(1)
    _species_var[135] = PFILE(149) * 1.66053872801495e-24;
    //V_LN.T16, mwc5695243_2e2d_457f_9ab7_72af8ada7f65, index: 136
    //Unit: mole^(1)
    _species_var[136] = PFILE(150) * 1.66053872801495e-24;
    //V_LN.nT17, mwbe0683ce_bbbc_4626_8989_d3831f332e71, index: 137
    //Unit: mole^(1)
    _species_var[137] = PFILE(151) * 1.66053872801495e-24;
    //V_LN.aT17, mw50683f94_425e_437c_991c_28bebf0f1f86, index: 138
    //Unit: mole^(1)
    _species_var[138] = PFILE(152) * 1.66053872801495e-24;
    //V_LN.T17, mw63fba837_1bc4_4dba_8c19_d47878e1f494, index: 139
    //Unit: mole^(1)
    _species_var[139] = PFILE(153) * 1.66053872801495e-24;
    //V_LN.APC, mw663f380f_3f97_4e09_bc80_16568ae946c9, index: 140
    //Unit: mole^(1)
    _species_var[140] = PFILE(154) * 1.66053872801495e-24;
    //V_LN.mAPC, mwb8bd1f0a_e1d7_4c91_a45c_377d1c692e7d, index: 141
    //Unit: mole^(1)
    _species_var[141] = PFILE(155) * 1.66053872801495e-24;
    //V_LN.P0, mwb5fecf35_53d3_4ab7_8918_6a417c9d6ab8, index: 142
    //Unit: mole^(1)metre^(-3)
    _species_var[142] = PFILE(156) * 999999999.9999999;
    //V_LN.P1, mw553f3a6f_d519_4103_8a71_cdb117aa394c, index: 143
    //Unit: mole^(1)metre^(-3)
    _species_var[143] = PFILE(157) * 999999999.9999999;
    //V_LN.P2, mw34a53628_9cf6_40a2_a46e_d5d02c0070ef, index: 144
    //Unit: mole^(1)metre^(-3)
    _species_var[144] = PFILE(158) * 999999999.9999999;
    //V_LN.P3, mw28230303_2f2f_402a_a7e7_f783a3ec7a47, index: 145
    //Unit: mole^(1)metre^(-3)
    _species_var[145] = PFILE(159) * 999999999.9999999;
    //V_LN.P4, mw63e3f4eb_2dda_43de_9d78_cf2e7df57b24, index: 146
    //Unit: mole^(1)metre^(-3)
    _species_var[146] = PFILE(160) * 999999999.9999999;
    //V_LN.P5, mw403f49d1_f960_4a88_a4d8_051c55c20aed, index: 147
    //Unit: mole^(1)metre^(-3)
    _species_var[147] = PFILE(161) * 999999999.9999999;
    //V_LN.P6, mwb670e6fc_0b3b_43d1_9845_78ecc071444b, index: 148
    //Unit: mole^(1)metre^(-3)
    _species_var[148] = PFILE(162) * 999999999.9999999;
    //V_LN.P7, mw368d8bee_5543_446b_9c5f_b6a1fb5a73be, index: 149
    //Unit: mole^(1)metre^(-3)
    _species_var[149] = PFILE(163) * 999999999.9999999;
    //V_LN.P8, mw1980ae68_0e95_4d11_8667_efb87a57c965, index: 150
    //Unit: mole^(1)metre^(-3)
    _species_var[150] = PFILE(164) * 999999999.9999999;
    //V_LN.P9, mw957d96d8_43ef_4e7b_a29c_1ace2ac3c8cf, index: 151
    //Unit: mole^(1)metre^(-3)
    _species_var[151] = PFILE(165) * 999999999.9999999;
    //V_LN.P10, mwe2ad5bed_ef1e_429b_a8d3_cff190650c3a, index: 152
    //Unit: mole^(1)metre^(-3)
    _species_var[152] = PFILE(166) * 999999999.9999999;
    //V_LN.P11, mw6a3d817d_d144_465c_8558_42b89f46d48f, index: 153
    //Unit: mole^(1)metre^(-3)
    _species_var[153] = PFILE(167) * 999999999.9999999;
    //V_LN.P12, mw0ecd1667_23cb_4254_ac30_a568a1567b99, index: 154
    //Unit: mole^(1)metre^(-3)
    _species_var[154] = PFILE(168) * 999999999.9999999;
    //V_LN.P13, mw867dbd51_a7ef_4aef_b330_2188d1fe78af, index: 155
    //Unit: mole^(1)metre^(-3)
    _species_var[155] = PFILE(169) * 999999999.9999999;
    //V_LN.P14, mw41dee060_98fe_4185_af26_b07f353361c7, index: 156
    //Unit: mole^(1)metre^(-3)
    _species_var[156] = PFILE(170) * 999999999.9999999;
    //V_LN.P15, mw54ce4b65_f2d5_4fc3_abcd_3a3ca6c14175, index: 157
    //Unit: mole^(1)metre^(-3)
    _species_var[157] = PFILE(171) * 999999999.9999999;
    //V_LN.P16, mwe5e01895_9e05_4ecc_a065_1a3de4f22ae3, index: 158
    //Unit: mole^(1)metre^(-3)
    _species_var[158] = PFILE(172) * 999999999.9999999;
    //V_LN.P17, mwe0a500db_2d88_4227_8638_d3987b7aeac8, index: 159
    //Unit: mole^(1)metre^(-3)
    _species_var[159] = PFILE(173) * 999999999.9999999;
    //V_LN.nivo, mw1d8daabb_f0b6_4d20_8334_c2d0b2191b2e, index: 160
    //Unit: mole^(1)metre^(-3)
    _species_var[160] = PFILE(174) * 999999999.9999999;
    //V_LN.durv, mw4efaa957_1491_4a52_a399_869b37f36e0e, index: 161
    //Unit: mole^(1)metre^(-3)
    _species_var[161] = PFILE(175) * 999999999.9999999;
    //V_LN.ipi, mw48441613_46e2_4485_abb0_564c42e78358, index: 162
    //Unit: mole^(1)metre^(-3)
    _species_var[162] = PFILE(176) * 999999999.9999999;
    //V_LN.ENT, mw514d8d00_53af_4738_9b64_aa90a62bf7ff, index: 163
    //Unit: mole^(1)metre^(-3)
    _species_var[163] = PFILE(177) * 999999999.9999999;
    //V_e.P0, mw030d747e_ba4c_4d54_bb88_6f75c66313ee, index: 164
    //Unit: mole^(1)metre^(-3)
    _species_var[164] = PFILE(178) * 999.9999999999999;
    //V_e.p0, mw34726bc3_9491_44d7_a973_d881fcbeb07d, index: 165
    //Unit: mole^(1)metre^(-3)
    _species_var[165] = PFILE(179) * 999.9999999999999;
    //V_e.P1, mw9e468cda_3810_400f_9e18_e48f6548bc8b, index: 166
    //Unit: mole^(1)metre^(-3)
    _species_var[166] = PFILE(180) * 999.9999999999999;
    //V_e.p1, mw32dfa816_e409_4edb_8e98_21c9f627c266, index: 167
    //Unit: mole^(1)metre^(-3)
    _species_var[167] = PFILE(181) * 999.9999999999999;
    //V_e.P2, mwacdfa785_39e4_432b_873e_b08df7e1a203, index: 168
    //Unit: mole^(1)metre^(-3)
    _species_var[168] = PFILE(182) * 999.9999999999999;
    //V_e.p2, mw258ee6d3_0cba_4708_ab32_83a5a8b4c1fd, index: 169
    //Unit: mole^(1)metre^(-3)
    _species_var[169] = PFILE(183) * 999.9999999999999;
    //V_e.P3, mw92476e2a_94d5_4d9a_8143_8ce235f25167, index: 170
    //Unit: mole^(1)metre^(-3)
    _species_var[170] = PFILE(184) * 999.9999999999999;
    //V_e.p3, mw85e5836b_e5c8_4255_bb8b_d17dce34310c, index: 171
    //Unit: mole^(1)metre^(-3)
    _species_var[171] = PFILE(185) * 999.9999999999999;
    //V_e.P4, mw939716e1_c9f1_4d46_9650_95ee33bd5d9f, index: 172
    //Unit: mole^(1)metre^(-3)
    _species_var[172] = PFILE(186) * 999.9999999999999;
    //V_e.p4, mw6992b072_8d47_45ca_8e63_22a75a8d8b04, index: 173
    //Unit: mole^(1)metre^(-3)
    _species_var[173] = PFILE(187) * 999.9999999999999;
    //V_e.P5, mwa081236b_8fd3_4dc1_8bec_fe043fc17e54, index: 174
    //Unit: mole^(1)metre^(-3)
    _species_var[174] = PFILE(188) * 999.9999999999999;
    //V_e.p5, mw2cbf8555_f5df_4221_bcdc_f18b25c4ac8f, index: 175
    //Unit: mole^(1)metre^(-3)
    _species_var[175] = PFILE(189) * 999.9999999999999;
    //V_e.P6, mwc9e768d3_db42_421e_8005_13e0d26b11a0, index: 176
    //Unit: mole^(1)metre^(-3)
    _species_var[176] = PFILE(190) * 999.9999999999999;
    //V_e.p6, mw0b05cf8f_3613_4af1_a769_95267fdfce7f, index: 177
    //Unit: mole^(1)metre^(-3)
    _species_var[177] = PFILE(191) * 999.9999999999999;
    //V_e.P7, mw591b9c0b_a429_4276_8f6f_850edbc6d6b8, index: 178
    //Unit: mole^(1)metre^(-3)
    _species_var[178] = PFILE(192) * 999.9999999999999;
    //V_e.p7, mw5e962ae5_19b1_4170_969b_393c00e5cf75, index: 179
    //Unit: mole^(1)metre^(-3)
    _species_var[179] = PFILE(193) * 999.9999999999999;
    //V_e.P8, mwc11d6511_0dc4_4fda_83e9_cd6216ba3639, index: 180
    //Unit: mole^(1)metre^(-3)
    _species_var[180] = PFILE(194) * 999.9999999999999;
    //V_e.p8, mw8653b9c2_168f_46a9_ae92_08511c1cf9c7, index: 181
    //Unit: mole^(1)metre^(-3)
    _species_var[181] = PFILE(195) * 999.9999999999999;
    //V_e.P9, mwb15d3c82_bc62_424c_98c1_7000bbf8b8b1, index: 182
    //Unit: mole^(1)metre^(-3)
    _species_var[182] = PFILE(196) * 999.9999999999999;
    //V_e.p9, mw215c00d2_7594_4d40_b507_e45094042aa5, index: 183
    //Unit: mole^(1)metre^(-3)
    _species_var[183] = PFILE(197) * 999.9999999999999;
    //V_e.P10, mw88413228_e4e0_4400_8150_abfc9ba6a26a, index: 184
    //Unit: mole^(1)metre^(-3)
    _species_var[184] = PFILE(198) * 999.9999999999999;
    //V_e.p10, mw498498a3_a8e8_4159_81b0_c9ac58441cdc, index: 185
    //Unit: mole^(1)metre^(-3)
    _species_var[185] = PFILE(199) * 999.9999999999999;
    //V_e.P11, mw2ac123a3_4511_4a31_9d20_bad9aac3c7d8, index: 186
    //Unit: mole^(1)metre^(-3)
    _species_var[186] = PFILE(200) * 999.9999999999999;
    //V_e.p11, mw47ec7d37_43e0_4a82_92c7_0e169d70ceab, index: 187
    //Unit: mole^(1)metre^(-3)
    _species_var[187] = PFILE(201) * 999.9999999999999;
    //V_e.P12, mwda99e494_9c81_497d_8424_cd3a3bc18102, index: 188
    //Unit: mole^(1)metre^(-3)
    _species_var[188] = PFILE(202) * 999.9999999999999;
    //V_e.p12, mw9f2b477f_9323_4f90_b454_61ea870fb275, index: 189
    //Unit: mole^(1)metre^(-3)
    _species_var[189] = PFILE(203) * 999.9999999999999;
    //V_e.P13, mw9776d384_703f_4464_9932_c5d9aaba0ff8, index: 190
    //Unit: mole^(1)metre^(-3)
    _species_var[190] = PFILE(204) * 999.9999999999999;
    //V_e.p13, mw194d17e9_3928_452b_bb16_bb3fe66a9621, index: 191
    //Unit: mole^(1)metre^(-3)
    _species_var[191] = PFILE(205) * 999.9999999999999;
    //V_e.P14, mw4099763b_0792_48c8_9dc5_a92fa73b3979, index: 192
    //Unit: mole^(1)metre^(-3)
    _species_var[192] = PFILE(206) * 999.9999999999999;
    //V_e.p14, mwb7eba224_7e64_487a_87bd_b73b219cd506, index: 193
    //Unit: mole^(1)metre^(-3)
    _species_var[193] = PFILE(207) * 999.9999999999999;
    //V_e.P15, mw17fc38d7_be9b_4786_ab6f_864cbd606f53, index: 194
    //Unit: mole^(1)metre^(-3)
    _species_var[194] = PFILE(208) * 999.9999999999999;
    //V_e.p15, mwb412fa67_e1c9_44a9_8b7b_5402ba78b6a9, index: 195
    //Unit: mole^(1)metre^(-3)
    _species_var[195] = PFILE(209) * 999.9999999999999;
    //V_e.P16, mw321cff9b_bb9a_421d_b06c_283474207107, index: 196
    //Unit: mole^(1)metre^(-3)
    _species_var[196] = PFILE(210) * 999.9999999999999;
    //V_e.p16, mw3f164f96_8fb1_412f_bdd6_a7a7cc5b45af, index: 197
    //Unit: mole^(1)metre^(-3)
    _species_var[197] = PFILE(211) * 999.9999999999999;
    //V_e.P17, mw1a60e594_2069_47d8_85d0_eb12b9e298e9, index: 198
    //Unit: mole^(1)metre^(-3)
    _species_var[198] = PFILE(212) * 999.9999999999999;
    //V_e.p17, mwcf4052d2_8902_42e0_b4d6_559a10de597b, index: 199
    //Unit: mole^(1)metre^(-3)
    _species_var[199] = PFILE(213) * 999.9999999999999;
    //A_e.M1, mwc38f9591_f266_441a_90f4_54778fcb19e4, index: 200
    //Unit: mole^(1)metre^(-2)
    _species_var[200] = PFILE(214) * 1.66053872801495e-12;
    //A_e.M1p0, mw7ef0d1d1_4819_4d82_addd_e034a414c19b, index: 201
    //Unit: mole^(1)metre^(-2)
    _species_var[201] = PFILE(215) * 1.66053872801495e-12;
    //A_e.M1p1, mw7a11e985_1284_4f96_93b6_b9cd2eee660a, index: 202
    //Unit: mole^(1)metre^(-2)
    _species_var[202] = PFILE(216) * 1.66053872801495e-12;
    //A_e.M1p2, mwad69c11f_9062_4a1f_9e1e_f833e982da13, index: 203
    //Unit: mole^(1)metre^(-2)
    _species_var[203] = PFILE(217) * 1.66053872801495e-12;
    //A_e.M1p3, mw1b35aced_867d_485e_a534_67a3c9e76842, index: 204
    //Unit: mole^(1)metre^(-2)
    _species_var[204] = PFILE(218) * 1.66053872801495e-12;
    //A_e.M1p4, mw5a5b188a_2e3f_4797_952f_2624de381ad0, index: 205
    //Unit: mole^(1)metre^(-2)
    _species_var[205] = PFILE(219) * 1.66053872801495e-12;
    //A_e.M1p5, mw31b46b0a_0d11_46c7_bed7_225d477fb3c4, index: 206
    //Unit: mole^(1)metre^(-2)
    _species_var[206] = PFILE(220) * 1.66053872801495e-12;
    //A_e.M1p6, mwcda87a39_c93f_4325_ada3_eaa7d943ead9, index: 207
    //Unit: mole^(1)metre^(-2)
    _species_var[207] = PFILE(221) * 1.66053872801495e-12;
    //A_e.M1p7, mwc9faa083_05ea_44ab_8eef_19a632929685, index: 208
    //Unit: mole^(1)metre^(-2)
    _species_var[208] = PFILE(222) * 1.66053872801495e-12;
    //A_e.M1p8, mw1bbfc86c_5124_44b8_bb1f_a1cfccd9d7f6, index: 209
    //Unit: mole^(1)metre^(-2)
    _species_var[209] = PFILE(223) * 1.66053872801495e-12;
    //A_e.M1p9, mwfd9689b7_eb2a_4fd1_a603_67a943a8a80b, index: 210
    //Unit: mole^(1)metre^(-2)
    _species_var[210] = PFILE(224) * 1.66053872801495e-12;
    //A_e.M1p10, mw4e5bdcd1_8641_48af_9d11_87c7ffd89d05, index: 211
    //Unit: mole^(1)metre^(-2)
    _species_var[211] = PFILE(225) * 1.66053872801495e-12;
    //A_e.M1p11, mw52cf88d0_f646_4f6f_a3a0_553144eafc5e, index: 212
    //Unit: mole^(1)metre^(-2)
    _species_var[212] = PFILE(226) * 1.66053872801495e-12;
    //A_e.M1p12, mwe0a877db_03d9_4b26_9d03_2b4f77fae16b, index: 213
    //Unit: mole^(1)metre^(-2)
    _species_var[213] = PFILE(227) * 1.66053872801495e-12;
    //A_e.M1p13, mwd4911096_ec90_45dd_aefa_0dc4f332b87c, index: 214
    //Unit: mole^(1)metre^(-2)
    _species_var[214] = PFILE(228) * 1.66053872801495e-12;
    //A_e.M1p14, mwcd103ae9_8e9d_4a91_9d79_56105f680016, index: 215
    //Unit: mole^(1)metre^(-2)
    _species_var[215] = PFILE(229) * 1.66053872801495e-12;
    //A_e.M1p15, mw42fc0a50_14f8_46ce_9a93_f075a8fead84, index: 216
    //Unit: mole^(1)metre^(-2)
    _species_var[216] = PFILE(230) * 1.66053872801495e-12;
    //A_e.M1p16, mw011515d9_7aa9_4107_b1c2_825a67583deb, index: 217
    //Unit: mole^(1)metre^(-2)
    _species_var[217] = PFILE(231) * 1.66053872801495e-12;
    //A_e.M1p17, mw3b2209e8_e7cf_4750_97f4_754ec135f003, index: 218
    //Unit: mole^(1)metre^(-2)
    _species_var[218] = PFILE(232) * 1.66053872801495e-12;
    //A_s.M1, mw0459e680_0784_4815_a106_0ee800c25c76, index: 219
    //Unit: mole^(1)metre^(-2)
    _species_var[219] = PFILE(233) * 1.66053872801495e-12;
    //A_s.M1p0, mwbedcc0de_536c_44d7_8bb5_59e09a9ae918, index: 220
    //Unit: mole^(1)metre^(-2)
    _species_var[220] = PFILE(234) * 1.66053872801495e-12;
    //A_s.M1p1, mwd9a17bbe_23a7_4890_9e91_f6913b8fb958, index: 221
    //Unit: mole^(1)metre^(-2)
    _species_var[221] = PFILE(235) * 1.66053872801495e-12;
    //A_s.M1p2, mw8db574d2_07b9_49a7_80bd_8a9b8c3007cd, index: 222
    //Unit: mole^(1)metre^(-2)
    _species_var[222] = PFILE(236) * 1.66053872801495e-12;
    //A_s.M1p3, mw254dba7a_5c2a_457a_a6a6_e99571cb917f, index: 223
    //Unit: mole^(1)metre^(-2)
    _species_var[223] = PFILE(237) * 1.66053872801495e-12;
    //A_s.M1p4, mwfb538002_d949_41a9_afd7_724719247696, index: 224
    //Unit: mole^(1)metre^(-2)
    _species_var[224] = PFILE(238) * 1.66053872801495e-12;
    //A_s.M1p5, mwfb301a02_5711_486b_a890_cf82759dfd1f, index: 225
    //Unit: mole^(1)metre^(-2)
    _species_var[225] = PFILE(239) * 1.66053872801495e-12;
    //A_s.M1p6, mw558a0b77_a222_4923_9f72_40d839b23b24, index: 226
    //Unit: mole^(1)metre^(-2)
    _species_var[226] = PFILE(240) * 1.66053872801495e-12;
    //A_s.M1p7, mw378fdc41_89af_47c2_8d4f_d695fd27b60a, index: 227
    //Unit: mole^(1)metre^(-2)
    _species_var[227] = PFILE(241) * 1.66053872801495e-12;
    //A_s.M1p8, mw4fcd88d9_3a0b_4a64_8ee3_ce40fe9716e1, index: 228
    //Unit: mole^(1)metre^(-2)
    _species_var[228] = PFILE(242) * 1.66053872801495e-12;
    //A_s.M1p9, mw5d5673f3_b5c8_4ce7_a2c2_b7ffabe9ab57, index: 229
    //Unit: mole^(1)metre^(-2)
    _species_var[229] = PFILE(243) * 1.66053872801495e-12;
    //A_s.M1p10, mw3378b345_ca89_4727_ad79_9df3e626f2f0, index: 230
    //Unit: mole^(1)metre^(-2)
    _species_var[230] = PFILE(244) * 1.66053872801495e-12;
    //A_s.M1p11, mwd9b127d1_0a1b_4aa4_af4f_6bc6bc42a40c, index: 231
    //Unit: mole^(1)metre^(-2)
    _species_var[231] = PFILE(245) * 1.66053872801495e-12;
    //A_s.M1p12, mw5859d1bf_e2df_405e_b359_613b7f8a2e56, index: 232
    //Unit: mole^(1)metre^(-2)
    _species_var[232] = PFILE(246) * 1.66053872801495e-12;
    //A_s.M1p13, mw1f392f4f_de53_4d7c_9fc0_f1a891424c85, index: 233
    //Unit: mole^(1)metre^(-2)
    _species_var[233] = PFILE(247) * 1.66053872801495e-12;
    //A_s.M1p14, mw29a131c7_cbbf_450b_9c19_c1bfac8a81dc, index: 234
    //Unit: mole^(1)metre^(-2)
    _species_var[234] = PFILE(248) * 1.66053872801495e-12;
    //A_s.M1p15, mw51f45325_47e9_43c3_81bc_5be4a71eb005, index: 235
    //Unit: mole^(1)metre^(-2)
    _species_var[235] = PFILE(249) * 1.66053872801495e-12;
    //A_s.M1p16, mw3e97bb76_ab1e_4213_be11_905007d958bf, index: 236
    //Unit: mole^(1)metre^(-2)
    _species_var[236] = PFILE(250) * 1.66053872801495e-12;
    //A_s.M1p17, mw87a5819c_cea3_4ab9_bbb3_2aa0f12d79d1, index: 237
    //Unit: mole^(1)metre^(-2)
    _species_var[237] = PFILE(251) * 1.66053872801495e-12;
    //syn_T_C1.PD1_PDL1, mw2db0fb9c_1cde_45f3_812c_211e3aaa12a7, index: 238
    //Unit: mole^(1)metre^(-2)
    _species_var[238] = PFILE(252) * 1.66053872801495e-12;
    //syn_T_C1.PD1_PDL2, mw38e21715_810f_4d22_a664_139304021014, index: 239
    //Unit: mole^(1)metre^(-2)
    _species_var[239] = PFILE(253) * 1.66053872801495e-12;
    //syn_T_C1.PD1, mw1bb203a4_6094_4b53_8e88_cb90f755b7b7, index: 240
    //Unit: mole^(1)metre^(-2)
    _species_var[240] = PFILE(254) * 1.66053872801495e-12;
    //syn_T_C1.PDL1, mw5aad9f11_4b01_41e1_bccb_c467f5816c91, index: 241
    //Unit: mole^(1)metre^(-2)
    _species_var[241] = PFILE(255) * 1.66053872801495e-12;
    //syn_T_C1.PDL2, mw5eb3e42d_364d_4ff3_8639_c3f341738226, index: 242
    //Unit: mole^(1)metre^(-2)
    _species_var[242] = PFILE(256) * 1.66053872801495e-12;
    //syn_T_C1.PD1_nivo, mw8bfd4c21_30de_4c8b_8212_1f927ad6c975, index: 243
    //Unit: mole^(1)metre^(-2)
    _species_var[243] = PFILE(257) * 1.66053872801495e-12;
    //syn_T_C1.PD1_nivo_PD1, mwb05e0618_492e_4c04_a739_eae47dd41cab, index: 244
    //Unit: mole^(1)metre^(-2)
    _species_var[244] = PFILE(258) * 1.66053872801495e-12;
    //syn_T_C1.PDL1_durv, mw2528a9a5_3237_405d_bd8a_d12808f82821, index: 245
    //Unit: mole^(1)metre^(-2)
    _species_var[245] = PFILE(259) * 1.66053872801495e-12;
    //syn_T_C1.PDL1_durv_PDL1, mw95adda89_7440_46de_8b62_28fd3911c372, index: 246
    //Unit: mole^(1)metre^(-2)
    _species_var[246] = PFILE(260) * 1.66053872801495e-12;
    //syn_T_C1.TPDL1, mw5144f0a3_b1f0_49ba_9d1b_41ad349725fd, index: 247
    //Unit: mole^(1)metre^(-2)
    _species_var[247] = PFILE(261) * 1.66053872801495e-12;
    //syn_T_C1.TPDL1_durv, mw3eeaaf14_9447_4a30_a69c_ac77aff34074, index: 248
    //Unit: mole^(1)metre^(-2)
    _species_var[248] = PFILE(262) * 1.66053872801495e-12;
    //syn_T_C1.TPDL1_durv_TPDL1, mwcfba9f77_0691_4a7a_a805_8b63c764bda9, index: 249
    //Unit: mole^(1)metre^(-2)
    _species_var[249] = PFILE(263) * 1.66053872801495e-12;
    //syn_T_C1.CD28_CD80, mw84eb780b_48c6_4fa5_b2be_4407efeafe58, index: 250
    //Unit: mole^(1)metre^(-2)
    _species_var[250] = PFILE(264) * 1.66053872801495e-12;
    //syn_T_C1.CD28_CD80_CD28, mwb87863cf_18e5_4971_80d1_47a71b69e0b4, index: 251
    //Unit: mole^(1)metre^(-2)
    _species_var[251] = PFILE(265) * 1.66053872801495e-12;
    //syn_T_C1.CD28_CD86, mw96eb5fed_a66b_4ac1_b3f5_f554b6c8a903, index: 252
    //Unit: mole^(1)metre^(-2)
    _species_var[252] = PFILE(266) * 1.66053872801495e-12;
    //syn_T_C1.CD80_CTLA4, mwf53ba71c_ccab_490a_b96d_47f5f8a14156, index: 253
    //Unit: mole^(1)metre^(-2)
    _species_var[253] = PFILE(267) * 1.66053872801495e-12;
    //syn_T_C1.CD80_CTLA4_CD80, mwe5786fd0_7ff6_430d_9c97_675a7ec3523f, index: 254
    //Unit: mole^(1)metre^(-2)
    _species_var[254] = PFILE(268) * 1.66053872801495e-12;
    //syn_T_C1.CTLA4_CD80_CTLA4, mw205dfb0d_b6bd_4a30_bc02_63afe6967196, index: 255
    //Unit: mole^(1)metre^(-2)
    _species_var[255] = PFILE(269) * 1.66053872801495e-12;
    //syn_T_C1.CD80_CTLA4_CD80_CTLA4, mw8d794e2d_2bb2_45f5_ab0a_b31db0602c78, index: 256
    //Unit: mole^(1)metre^(-2)
    _species_var[256] = PFILE(270) * 1.66053872801495e-12;
    //syn_T_C1.CD86_CTLA4, mw65785cb3_d21f_45ed_8fb3_446008197d69, index: 257
    //Unit: mole^(1)metre^(-2)
    _species_var[257] = PFILE(271) * 1.66053872801495e-12;
    //syn_T_C1.CD86_CTLA4_CD86, mwd6f8527b_0c08_4dde_9b71_4ada3422e01a, index: 258
    //Unit: mole^(1)metre^(-2)
    _species_var[258] = PFILE(272) * 1.66053872801495e-12;
    //syn_T_C1.TPDL1_CD80, mw85d41592_ecca_48a1_8e05_a25cfc5eab98, index: 259
    //Unit: mole^(1)metre^(-2)
    _species_var[259] = PFILE(273) * 1.66053872801495e-12;
    //syn_T_C1.TPDL1_CD80_TPDL1, mw89f0eb68_ffa3_40bd_9a0b_1d0608abd28c, index: 260
    //Unit: mole^(1)metre^(-2)
    _species_var[260] = PFILE(274) * 1.66053872801495e-12;
    //syn_T_C1.CD28, mw0847004a_1925_4aaa_9cba_3b724f25ad03, index: 261
    //Unit: mole^(1)metre^(-2)
    _species_var[261] = PFILE(275) * 1.66053872801495e-12;
    //syn_T_C1.CTLA4, mw15e77317_b36e_4bcb_81dc_7324a20f595b, index: 262
    //Unit: mole^(1)metre^(-2)
    _species_var[262] = PFILE(276) * 1.66053872801495e-12;
    //syn_T_C1.CD80, mwb01d4ed9_c7ba_42aa_a242_edb1fe194838, index: 263
    //Unit: mole^(1)metre^(-2)
    _species_var[263] = PFILE(277) * 1.66053872801495e-12;
    //syn_T_C1.CD86, mw1c2cd67b_fbac_4521_a619_ad14889ce948, index: 264
    //Unit: mole^(1)metre^(-2)
    _species_var[264] = PFILE(278) * 1.66053872801495e-12;
    //syn_T_C1.CTLA4_ipi, mw96ed0142_e96c_492b_93cf_acc179ee8316, index: 265
    //Unit: mole^(1)metre^(-2)
    _species_var[265] = PFILE(279) * 1.66053872801495e-12;
    //syn_T_C1.CTLA4_ipi_CTLA4, mwe822dcf7_e775_4a8a_b387_596c6b082128, index: 266
    //Unit: mole^(1)metre^(-2)
    _species_var[266] = PFILE(280) * 1.66053872801495e-12;
    //syn_T_APC.PD1_PDL1, mw41cbe269_f910_4bc7_9081_c272aea06d91, index: 267
    //Unit: mole^(1)metre^(-2)
    _species_var[267] = PFILE(281) * 1.66053872801495e-12;
    //syn_T_APC.PD1_PDL2, mw198a44cb_8ad7_4822_b12f_38b806b74f22, index: 268
    //Unit: mole^(1)metre^(-2)
    _species_var[268] = PFILE(282) * 1.66053872801495e-12;
    //syn_T_APC.PD1, mwf953b969_4cb1_4dd5_9b47_ea562a58fbc7, index: 269
    //Unit: mole^(1)metre^(-2)
    _species_var[269] = PFILE(283) * 1.66053872801495e-12;
    //syn_T_APC.PDL1, mw66b96129_5d7d_487f_b151_b13edce7ff64, index: 270
    //Unit: mole^(1)metre^(-2)
    _species_var[270] = PFILE(284) * 1.66053872801495e-12;
    //syn_T_APC.PDL2, mw55b4f1e0_fc97_4a8c_8580_63e478b405ed, index: 271
    //Unit: mole^(1)metre^(-2)
    _species_var[271] = PFILE(285) * 1.66053872801495e-12;
    //syn_T_APC.PD1_nivo, mwe446d0dd_738d_4bdb_ab45_b4516e41c288, index: 272
    //Unit: mole^(1)metre^(-2)
    _species_var[272] = PFILE(286) * 1.66053872801495e-12;
    //syn_T_APC.PD1_nivo_PD1, mw15f33334_f989_406a_a18e_89df7f8bf1aa, index: 273
    //Unit: mole^(1)metre^(-2)
    _species_var[273] = PFILE(287) * 1.66053872801495e-12;
    //syn_T_APC.PDL1_durv, mw4c4abe20_c04a_4125_8b8f_4685207cfc98, index: 274
    //Unit: mole^(1)metre^(-2)
    _species_var[274] = PFILE(288) * 1.66053872801495e-12;
    //syn_T_APC.PDL1_durv_PDL1, mwf10e7254_4a9c_42ae_8531_c111c5cf3b60, index: 275
    //Unit: mole^(1)metre^(-2)
    _species_var[275] = PFILE(289) * 1.66053872801495e-12;
    //syn_T_APC.TPDL1, mw31d01d0a_0f0b_4909_b89b_f5224b433898, index: 276
    //Unit: mole^(1)metre^(-2)
    _species_var[276] = PFILE(290) * 1.66053872801495e-12;
    //syn_T_APC.TPDL1_durv, mw139a37a6_fbae_433a_880b_8988b338bbe5, index: 277
    //Unit: mole^(1)metre^(-2)
    _species_var[277] = PFILE(291) * 1.66053872801495e-12;
    //syn_T_APC.TPDL1_durv_TPDL1, mweed652b9_58c2_4a40_9a97_7ea6ffdcefb4, index: 278
    //Unit: mole^(1)metre^(-2)
    _species_var[278] = PFILE(292) * 1.66053872801495e-12;
    //syn_T_APC.CD28_CD80, mw3d8bacd7_8629_457d_b2b0_e79f6b757232, index: 279
    //Unit: mole^(1)metre^(-2)
    _species_var[279] = PFILE(293) * 1.66053872801495e-12;
    //syn_T_APC.CD28_CD80_CD28, mwca557913_acda_479f_9ea2_24864b698cca, index: 280
    //Unit: mole^(1)metre^(-2)
    _species_var[280] = PFILE(294) * 1.66053872801495e-12;
    //syn_T_APC.CD28_CD86, mw681de61f_cd9e_4e2b_907c_5ad250aa048d, index: 281
    //Unit: mole^(1)metre^(-2)
    _species_var[281] = PFILE(295) * 1.66053872801495e-12;
    //syn_T_APC.CD80_CTLA4, mw22436f8a_23b4_4b01_b493_0b7d02740f1f, index: 282
    //Unit: mole^(1)metre^(-2)
    _species_var[282] = PFILE(296) * 1.66053872801495e-12;
    //syn_T_APC.CD80_CTLA4_CD80, mw87451fe8_c271_410f_bc50_83d57c6d881b, index: 283
    //Unit: mole^(1)metre^(-2)
    _species_var[283] = PFILE(297) * 1.66053872801495e-12;
    //syn_T_APC.CTLA4_CD80_CTLA4, mwa6ed6cd0_2cbe_4d33_85c9_71b137f0d1c1, index: 284
    //Unit: mole^(1)metre^(-2)
    _species_var[284] = PFILE(298) * 1.66053872801495e-12;
    //syn_T_APC.CD80_CTLA4_CD80_CTLA4, mwcb89e827_e971_4d66_9fab_1cdb06632993, index: 285
    //Unit: mole^(1)metre^(-2)
    _species_var[285] = PFILE(299) * 1.66053872801495e-12;
    //syn_T_APC.CD86_CTLA4, mwdb1bb2dd_1e72_4b0a_b20d_624a8f7190c3, index: 286
    //Unit: mole^(1)metre^(-2)
    _species_var[286] = PFILE(300) * 1.66053872801495e-12;
    //syn_T_APC.CD86_CTLA4_CD86, mw59a6d489_3420_48ce_ba0a_3878584e415c, index: 287
    //Unit: mole^(1)metre^(-2)
    _species_var[287] = PFILE(301) * 1.66053872801495e-12;
    //syn_T_APC.TPDL1_CD80, mw396cd5ad_4a36_4c0f_9981_84fbeafe0e9e, index: 288
    //Unit: mole^(1)metre^(-2)
    _species_var[288] = PFILE(302) * 1.66053872801495e-12;
    //syn_T_APC.TPDL1_CD80_TPDL1, mwb92e6ac5_9877_41de_b75a_0ca631ba7a3f, index: 289
    //Unit: mole^(1)metre^(-2)
    _species_var[289] = PFILE(303) * 1.66053872801495e-12;
    //syn_T_APC.CD28, mwbfebef30_7586_478a_82d9_dfb31937e6c2, index: 290
    //Unit: mole^(1)metre^(-2)
    _species_var[290] = PFILE(304) * 1.66053872801495e-12;
    //syn_T_APC.CTLA4, mw41e8e237_01f8_4f44_9a07_28de8d6e726f, index: 291
    //Unit: mole^(1)metre^(-2)
    _species_var[291] = PFILE(305) * 1.66053872801495e-12;
    //syn_T_APC.CD80, mwfdfbdf68_edf3_4119_8543_147010f24869, index: 292
    //Unit: mole^(1)metre^(-2)
    _species_var[292] = PFILE(306) * 1.66053872801495e-12;
    //syn_T_APC.CD86, mwabfd1a4e_d805_4c43_9e7b_eadec524f6d7, index: 293
    //Unit: mole^(1)metre^(-2)
    _species_var[293] = PFILE(307) * 1.66053872801495e-12;
    //syn_T_APC.CTLA4_ipi, mw02457974_d325_4851_b784_80cb51167a10, index: 294
    //Unit: mole^(1)metre^(-2)
    _species_var[294] = PFILE(308) * 1.66053872801495e-12;
    //syn_T_APC.CTLA4_ipi_CTLA4, mw56013dcf_c76b_4604_a5d3_a85937494906, index: 295
    //Unit: mole^(1)metre^(-2)
    _species_var[295] = PFILE(309) * 1.66053872801495e-12;
    
    return;
}

void ODE_system::setup_instance_tolerance(Param& param){

    //Tolerance
    realtype reltol = PFILE(3);
    realtype abstol_base = PFILE(4);
    N_Vector abstol = N_VNew_Serial(_neq);

    for (size_t i = 0; i < 296; i++)
    {
        NV_DATA_S(abstol)[i] = abstol_base * get_unit_conversion_species(i);
    }
    int flag = CVodeSVtolerances(_cvode_mem, reltol, abstol);
    check_flag(&flag, "CVodeSVtolerances", 1);

    
    return;
}

void ODE_system::eval_init_assignment(void){
    //Assignment Rules required before IA
    //InitialAssignment
    _species_var[240] = _class_parameter[562] / _class_parameter[264];
    _species_var[261] = _class_parameter[563] / _class_parameter[264];
    _species_var[262] = _class_parameter[564] / _class_parameter[263];
    _species_var[247] = _class_parameter[565] / _class_parameter[264];
    _species_var[241] = _class_parameter[566] / _class_parameter[265];
    _species_var[242] = _class_parameter[567] / _class_parameter[265];
    _species_var[263] = _class_parameter[568] / _class_parameter[265];
    _species_var[264] = _class_parameter[569] / _class_parameter[265];
    _species_var[269] = _class_parameter[562] / _class_parameter[264];
    _species_var[290] = _class_parameter[563] / _class_parameter[264];
    _species_var[291] = _class_parameter[564] / _class_parameter[263];
    _species_var[276] = _class_parameter[565] / _class_parameter[264];
    _species_var[270] = _class_parameter[570] / _class_parameter[266];
    _species_var[271] = _class_parameter[571] / _class_parameter[266];
    _species_var[292] = _class_parameter[572] / _class_parameter[266];
    _species_var[293] = _class_parameter[573] / _class_parameter[266];
    _species_var[77] = _class_parameter[574];
    _species_var[46] = _class_parameter[574];

    updateVar();
    
    return;
}
void ODE_system::setupEvents(void){

    _nevent = 1;
    _nroot = 1;

    _trigger_element_type = std::vector<EVENT_TRIGGER_ELEM_TYPE>(_nroot, TRIGGER_NON_INSTANT);
    _trigger_element_satisfied = std::vector<bool>(_nroot, false);
    _event_triggered = std::vector<bool>(_nevent, false);

    //V_T.C1 < (0.9 * cell)
    _trigger_element_type[0] = TRIGGER_NON_INSTANT;

    _event_triggered[0] = true;

    return;
}
int ODE_system::f(realtype t, N_Vector y, N_Vector ydot, void *user_data){

    ODE_system* ptrOde = static_cast<ODE_system*>(user_data);

    //Assignment rules:

    realtype AUX_VAR_V_T = PARAM(13) + PARAM(11) * SPVAR(50) + PARAM(12) * SPVAR(51) + PARAM(11) * SPVAR(52) + PARAM(12) * SPVAR(53) + PARAM(12) * SPVAR(54) + PARAM(12) * SPVAR(55) + PARAM(12) * SPVAR(56) + PARAM(12) * SPVAR(57) + PARAM(12) * SPVAR(58) + PARAM(12) * SPVAR(59) + PARAM(12) * SPVAR(60) + PARAM(12) * SPVAR(61) + PARAM(12) * SPVAR(62) + PARAM(12) * SPVAR(63) + PARAM(12) * SPVAR(64) + PARAM(12) * SPVAR(65) + PARAM(12) * SPVAR(66) + PARAM(12) * SPVAR(67) + PARAM(12) * SPVAR(68) + PARAM(12) * SPVAR(69) + PARAM(12) * SPVAR(70);

    realtype AUX_VAR_C_total = 0.0 * PARAM(9) + SPVAR(52);

    realtype AUX_VAR_T_total = 0.0 * PARAM(9) + SPVAR(53) + SPVAR(54) + SPVAR(55) + SPVAR(56) + SPVAR(57) + SPVAR(58) + SPVAR(59) + SPVAR(60) + SPVAR(61) + SPVAR(62) + SPVAR(63) + SPVAR(64) + SPVAR(65) + SPVAR(66) + SPVAR(67) + SPVAR(68) + SPVAR(69) + SPVAR(70);

    realtype AUX_VAR_T_total_LN = 0.0 * PARAM(9) + SPVAR(91) + SPVAR(94) + SPVAR(97) + SPVAR(100) + SPVAR(103) + SPVAR(106) + SPVAR(109) + SPVAR(112) + SPVAR(115) + SPVAR(118) + SPVAR(121) + SPVAR(124) + SPVAR(127) + SPVAR(130) + SPVAR(133) + SPVAR(136) + SPVAR(139);

    realtype AUX_VAR_H_PD1_C1 = std::pow((SPVAR(238) + SPVAR(239)) / PARAM(558), PARAM(559)) / (std::pow((SPVAR(238) + SPVAR(239)) / PARAM(558), PARAM(559)) + 1.0);

    realtype AUX_VAR_H_MDSC_C1 = 1.0 - (1.0 - SPVAR(82) / (PARAM(590) + SPVAR(82))) * (1.0 - SPVAR(83) / (PARAM(589) + SPVAR(83))) * (1.0 - SPVAR(84) / (SPVAR(84) + PARAM(597)));

    realtype AUX_VAR_Tregs_ = SPVAR(53);

    realtype AUX_VAR_H_CD28_APC = std::pow((SPVAR(279) + SPVAR(281) + 2.0 * SPVAR(280)) / PARAM(560), PARAM(561)) / (std::pow((SPVAR(279) + SPVAR(281) + 2.0 * SPVAR(280)) / PARAM(560), PARAM(561)) + 1.0);

    realtype AUX_VAR_pTCR_p0_MHC_tot = 0.5 * (SPVAR(220) / PARAM(19) + PARAM(269) + PARAM(268) / PARAM(267) - PARAM(269) * std::pow(std::pow((SPVAR(220) / PARAM(19) + PARAM(269) + PARAM(268) / PARAM(267)) / PARAM(269), 2.0) - 4.0 * SPVAR(220) / PARAM(19) / PARAM(269), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p1_MHC_tot = PARAM(279) / (PARAM(279) + PARAM(281)) * std::pow(PARAM(280) / (PARAM(279) + PARAM(280)), PARAM(282)) * 0.5 * (SPVAR(221) / PARAM(38) + PARAM(283) + PARAM(279) / PARAM(278) - PARAM(283) * std::pow(std::pow((SPVAR(221) / PARAM(38) + PARAM(283) + PARAM(279) / PARAM(278)) / PARAM(283), 2.0) - 4.0 * SPVAR(221) / PARAM(38) / PARAM(283), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p2_MHC_tot = PARAM(293) / (PARAM(293) + PARAM(295)) * std::pow(PARAM(294) / (PARAM(293) + PARAM(294)), PARAM(296)) * 0.5 * (SPVAR(222) / PARAM(50) + PARAM(297) + PARAM(293) / PARAM(292) - PARAM(297) * std::pow(std::pow((SPVAR(222) / PARAM(50) + PARAM(297) + PARAM(293) / PARAM(292)) / PARAM(297), 2.0) - 4.0 * SPVAR(222) / PARAM(50) / PARAM(297), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p3_MHC_tot = PARAM(307) / (PARAM(307) + PARAM(309)) * std::pow(PARAM(308) / (PARAM(307) + PARAM(308)), PARAM(310)) * 0.5 * (SPVAR(223) / PARAM(62) + PARAM(311) + PARAM(307) / PARAM(306) - PARAM(311) * std::pow(std::pow((SPVAR(223) / PARAM(62) + PARAM(311) + PARAM(307) / PARAM(306)) / PARAM(311), 2.0) - 4.0 * SPVAR(223) / PARAM(62) / PARAM(311), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p4_MHC_tot = PARAM(321) / (PARAM(321) + PARAM(323)) * std::pow(PARAM(322) / (PARAM(321) + PARAM(322)), PARAM(324)) * 0.5 * (SPVAR(224) / PARAM(74) + PARAM(325) + PARAM(321) / PARAM(320) - PARAM(325) * std::pow(std::pow((SPVAR(224) / PARAM(74) + PARAM(325) + PARAM(321) / PARAM(320)) / PARAM(325), 2.0) - 4.0 * SPVAR(224) / PARAM(74) / PARAM(325), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p5_MHC_tot = PARAM(335) / (PARAM(335) + PARAM(337)) * std::pow(PARAM(336) / (PARAM(335) + PARAM(336)), PARAM(338)) * 0.5 * (SPVAR(225) / PARAM(86) + PARAM(339) + PARAM(335) / PARAM(334) - PARAM(339) * std::pow(std::pow((SPVAR(225) / PARAM(86) + PARAM(339) + PARAM(335) / PARAM(334)) / PARAM(339), 2.0) - 4.0 * SPVAR(225) / PARAM(86) / PARAM(339), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p6_MHC_tot = PARAM(349) / (PARAM(349) + PARAM(351)) * std::pow(PARAM(350) / (PARAM(349) + PARAM(350)), PARAM(352)) * 0.5 * (SPVAR(226) / PARAM(98) + PARAM(353) + PARAM(349) / PARAM(348) - PARAM(353) * std::pow(std::pow((SPVAR(226) / PARAM(98) + PARAM(353) + PARAM(349) / PARAM(348)) / PARAM(353), 2.0) - 4.0 * SPVAR(226) / PARAM(98) / PARAM(353), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p7_MHC_tot = PARAM(363) / (PARAM(363) + PARAM(365)) * std::pow(PARAM(364) / (PARAM(363) + PARAM(364)), PARAM(366)) * 0.5 * (SPVAR(227) / PARAM(110) + PARAM(367) + PARAM(363) / PARAM(362) - PARAM(367) * std::pow(std::pow((SPVAR(227) / PARAM(110) + PARAM(367) + PARAM(363) / PARAM(362)) / PARAM(367), 2.0) - 4.0 * SPVAR(227) / PARAM(110) / PARAM(367), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p8_MHC_tot = PARAM(377) / (PARAM(377) + PARAM(379)) * std::pow(PARAM(378) / (PARAM(377) + PARAM(378)), PARAM(380)) * 0.5 * (SPVAR(228) / PARAM(122) + PARAM(381) + PARAM(377) / PARAM(376) - PARAM(381) * std::pow(std::pow((SPVAR(228) / PARAM(122) + PARAM(381) + PARAM(377) / PARAM(376)) / PARAM(381), 2.0) - 4.0 * SPVAR(228) / PARAM(122) / PARAM(381), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p9_MHC_tot = PARAM(391) / (PARAM(391) + PARAM(393)) * std::pow(PARAM(392) / (PARAM(391) + PARAM(392)), PARAM(394)) * 0.5 * (SPVAR(229) / PARAM(134) + PARAM(395) + PARAM(391) / PARAM(390) - PARAM(395) * std::pow(std::pow((SPVAR(229) / PARAM(134) + PARAM(395) + PARAM(391) / PARAM(390)) / PARAM(395), 2.0) - 4.0 * SPVAR(229) / PARAM(134) / PARAM(395), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p10_MHC_tot = PARAM(405) / (PARAM(405) + PARAM(407)) * std::pow(PARAM(406) / (PARAM(405) + PARAM(406)), PARAM(408)) * 0.5 * (SPVAR(230) / PARAM(146) + PARAM(409) + PARAM(405) / PARAM(404) - PARAM(409) * std::pow(std::pow((SPVAR(230) / PARAM(146) + PARAM(409) + PARAM(405) / PARAM(404)) / PARAM(409), 2.0) - 4.0 * SPVAR(230) / PARAM(146) / PARAM(409), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p11_MHC_tot = PARAM(419) / (PARAM(419) + PARAM(421)) * std::pow(PARAM(420) / (PARAM(419) + PARAM(420)), PARAM(422)) * 0.5 * (SPVAR(231) / PARAM(158) + PARAM(423) + PARAM(419) / PARAM(418) - PARAM(423) * std::pow(std::pow((SPVAR(231) / PARAM(158) + PARAM(423) + PARAM(419) / PARAM(418)) / PARAM(423), 2.0) - 4.0 * SPVAR(231) / PARAM(158) / PARAM(423), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p12_MHC_tot = PARAM(433) / (PARAM(433) + PARAM(435)) * std::pow(PARAM(434) / (PARAM(433) + PARAM(434)), PARAM(436)) * 0.5 * (SPVAR(232) / PARAM(170) + PARAM(437) + PARAM(433) / PARAM(432) - PARAM(437) * std::pow(std::pow((SPVAR(232) / PARAM(170) + PARAM(437) + PARAM(433) / PARAM(432)) / PARAM(437), 2.0) - 4.0 * SPVAR(232) / PARAM(170) / PARAM(437), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p13_MHC_tot = PARAM(447) / (PARAM(447) + PARAM(449)) * std::pow(PARAM(448) / (PARAM(447) + PARAM(448)), PARAM(450)) * 0.5 * (SPVAR(233) / PARAM(182) + PARAM(451) + PARAM(447) / PARAM(446) - PARAM(451) * std::pow(std::pow((SPVAR(233) / PARAM(182) + PARAM(451) + PARAM(447) / PARAM(446)) / PARAM(451), 2.0) - 4.0 * SPVAR(233) / PARAM(182) / PARAM(451), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p14_MHC_tot = PARAM(461) / (PARAM(461) + PARAM(463)) * std::pow(PARAM(462) / (PARAM(461) + PARAM(462)), PARAM(464)) * 0.5 * (SPVAR(234) / PARAM(194) + PARAM(465) + PARAM(461) / PARAM(460) - PARAM(465) * std::pow(std::pow((SPVAR(234) / PARAM(194) + PARAM(465) + PARAM(461) / PARAM(460)) / PARAM(465), 2.0) - 4.0 * SPVAR(234) / PARAM(194) / PARAM(465), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p15_MHC_tot = PARAM(475) / (PARAM(475) + PARAM(477)) * std::pow(PARAM(476) / (PARAM(475) + PARAM(476)), PARAM(478)) * 0.5 * (SPVAR(235) / PARAM(206) + PARAM(479) + PARAM(475) / PARAM(474) - PARAM(479) * std::pow(std::pow((SPVAR(235) / PARAM(206) + PARAM(479) + PARAM(475) / PARAM(474)) / PARAM(479), 2.0) - 4.0 * SPVAR(235) / PARAM(206) / PARAM(479), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p16_MHC_tot = PARAM(489) / (PARAM(489) + PARAM(491)) * std::pow(PARAM(490) / (PARAM(489) + PARAM(490)), PARAM(492)) * 0.5 * (SPVAR(236) / PARAM(218) + PARAM(493) + PARAM(489) / PARAM(488) - PARAM(493) * std::pow(std::pow((SPVAR(236) / PARAM(218) + PARAM(493) + PARAM(489) / PARAM(488)) / PARAM(493), 2.0) - 4.0 * SPVAR(236) / PARAM(218) / PARAM(493), 1.0 / 2.0));

    realtype AUX_VAR_pTCR_p17_MHC_tot = PARAM(503) / (PARAM(503) + PARAM(505)) * std::pow(PARAM(504) / (PARAM(503) + PARAM(504)), PARAM(506)) * 0.5 * (SPVAR(237) / PARAM(230) + PARAM(507) + PARAM(503) / PARAM(502) - PARAM(507) * std::pow(std::pow((SPVAR(237) / PARAM(230) + PARAM(507) + PARAM(503) / PARAM(502)) / PARAM(507), 2.0) - 4.0 * SPVAR(237) / PARAM(230) / PARAM(507), 1.0 / 2.0));

    realtype AUX_VAR_H_CD28_C1 = std::pow((SPVAR(250) + SPVAR(252) + 2.0 * SPVAR(251)) / PARAM(560), PARAM(561)) / (std::pow((SPVAR(250) + SPVAR(252) + 2.0 * SPVAR(251)) / PARAM(560), PARAM(561)) + 1.0);

    realtype AUX_VAR_H_PD1_APC = std::pow((SPVAR(267) + SPVAR(268)) / PARAM(558), PARAM(559)) / (std::pow((SPVAR(267) + SPVAR(268)) / PARAM(558), PARAM(559)) + 1.0);

    realtype AUX_VAR_H_Treg_T = std::pow((SPVAR(78) + 2.0 * SPVAR(79)) / PARAM(575), PARAM(576)) / (std::pow((SPVAR(78) + 2.0 * SPVAR(79)) / PARAM(575), PARAM(576)) + 1.0);

    realtype AUX_VAR_H_Treg_P = std::pow((SPVAR(47) + 2.0 * SPVAR(48)) / PARAM(575), PARAM(576)) / (std::pow((SPVAR(47) + 2.0 * SPVAR(48)) / PARAM(575), PARAM(576)) + 1.0);

    realtype AUX_VAR_H_ENT_C1 = SPVAR(84) / (SPVAR(84) + PARAM(580));

    realtype AUX_VAR_H_APC = PARAM(252) * SPVAR(140) / (PARAM(252) * SPVAR(140) + AUX_VAR_T_total_LN + PARAM(9));

    realtype AUX_VAR_H_mAPC = PARAM(252) * SPVAR(141) / (PARAM(252) * SPVAR(141) + AUX_VAR_T_total_LN + PARAM(9));

    realtype AUX_VAR_R_Tcell = 0.0 * PARAM(9) / PARAM(10) + PARAM(49) * SPVAR(54) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(61) * SPVAR(55) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(73) * SPVAR(56) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(85) * SPVAR(57) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(97) * SPVAR(58) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(109) * SPVAR(59) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(121) * SPVAR(60) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(133) * SPVAR(61) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(145) * SPVAR(62) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(157) * SPVAR(63) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(169) * SPVAR(64) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(181) * SPVAR(65) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(193) * SPVAR(66) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(205) * SPVAR(67) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(217) * SPVAR(68) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(229) * SPVAR(69) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1) + PARAM(241) * SPVAR(70) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype AUX_VAR_N_aT = PARAM(34) + PARAM(35) * AUX_VAR_H_CD28_APC + PARAM(36) * SPVAR(88) / (PARAM(32) + SPVAR(88));

    realtype AUX_VAR_H_P0 = AUX_VAR_pTCR_p0_MHC_tot / (AUX_VAR_pTCR_p0_MHC_tot + PARAM(261));

    realtype AUX_VAR_H_P1 = AUX_VAR_pTCR_p1_MHC_tot / (AUX_VAR_pTCR_p1_MHC_tot + PARAM(276));

    realtype AUX_VAR_H_P2 = AUX_VAR_pTCR_p2_MHC_tot / (AUX_VAR_pTCR_p2_MHC_tot + PARAM(290));

    realtype AUX_VAR_H_P3 = AUX_VAR_pTCR_p3_MHC_tot / (AUX_VAR_pTCR_p3_MHC_tot + PARAM(304));

    realtype AUX_VAR_H_P4 = AUX_VAR_pTCR_p4_MHC_tot / (AUX_VAR_pTCR_p4_MHC_tot + PARAM(318));

    realtype AUX_VAR_H_P5 = AUX_VAR_pTCR_p5_MHC_tot / (AUX_VAR_pTCR_p5_MHC_tot + PARAM(332));

    realtype AUX_VAR_H_P6 = AUX_VAR_pTCR_p6_MHC_tot / (AUX_VAR_pTCR_p6_MHC_tot + PARAM(346));

    realtype AUX_VAR_H_P7 = AUX_VAR_pTCR_p7_MHC_tot / (AUX_VAR_pTCR_p7_MHC_tot + PARAM(360));

    realtype AUX_VAR_H_P8 = AUX_VAR_pTCR_p8_MHC_tot / (AUX_VAR_pTCR_p8_MHC_tot + PARAM(374));

    realtype AUX_VAR_H_P9 = AUX_VAR_pTCR_p9_MHC_tot / (AUX_VAR_pTCR_p9_MHC_tot + PARAM(388));

    realtype AUX_VAR_H_P10 = AUX_VAR_pTCR_p10_MHC_tot / (AUX_VAR_pTCR_p10_MHC_tot + PARAM(402));

    realtype AUX_VAR_H_P11 = AUX_VAR_pTCR_p11_MHC_tot / (AUX_VAR_pTCR_p11_MHC_tot + PARAM(416));

    realtype AUX_VAR_H_P12 = AUX_VAR_pTCR_p12_MHC_tot / (AUX_VAR_pTCR_p12_MHC_tot + PARAM(430));

    realtype AUX_VAR_H_P13 = AUX_VAR_pTCR_p13_MHC_tot / (AUX_VAR_pTCR_p13_MHC_tot + PARAM(444));

    realtype AUX_VAR_H_P14 = AUX_VAR_pTCR_p14_MHC_tot / (AUX_VAR_pTCR_p14_MHC_tot + PARAM(458));

    realtype AUX_VAR_H_P15 = AUX_VAR_pTCR_p15_MHC_tot / (AUX_VAR_pTCR_p15_MHC_tot + PARAM(472));

    realtype AUX_VAR_H_P16 = AUX_VAR_pTCR_p16_MHC_tot / (AUX_VAR_pTCR_p16_MHC_tot + PARAM(486));

    realtype AUX_VAR_H_P17 = AUX_VAR_pTCR_p17_MHC_tot / (AUX_VAR_pTCR_p17_MHC_tot + PARAM(500));

    //Reaction fluxes:

    realtype ReactionFlux1 = PARAM(8) * SPVAR(50);

    realtype ReactionFlux2 = PARAM(8) * SPVAR(51);

    realtype ReactionFlux3 = PARAM(14) * SPVAR(52) * (1.0 - AUX_VAR_C_total / PARAM(15)) * (1.0 - AUX_VAR_H_ENT_C1);

    realtype ReactionFlux4 = PARAM(16) * SPVAR(52);

    realtype ReactionFlux5 = PARAM(20) * PARAM(19);

    realtype ReactionFlux6 = PARAM(21) * SPVAR(85);

    realtype ReactionFlux7 = PARAM(22) * AUX_VAR_H_APC * AUX_VAR_H_P0 * SPVAR(85);

    realtype ReactionFlux8 = PARAM(23) / AUX_VAR_N_aT * SPVAR(86);

    realtype ReactionFlux9 = PARAM(23) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(86);

    realtype ReactionFlux10 = PARAM(24) * SPVAR(0);

    realtype ReactionFlux11 = PARAM(24) * SPVAR(25);

    realtype ReactionFlux12 = PARAM(24) * SPVAR(53);

    realtype ReactionFlux13 = PARAM(24) * SPVAR(87);

    realtype ReactionFlux14 = PARAM(25) * SPVAR(0);

    realtype ReactionFlux15 = PARAM(26) * SPVAR(25);

    realtype ReactionFlux16 = PARAM(27) * AUX_VAR_V_T * SPVAR(0) * (1.0 - SPVAR(53) / (PARAM(594) * AUX_VAR_V_T));

    realtype ReactionFlux17 = PARAM(28) * SPVAR(87);

    realtype ReactionFlux18 = PARAM(29) * SPVAR(88) * PARAM(2);

    realtype ReactionFlux19 = PARAM(30) * AUX_VAR_T_total_LN * SPVAR(88) / (PARAM(32) + SPVAR(88));

    realtype ReactionFlux20 = PARAM(30) * SPVAR(87) * SPVAR(88) / (PARAM(33) + SPVAR(88));

    realtype ReactionFlux21 = PARAM(39) * PARAM(38);

    realtype ReactionFlux22 = PARAM(40) * SPVAR(89);

    realtype ReactionFlux23 = PARAM(41) * AUX_VAR_H_mAPC * AUX_VAR_H_P1 * SPVAR(89);

    realtype ReactionFlux24 = PARAM(42) / AUX_VAR_N_aT * SPVAR(90);

    realtype ReactionFlux25 = PARAM(42) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(90);

    realtype ReactionFlux26 = PARAM(43) * SPVAR(1);

    realtype ReactionFlux27 = PARAM(43) * SPVAR(26);

    realtype ReactionFlux28 = PARAM(43) * SPVAR(54);

    realtype ReactionFlux29 = PARAM(43) * SPVAR(91);

    realtype ReactionFlux30 = PARAM(37) * SPVAR(54) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux31 = PARAM(48) * SPVAR(54) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux32 = PARAM(44) * SPVAR(1);

    realtype ReactionFlux33 = PARAM(45) * SPVAR(26);

    realtype ReactionFlux34 = PARAM(46) * AUX_VAR_V_T * SPVAR(1);

    realtype ReactionFlux35 = PARAM(47) * SPVAR(91);

    realtype ReactionFlux36 = PARAM(31) * SPVAR(90);

    realtype ReactionFlux37 = PARAM(49) * SPVAR(54) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux38 = PARAM(51) * PARAM(50);

    realtype ReactionFlux39 = PARAM(52) * SPVAR(92);

    realtype ReactionFlux40 = PARAM(53) * AUX_VAR_H_mAPC * AUX_VAR_H_P2 * SPVAR(92);

    realtype ReactionFlux41 = PARAM(54) / AUX_VAR_N_aT * SPVAR(93);

    realtype ReactionFlux42 = PARAM(54) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(93);

    realtype ReactionFlux43 = PARAM(55) * SPVAR(2);

    realtype ReactionFlux44 = PARAM(55) * SPVAR(27);

    realtype ReactionFlux45 = PARAM(55) * SPVAR(55);

    realtype ReactionFlux46 = PARAM(55) * SPVAR(94);

    realtype ReactionFlux47 = PARAM(37) * SPVAR(55) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux48 = PARAM(60) * SPVAR(55) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux49 = PARAM(56) * SPVAR(2);

    realtype ReactionFlux50 = PARAM(57) * SPVAR(27);

    realtype ReactionFlux51 = PARAM(58) * AUX_VAR_V_T * SPVAR(2);

    realtype ReactionFlux52 = PARAM(59) * SPVAR(94);

    realtype ReactionFlux53 = PARAM(31) * SPVAR(93);

    realtype ReactionFlux54 = PARAM(61) * SPVAR(55) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux55 = PARAM(63) * PARAM(62);

    realtype ReactionFlux56 = PARAM(64) * SPVAR(95);

    realtype ReactionFlux57 = PARAM(65) * AUX_VAR_H_mAPC * AUX_VAR_H_P3 * SPVAR(95);

    realtype ReactionFlux58 = PARAM(66) / AUX_VAR_N_aT * SPVAR(96);

    realtype ReactionFlux59 = PARAM(66) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(96);

    realtype ReactionFlux60 = PARAM(67) * SPVAR(3);

    realtype ReactionFlux61 = PARAM(67) * SPVAR(28);

    realtype ReactionFlux62 = PARAM(67) * SPVAR(56);

    realtype ReactionFlux63 = PARAM(67) * SPVAR(97);

    realtype ReactionFlux64 = PARAM(37) * SPVAR(56) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux65 = PARAM(72) * SPVAR(56) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux66 = PARAM(68) * SPVAR(3);

    realtype ReactionFlux67 = PARAM(69) * SPVAR(28);

    realtype ReactionFlux68 = PARAM(70) * AUX_VAR_V_T * SPVAR(3);

    realtype ReactionFlux69 = PARAM(71) * SPVAR(97);

    realtype ReactionFlux70 = PARAM(31) * SPVAR(96);

    realtype ReactionFlux71 = PARAM(73) * SPVAR(56) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux72 = PARAM(75) * PARAM(74);

    realtype ReactionFlux73 = PARAM(76) * SPVAR(98);

    realtype ReactionFlux74 = PARAM(77) * AUX_VAR_H_mAPC * AUX_VAR_H_P4 * SPVAR(98);

    realtype ReactionFlux75 = PARAM(78) / AUX_VAR_N_aT * SPVAR(99);

    realtype ReactionFlux76 = PARAM(78) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(99);

    realtype ReactionFlux77 = PARAM(79) * SPVAR(4);

    realtype ReactionFlux78 = PARAM(79) * SPVAR(29);

    realtype ReactionFlux79 = PARAM(79) * SPVAR(57);

    realtype ReactionFlux80 = PARAM(79) * SPVAR(100);

    realtype ReactionFlux81 = PARAM(37) * SPVAR(57) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux82 = PARAM(84) * SPVAR(57) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux83 = PARAM(80) * SPVAR(4);

    realtype ReactionFlux84 = PARAM(81) * SPVAR(29);

    realtype ReactionFlux85 = PARAM(82) * AUX_VAR_V_T * SPVAR(4);

    realtype ReactionFlux86 = PARAM(83) * SPVAR(100);

    realtype ReactionFlux87 = PARAM(31) * SPVAR(99);

    realtype ReactionFlux88 = PARAM(85) * SPVAR(57) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux89 = PARAM(87) * PARAM(86);

    realtype ReactionFlux90 = PARAM(88) * SPVAR(101);

    realtype ReactionFlux91 = PARAM(89) * AUX_VAR_H_mAPC * AUX_VAR_H_P5 * SPVAR(101);

    realtype ReactionFlux92 = PARAM(90) / AUX_VAR_N_aT * SPVAR(102);

    realtype ReactionFlux93 = PARAM(90) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(102);

    realtype ReactionFlux94 = PARAM(91) * SPVAR(5);

    realtype ReactionFlux95 = PARAM(91) * SPVAR(30);

    realtype ReactionFlux96 = PARAM(91) * SPVAR(58);

    realtype ReactionFlux97 = PARAM(91) * SPVAR(103);

    realtype ReactionFlux98 = PARAM(37) * SPVAR(58) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux99 = PARAM(96) * SPVAR(58) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux100 = PARAM(92) * SPVAR(5);

    realtype ReactionFlux101 = PARAM(93) * SPVAR(30);

    realtype ReactionFlux102 = PARAM(94) * AUX_VAR_V_T * SPVAR(5);

    realtype ReactionFlux103 = PARAM(95) * SPVAR(103);

    realtype ReactionFlux104 = PARAM(31) * SPVAR(102);

    realtype ReactionFlux105 = PARAM(97) * SPVAR(58) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux106 = PARAM(99) * PARAM(98);

    realtype ReactionFlux107 = PARAM(100) * SPVAR(104);

    realtype ReactionFlux108 = PARAM(101) * AUX_VAR_H_mAPC * AUX_VAR_H_P6 * SPVAR(104);

    realtype ReactionFlux109 = PARAM(102) / AUX_VAR_N_aT * SPVAR(105);

    realtype ReactionFlux110 = PARAM(102) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(105);

    realtype ReactionFlux111 = PARAM(103) * SPVAR(6);

    realtype ReactionFlux112 = PARAM(103) * SPVAR(31);

    realtype ReactionFlux113 = PARAM(103) * SPVAR(59);

    realtype ReactionFlux114 = PARAM(103) * SPVAR(106);

    realtype ReactionFlux115 = PARAM(37) * SPVAR(59) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux116 = PARAM(108) * SPVAR(59) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux117 = PARAM(104) * SPVAR(6);

    realtype ReactionFlux118 = PARAM(105) * SPVAR(31);

    realtype ReactionFlux119 = PARAM(106) * AUX_VAR_V_T * SPVAR(6);

    realtype ReactionFlux120 = PARAM(107) * SPVAR(106);

    realtype ReactionFlux121 = PARAM(31) * SPVAR(105);

    realtype ReactionFlux122 = PARAM(109) * SPVAR(59) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux123 = PARAM(111) * PARAM(110);

    realtype ReactionFlux124 = PARAM(112) * SPVAR(107);

    realtype ReactionFlux125 = PARAM(113) * AUX_VAR_H_mAPC * AUX_VAR_H_P7 * SPVAR(107);

    realtype ReactionFlux126 = PARAM(114) / AUX_VAR_N_aT * SPVAR(108);

    realtype ReactionFlux127 = PARAM(114) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(108);

    realtype ReactionFlux128 = PARAM(115) * SPVAR(7);

    realtype ReactionFlux129 = PARAM(115) * SPVAR(32);

    realtype ReactionFlux130 = PARAM(115) * SPVAR(60);

    realtype ReactionFlux131 = PARAM(115) * SPVAR(109);

    realtype ReactionFlux132 = PARAM(37) * SPVAR(60) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux133 = PARAM(120) * SPVAR(60) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux134 = PARAM(116) * SPVAR(7);

    realtype ReactionFlux135 = PARAM(117) * SPVAR(32);

    realtype ReactionFlux136 = PARAM(118) * AUX_VAR_V_T * SPVAR(7);

    realtype ReactionFlux137 = PARAM(119) * SPVAR(109);

    realtype ReactionFlux138 = PARAM(31) * SPVAR(108);

    realtype ReactionFlux139 = PARAM(121) * SPVAR(60) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux140 = PARAM(123) * PARAM(122);

    realtype ReactionFlux141 = PARAM(124) * SPVAR(110);

    realtype ReactionFlux142 = PARAM(125) * AUX_VAR_H_mAPC * AUX_VAR_H_P8 * SPVAR(110);

    realtype ReactionFlux143 = PARAM(126) / AUX_VAR_N_aT * SPVAR(111);

    realtype ReactionFlux144 = PARAM(126) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(111);

    realtype ReactionFlux145 = PARAM(127) * SPVAR(8);

    realtype ReactionFlux146 = PARAM(127) * SPVAR(33);

    realtype ReactionFlux147 = PARAM(127) * SPVAR(61);

    realtype ReactionFlux148 = PARAM(127) * SPVAR(112);

    realtype ReactionFlux149 = PARAM(37) * SPVAR(61) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux150 = PARAM(132) * SPVAR(61) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux151 = PARAM(128) * SPVAR(8);

    realtype ReactionFlux152 = PARAM(129) * SPVAR(33);

    realtype ReactionFlux153 = PARAM(130) * AUX_VAR_V_T * SPVAR(8);

    realtype ReactionFlux154 = PARAM(131) * SPVAR(112);

    realtype ReactionFlux155 = PARAM(31) * SPVAR(111);

    realtype ReactionFlux156 = PARAM(133) * SPVAR(61) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux157 = PARAM(135) * PARAM(134);

    realtype ReactionFlux158 = PARAM(136) * SPVAR(113);

    realtype ReactionFlux159 = PARAM(137) * AUX_VAR_H_mAPC * AUX_VAR_H_P9 * SPVAR(113);

    realtype ReactionFlux160 = PARAM(138) / AUX_VAR_N_aT * SPVAR(114);

    realtype ReactionFlux161 = PARAM(138) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(114);

    realtype ReactionFlux162 = PARAM(139) * SPVAR(9);

    realtype ReactionFlux163 = PARAM(139) * SPVAR(34);

    realtype ReactionFlux164 = PARAM(139) * SPVAR(62);

    realtype ReactionFlux165 = PARAM(139) * SPVAR(115);

    realtype ReactionFlux166 = PARAM(37) * SPVAR(62) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux167 = PARAM(144) * SPVAR(62) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux168 = PARAM(140) * SPVAR(9);

    realtype ReactionFlux169 = PARAM(141) * SPVAR(34);

    realtype ReactionFlux170 = PARAM(142) * AUX_VAR_V_T * SPVAR(9);

    realtype ReactionFlux171 = PARAM(143) * SPVAR(115);

    realtype ReactionFlux172 = PARAM(31) * SPVAR(114);

    realtype ReactionFlux173 = PARAM(145) * SPVAR(62) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux174 = PARAM(147) * PARAM(146);

    realtype ReactionFlux175 = PARAM(148) * SPVAR(116);

    realtype ReactionFlux176 = PARAM(149) * AUX_VAR_H_mAPC * AUX_VAR_H_P10 * SPVAR(116);

    realtype ReactionFlux177 = PARAM(150) / AUX_VAR_N_aT * SPVAR(117);

    realtype ReactionFlux178 = PARAM(150) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(117);

    realtype ReactionFlux179 = PARAM(151) * SPVAR(10);

    realtype ReactionFlux180 = PARAM(151) * SPVAR(35);

    realtype ReactionFlux181 = PARAM(151) * SPVAR(63);

    realtype ReactionFlux182 = PARAM(151) * SPVAR(118);

    realtype ReactionFlux183 = PARAM(37) * SPVAR(63) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux184 = PARAM(156) * SPVAR(63) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux185 = PARAM(152) * SPVAR(10);

    realtype ReactionFlux186 = PARAM(153) * SPVAR(35);

    realtype ReactionFlux187 = PARAM(154) * AUX_VAR_V_T * SPVAR(10);

    realtype ReactionFlux188 = PARAM(155) * SPVAR(118);

    realtype ReactionFlux189 = PARAM(31) * SPVAR(117);

    realtype ReactionFlux190 = PARAM(157) * SPVAR(63) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux191 = PARAM(159) * PARAM(158);

    realtype ReactionFlux192 = PARAM(160) * SPVAR(119);

    realtype ReactionFlux193 = PARAM(161) * AUX_VAR_H_mAPC * AUX_VAR_H_P11 * SPVAR(119);

    realtype ReactionFlux194 = PARAM(162) / AUX_VAR_N_aT * SPVAR(120);

    realtype ReactionFlux195 = PARAM(162) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(120);

    realtype ReactionFlux196 = PARAM(163) * SPVAR(11);

    realtype ReactionFlux197 = PARAM(163) * SPVAR(36);

    realtype ReactionFlux198 = PARAM(163) * SPVAR(64);

    realtype ReactionFlux199 = PARAM(163) * SPVAR(121);

    realtype ReactionFlux200 = PARAM(37) * SPVAR(64) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux201 = PARAM(168) * SPVAR(64) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux202 = PARAM(164) * SPVAR(11);

    realtype ReactionFlux203 = PARAM(165) * SPVAR(36);

    realtype ReactionFlux204 = PARAM(166) * AUX_VAR_V_T * SPVAR(11);

    realtype ReactionFlux205 = PARAM(167) * SPVAR(121);

    realtype ReactionFlux206 = PARAM(31) * SPVAR(120);

    realtype ReactionFlux207 = PARAM(169) * SPVAR(64) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux208 = PARAM(171) * PARAM(170);

    realtype ReactionFlux209 = PARAM(172) * SPVAR(122);

    realtype ReactionFlux210 = PARAM(173) * AUX_VAR_H_mAPC * AUX_VAR_H_P12 * SPVAR(122);

    realtype ReactionFlux211 = PARAM(174) / AUX_VAR_N_aT * SPVAR(123);

    realtype ReactionFlux212 = PARAM(174) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(123);

    realtype ReactionFlux213 = PARAM(175) * SPVAR(12);

    realtype ReactionFlux214 = PARAM(175) * SPVAR(37);

    realtype ReactionFlux215 = PARAM(175) * SPVAR(65);

    realtype ReactionFlux216 = PARAM(175) * SPVAR(124);

    realtype ReactionFlux217 = PARAM(37) * SPVAR(65) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux218 = PARAM(180) * SPVAR(65) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux219 = PARAM(176) * SPVAR(12);

    realtype ReactionFlux220 = PARAM(177) * SPVAR(37);

    realtype ReactionFlux221 = PARAM(178) * AUX_VAR_V_T * SPVAR(12);

    realtype ReactionFlux222 = PARAM(179) * SPVAR(124);

    realtype ReactionFlux223 = PARAM(31) * SPVAR(123);

    realtype ReactionFlux224 = PARAM(181) * SPVAR(65) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux225 = PARAM(183) * PARAM(182);

    realtype ReactionFlux226 = PARAM(184) * SPVAR(125);

    realtype ReactionFlux227 = PARAM(185) * AUX_VAR_H_mAPC * AUX_VAR_H_P13 * SPVAR(125);

    realtype ReactionFlux228 = PARAM(186) / AUX_VAR_N_aT * SPVAR(126);

    realtype ReactionFlux229 = PARAM(186) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(126);

    realtype ReactionFlux230 = PARAM(187) * SPVAR(13);

    realtype ReactionFlux231 = PARAM(187) * SPVAR(38);

    realtype ReactionFlux232 = PARAM(187) * SPVAR(66);

    realtype ReactionFlux233 = PARAM(187) * SPVAR(127);

    realtype ReactionFlux234 = PARAM(37) * SPVAR(66) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux235 = PARAM(192) * SPVAR(66) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux236 = PARAM(188) * SPVAR(13);

    realtype ReactionFlux237 = PARAM(189) * SPVAR(38);

    realtype ReactionFlux238 = PARAM(190) * AUX_VAR_V_T * SPVAR(13);

    realtype ReactionFlux239 = PARAM(191) * SPVAR(127);

    realtype ReactionFlux240 = PARAM(31) * SPVAR(126);

    realtype ReactionFlux241 = PARAM(193) * SPVAR(66) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux242 = PARAM(195) * PARAM(194);

    realtype ReactionFlux243 = PARAM(196) * SPVAR(128);

    realtype ReactionFlux244 = PARAM(197) * AUX_VAR_H_mAPC * AUX_VAR_H_P14 * SPVAR(128);

    realtype ReactionFlux245 = PARAM(198) / AUX_VAR_N_aT * SPVAR(129);

    realtype ReactionFlux246 = PARAM(198) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(129);

    realtype ReactionFlux247 = PARAM(199) * SPVAR(14);

    realtype ReactionFlux248 = PARAM(199) * SPVAR(39);

    realtype ReactionFlux249 = PARAM(199) * SPVAR(67);

    realtype ReactionFlux250 = PARAM(199) * SPVAR(130);

    realtype ReactionFlux251 = PARAM(37) * SPVAR(67) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux252 = PARAM(204) * SPVAR(67) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux253 = PARAM(200) * SPVAR(14);

    realtype ReactionFlux254 = PARAM(201) * SPVAR(39);

    realtype ReactionFlux255 = PARAM(202) * AUX_VAR_V_T * SPVAR(14);

    realtype ReactionFlux256 = PARAM(203) * SPVAR(130);

    realtype ReactionFlux257 = PARAM(31) * SPVAR(129);

    realtype ReactionFlux258 = PARAM(205) * SPVAR(67) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux259 = PARAM(207) * PARAM(206);

    realtype ReactionFlux260 = PARAM(208) * SPVAR(131);

    realtype ReactionFlux261 = PARAM(209) * AUX_VAR_H_mAPC * AUX_VAR_H_P15 * SPVAR(131);

    realtype ReactionFlux262 = PARAM(210) / AUX_VAR_N_aT * SPVAR(132);

    realtype ReactionFlux263 = PARAM(210) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(132);

    realtype ReactionFlux264 = PARAM(211) * SPVAR(15);

    realtype ReactionFlux265 = PARAM(211) * SPVAR(40);

    realtype ReactionFlux266 = PARAM(211) * SPVAR(68);

    realtype ReactionFlux267 = PARAM(211) * SPVAR(133);

    realtype ReactionFlux268 = PARAM(37) * SPVAR(68) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux269 = PARAM(216) * SPVAR(68) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux270 = PARAM(212) * SPVAR(15);

    realtype ReactionFlux271 = PARAM(213) * SPVAR(40);

    realtype ReactionFlux272 = PARAM(214) * AUX_VAR_V_T * SPVAR(15);

    realtype ReactionFlux273 = PARAM(215) * SPVAR(133);

    realtype ReactionFlux274 = PARAM(31) * SPVAR(132);

    realtype ReactionFlux275 = PARAM(217) * SPVAR(68) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux276 = PARAM(219) * PARAM(218);

    realtype ReactionFlux277 = PARAM(220) * SPVAR(134);

    realtype ReactionFlux278 = PARAM(221) * AUX_VAR_H_mAPC * AUX_VAR_H_P16 * SPVAR(134);

    realtype ReactionFlux279 = PARAM(222) / AUX_VAR_N_aT * SPVAR(135);

    realtype ReactionFlux280 = PARAM(222) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(135);

    realtype ReactionFlux281 = PARAM(223) * SPVAR(16);

    realtype ReactionFlux282 = PARAM(223) * SPVAR(41);

    realtype ReactionFlux283 = PARAM(223) * SPVAR(69);

    realtype ReactionFlux284 = PARAM(223) * SPVAR(136);

    realtype ReactionFlux285 = PARAM(37) * SPVAR(69) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux286 = PARAM(228) * SPVAR(69) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux287 = PARAM(224) * SPVAR(16);

    realtype ReactionFlux288 = PARAM(225) * SPVAR(41);

    realtype ReactionFlux289 = PARAM(226) * AUX_VAR_V_T * SPVAR(16);

    realtype ReactionFlux290 = PARAM(227) * SPVAR(136);

    realtype ReactionFlux291 = PARAM(31) * SPVAR(135);

    realtype ReactionFlux292 = PARAM(229) * SPVAR(69) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux293 = PARAM(231) * PARAM(230);

    realtype ReactionFlux294 = PARAM(232) * SPVAR(137);

    realtype ReactionFlux295 = PARAM(233) * AUX_VAR_H_mAPC * AUX_VAR_H_P17 * SPVAR(137);

    realtype ReactionFlux296 = PARAM(234) / AUX_VAR_N_aT * SPVAR(138);

    realtype ReactionFlux297 = PARAM(234) / AUX_VAR_N_aT * (std::pow(2.0, AUX_VAR_N_aT) - 1.0) * SPVAR(138);

    realtype ReactionFlux298 = PARAM(235) * SPVAR(17);

    realtype ReactionFlux299 = PARAM(235) * SPVAR(42);

    realtype ReactionFlux300 = PARAM(235) * SPVAR(70);

    realtype ReactionFlux301 = PARAM(235) * SPVAR(139);

    realtype ReactionFlux302 = PARAM(37) * SPVAR(70) * AUX_VAR_Tregs_ / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9));

    realtype ReactionFlux303 = PARAM(240) * SPVAR(70) * AUX_VAR_C_total / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * AUX_VAR_H_PD1_C1;

    realtype ReactionFlux304 = PARAM(236) * SPVAR(17);

    realtype ReactionFlux305 = PARAM(237) * SPVAR(42);

    realtype ReactionFlux306 = PARAM(238) * AUX_VAR_V_T * SPVAR(17);

    realtype ReactionFlux307 = PARAM(239) * SPVAR(139);

    realtype ReactionFlux308 = PARAM(31) * SPVAR(138);

    realtype ReactionFlux309 = PARAM(241) * SPVAR(70) * SPVAR(52) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1) * (1.0 - AUX_VAR_H_MDSC_C1);

    realtype ReactionFlux310 = PARAM(244) * (PARAM(246) * AUX_VAR_V_T - SPVAR(71));

    realtype ReactionFlux311 = PARAM(244) * (PARAM(247) * PARAM(2) - SPVAR(140));

    realtype ReactionFlux312 = PARAM(242) * SPVAR(73) / (SPVAR(73) + PARAM(250)) * SPVAR(71);

    realtype ReactionFlux313 = PARAM(243) * SPVAR(72);

    realtype ReactionFlux314 = PARAM(245) * SPVAR(72);

    realtype ReactionFlux315 = PARAM(245) * SPVAR(141);

    realtype ReactionFlux316 = PARAM(248) * (PARAM(249) - SPVAR(73)) * AUX_VAR_V_T;

    realtype ReactionFlux317 = AUX_VAR_R_Tcell * PARAM(251);

    realtype ReactionFlux318 = PARAM(254) * SPVAR(200) * PARAM(4) - PARAM(253) * SPVAR(219) * PARAM(5);

    realtype ReactionFlux319 = PARAM(19) * PARAM(262) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux320 = PARAM(256) * SPVAR(142) * PARAM(2);

    realtype ReactionFlux321 = PARAM(255) * SPVAR(141) * SPVAR(142) * PARAM(2);

    realtype ReactionFlux322 = PARAM(255) * PARAM(9) * SPVAR(142) * PARAM(3);

    realtype ReactionFlux323 = PARAM(257) * SPVAR(164) * PARAM(3);

    realtype ReactionFlux324 = PARAM(258) * SPVAR(165) * PARAM(3);

    realtype ReactionFlux325 = PARAM(259) * SPVAR(165) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux326 = PARAM(260) * PARAM(259) * SPVAR(201) * PARAM(4);

    realtype ReactionFlux327 = PARAM(260) * PARAM(259) * SPVAR(220) * PARAM(5);

    realtype ReactionFlux328 = PARAM(254) * SPVAR(201) * PARAM(4);

    realtype ReactionFlux329 = PARAM(38) * PARAM(277) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux330 = PARAM(271) * SPVAR(143) * PARAM(2);

    realtype ReactionFlux331 = PARAM(270) * SPVAR(141) * SPVAR(143) * PARAM(2);

    realtype ReactionFlux332 = PARAM(270) * PARAM(9) * SPVAR(143) * PARAM(3);

    realtype ReactionFlux333 = PARAM(272) * SPVAR(166) * PARAM(3);

    realtype ReactionFlux334 = PARAM(273) * SPVAR(167) * PARAM(3);

    realtype ReactionFlux335 = PARAM(274) * SPVAR(167) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux336 = PARAM(275) * PARAM(274) * SPVAR(202) * PARAM(4);

    realtype ReactionFlux337 = PARAM(275) * PARAM(274) * SPVAR(221) * PARAM(5);

    realtype ReactionFlux338 = PARAM(254) * SPVAR(202) * PARAM(4);

    realtype ReactionFlux339 = PARAM(50) * PARAM(291) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux340 = PARAM(285) * SPVAR(144) * PARAM(2);

    realtype ReactionFlux341 = PARAM(284) * SPVAR(141) * SPVAR(144) * PARAM(2);

    realtype ReactionFlux342 = PARAM(284) * PARAM(9) * SPVAR(144) * PARAM(3);

    realtype ReactionFlux343 = PARAM(286) * SPVAR(168) * PARAM(3);

    realtype ReactionFlux344 = PARAM(287) * SPVAR(169) * PARAM(3);

    realtype ReactionFlux345 = PARAM(288) * SPVAR(169) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux346 = PARAM(289) * PARAM(288) * SPVAR(203) * PARAM(4);

    realtype ReactionFlux347 = PARAM(289) * PARAM(288) * SPVAR(222) * PARAM(5);

    realtype ReactionFlux348 = PARAM(254) * SPVAR(203) * PARAM(4);

    realtype ReactionFlux349 = PARAM(62) * PARAM(305) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux350 = PARAM(299) * SPVAR(145) * PARAM(2);

    realtype ReactionFlux351 = PARAM(298) * SPVAR(141) * SPVAR(145) * PARAM(2);

    realtype ReactionFlux352 = PARAM(298) * PARAM(9) * SPVAR(145) * PARAM(3);

    realtype ReactionFlux353 = PARAM(300) * SPVAR(170) * PARAM(3);

    realtype ReactionFlux354 = PARAM(301) * SPVAR(171) * PARAM(3);

    realtype ReactionFlux355 = PARAM(302) * SPVAR(171) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux356 = PARAM(303) * PARAM(302) * SPVAR(204) * PARAM(4);

    realtype ReactionFlux357 = PARAM(303) * PARAM(302) * SPVAR(223) * PARAM(5);

    realtype ReactionFlux358 = PARAM(254) * SPVAR(204) * PARAM(4);

    realtype ReactionFlux359 = PARAM(74) * PARAM(319) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux360 = PARAM(313) * SPVAR(146) * PARAM(2);

    realtype ReactionFlux361 = PARAM(312) * SPVAR(141) * SPVAR(146) * PARAM(2);

    realtype ReactionFlux362 = PARAM(312) * PARAM(9) * SPVAR(146) * PARAM(3);

    realtype ReactionFlux363 = PARAM(314) * SPVAR(172) * PARAM(3);

    realtype ReactionFlux364 = PARAM(315) * SPVAR(173) * PARAM(3);

    realtype ReactionFlux365 = PARAM(316) * SPVAR(173) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux366 = PARAM(317) * PARAM(316) * SPVAR(205) * PARAM(4);

    realtype ReactionFlux367 = PARAM(317) * PARAM(316) * SPVAR(224) * PARAM(5);

    realtype ReactionFlux368 = PARAM(254) * SPVAR(205) * PARAM(4);

    realtype ReactionFlux369 = PARAM(86) * PARAM(333) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux370 = PARAM(327) * SPVAR(147) * PARAM(2);

    realtype ReactionFlux371 = PARAM(326) * SPVAR(141) * SPVAR(147) * PARAM(2);

    realtype ReactionFlux372 = PARAM(326) * PARAM(9) * SPVAR(147) * PARAM(3);

    realtype ReactionFlux373 = PARAM(328) * SPVAR(174) * PARAM(3);

    realtype ReactionFlux374 = PARAM(329) * SPVAR(175) * PARAM(3);

    realtype ReactionFlux375 = PARAM(330) * SPVAR(175) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux376 = PARAM(331) * PARAM(330) * SPVAR(206) * PARAM(4);

    realtype ReactionFlux377 = PARAM(331) * PARAM(330) * SPVAR(225) * PARAM(5);

    realtype ReactionFlux378 = PARAM(254) * SPVAR(206) * PARAM(4);

    realtype ReactionFlux379 = PARAM(98) * PARAM(347) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux380 = PARAM(341) * SPVAR(148) * PARAM(2);

    realtype ReactionFlux381 = PARAM(340) * SPVAR(141) * SPVAR(148) * PARAM(2);

    realtype ReactionFlux382 = PARAM(340) * PARAM(9) * SPVAR(148) * PARAM(3);

    realtype ReactionFlux383 = PARAM(342) * SPVAR(176) * PARAM(3);

    realtype ReactionFlux384 = PARAM(343) * SPVAR(177) * PARAM(3);

    realtype ReactionFlux385 = PARAM(344) * SPVAR(177) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux386 = PARAM(345) * PARAM(344) * SPVAR(207) * PARAM(4);

    realtype ReactionFlux387 = PARAM(345) * PARAM(344) * SPVAR(226) * PARAM(5);

    realtype ReactionFlux388 = PARAM(254) * SPVAR(207) * PARAM(4);

    realtype ReactionFlux389 = PARAM(110) * PARAM(361) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux390 = PARAM(355) * SPVAR(149) * PARAM(2);

    realtype ReactionFlux391 = PARAM(354) * SPVAR(141) * SPVAR(149) * PARAM(2);

    realtype ReactionFlux392 = PARAM(354) * PARAM(9) * SPVAR(149) * PARAM(3);

    realtype ReactionFlux393 = PARAM(356) * SPVAR(178) * PARAM(3);

    realtype ReactionFlux394 = PARAM(357) * SPVAR(179) * PARAM(3);

    realtype ReactionFlux395 = PARAM(358) * SPVAR(179) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux396 = PARAM(359) * PARAM(358) * SPVAR(208) * PARAM(4);

    realtype ReactionFlux397 = PARAM(359) * PARAM(358) * SPVAR(227) * PARAM(5);

    realtype ReactionFlux398 = PARAM(254) * SPVAR(208) * PARAM(4);

    realtype ReactionFlux399 = PARAM(122) * PARAM(375) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux400 = PARAM(369) * SPVAR(150) * PARAM(2);

    realtype ReactionFlux401 = PARAM(368) * SPVAR(141) * SPVAR(150) * PARAM(2);

    realtype ReactionFlux402 = PARAM(368) * PARAM(9) * SPVAR(150) * PARAM(3);

    realtype ReactionFlux403 = PARAM(370) * SPVAR(180) * PARAM(3);

    realtype ReactionFlux404 = PARAM(371) * SPVAR(181) * PARAM(3);

    realtype ReactionFlux405 = PARAM(372) * SPVAR(181) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux406 = PARAM(373) * PARAM(372) * SPVAR(209) * PARAM(4);

    realtype ReactionFlux407 = PARAM(373) * PARAM(372) * SPVAR(228) * PARAM(5);

    realtype ReactionFlux408 = PARAM(254) * SPVAR(209) * PARAM(4);

    realtype ReactionFlux409 = PARAM(134) * PARAM(389) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux410 = PARAM(383) * SPVAR(151) * PARAM(2);

    realtype ReactionFlux411 = PARAM(382) * SPVAR(141) * SPVAR(151) * PARAM(2);

    realtype ReactionFlux412 = PARAM(382) * PARAM(9) * SPVAR(151) * PARAM(3);

    realtype ReactionFlux413 = PARAM(384) * SPVAR(182) * PARAM(3);

    realtype ReactionFlux414 = PARAM(385) * SPVAR(183) * PARAM(3);

    realtype ReactionFlux415 = PARAM(386) * SPVAR(183) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux416 = PARAM(387) * PARAM(386) * SPVAR(210) * PARAM(4);

    realtype ReactionFlux417 = PARAM(387) * PARAM(386) * SPVAR(229) * PARAM(5);

    realtype ReactionFlux418 = PARAM(254) * SPVAR(210) * PARAM(4);

    realtype ReactionFlux419 = PARAM(146) * PARAM(403) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux420 = PARAM(397) * SPVAR(152) * PARAM(2);

    realtype ReactionFlux421 = PARAM(396) * SPVAR(141) * SPVAR(152) * PARAM(2);

    realtype ReactionFlux422 = PARAM(396) * PARAM(9) * SPVAR(152) * PARAM(3);

    realtype ReactionFlux423 = PARAM(398) * SPVAR(184) * PARAM(3);

    realtype ReactionFlux424 = PARAM(399) * SPVAR(185) * PARAM(3);

    realtype ReactionFlux425 = PARAM(400) * SPVAR(185) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux426 = PARAM(401) * PARAM(400) * SPVAR(211) * PARAM(4);

    realtype ReactionFlux427 = PARAM(401) * PARAM(400) * SPVAR(230) * PARAM(5);

    realtype ReactionFlux428 = PARAM(254) * SPVAR(211) * PARAM(4);

    realtype ReactionFlux429 = PARAM(158) * PARAM(417) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux430 = PARAM(411) * SPVAR(153) * PARAM(2);

    realtype ReactionFlux431 = PARAM(410) * SPVAR(141) * SPVAR(153) * PARAM(2);

    realtype ReactionFlux432 = PARAM(410) * PARAM(9) * SPVAR(153) * PARAM(3);

    realtype ReactionFlux433 = PARAM(412) * SPVAR(186) * PARAM(3);

    realtype ReactionFlux434 = PARAM(413) * SPVAR(187) * PARAM(3);

    realtype ReactionFlux435 = PARAM(414) * SPVAR(187) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux436 = PARAM(415) * PARAM(414) * SPVAR(212) * PARAM(4);

    realtype ReactionFlux437 = PARAM(415) * PARAM(414) * SPVAR(231) * PARAM(5);

    realtype ReactionFlux438 = PARAM(254) * SPVAR(212) * PARAM(4);

    realtype ReactionFlux439 = PARAM(170) * PARAM(431) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux440 = PARAM(425) * SPVAR(154) * PARAM(2);

    realtype ReactionFlux441 = PARAM(424) * SPVAR(141) * SPVAR(154) * PARAM(2);

    realtype ReactionFlux442 = PARAM(424) * PARAM(9) * SPVAR(154) * PARAM(3);

    realtype ReactionFlux443 = PARAM(426) * SPVAR(188) * PARAM(3);

    realtype ReactionFlux444 = PARAM(427) * SPVAR(189) * PARAM(3);

    realtype ReactionFlux445 = PARAM(428) * SPVAR(189) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux446 = PARAM(429) * PARAM(428) * SPVAR(213) * PARAM(4);

    realtype ReactionFlux447 = PARAM(429) * PARAM(428) * SPVAR(232) * PARAM(5);

    realtype ReactionFlux448 = PARAM(254) * SPVAR(213) * PARAM(4);

    realtype ReactionFlux449 = PARAM(182) * PARAM(445) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux450 = PARAM(439) * SPVAR(155) * PARAM(2);

    realtype ReactionFlux451 = PARAM(438) * SPVAR(141) * SPVAR(155) * PARAM(2);

    realtype ReactionFlux452 = PARAM(438) * PARAM(9) * SPVAR(155) * PARAM(3);

    realtype ReactionFlux453 = PARAM(440) * SPVAR(190) * PARAM(3);

    realtype ReactionFlux454 = PARAM(441) * SPVAR(191) * PARAM(3);

    realtype ReactionFlux455 = PARAM(442) * SPVAR(191) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux456 = PARAM(443) * PARAM(442) * SPVAR(214) * PARAM(4);

    realtype ReactionFlux457 = PARAM(443) * PARAM(442) * SPVAR(233) * PARAM(5);

    realtype ReactionFlux458 = PARAM(254) * SPVAR(214) * PARAM(4);

    realtype ReactionFlux459 = PARAM(194) * PARAM(459) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux460 = PARAM(453) * SPVAR(156) * PARAM(2);

    realtype ReactionFlux461 = PARAM(452) * SPVAR(141) * SPVAR(156) * PARAM(2);

    realtype ReactionFlux462 = PARAM(452) * PARAM(9) * SPVAR(156) * PARAM(3);

    realtype ReactionFlux463 = PARAM(454) * SPVAR(192) * PARAM(3);

    realtype ReactionFlux464 = PARAM(455) * SPVAR(193) * PARAM(3);

    realtype ReactionFlux465 = PARAM(456) * SPVAR(193) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux466 = PARAM(457) * PARAM(456) * SPVAR(215) * PARAM(4);

    realtype ReactionFlux467 = PARAM(457) * PARAM(456) * SPVAR(234) * PARAM(5);

    realtype ReactionFlux468 = PARAM(254) * SPVAR(215) * PARAM(4);

    realtype ReactionFlux469 = PARAM(206) * PARAM(473) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux470 = PARAM(467) * SPVAR(157) * PARAM(2);

    realtype ReactionFlux471 = PARAM(466) * SPVAR(141) * SPVAR(157) * PARAM(2);

    realtype ReactionFlux472 = PARAM(466) * PARAM(9) * SPVAR(157) * PARAM(3);

    realtype ReactionFlux473 = PARAM(468) * SPVAR(194) * PARAM(3);

    realtype ReactionFlux474 = PARAM(469) * SPVAR(195) * PARAM(3);

    realtype ReactionFlux475 = PARAM(470) * SPVAR(195) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux476 = PARAM(471) * PARAM(470) * SPVAR(216) * PARAM(4);

    realtype ReactionFlux477 = PARAM(471) * PARAM(470) * SPVAR(235) * PARAM(5);

    realtype ReactionFlux478 = PARAM(254) * SPVAR(216) * PARAM(4);

    realtype ReactionFlux479 = PARAM(218) * PARAM(487) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux480 = PARAM(481) * SPVAR(158) * PARAM(2);

    realtype ReactionFlux481 = PARAM(480) * SPVAR(141) * SPVAR(158) * PARAM(2);

    realtype ReactionFlux482 = PARAM(480) * PARAM(9) * SPVAR(158) * PARAM(3);

    realtype ReactionFlux483 = PARAM(482) * SPVAR(196) * PARAM(3);

    realtype ReactionFlux484 = PARAM(483) * SPVAR(197) * PARAM(3);

    realtype ReactionFlux485 = PARAM(484) * SPVAR(197) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux486 = PARAM(485) * PARAM(484) * SPVAR(217) * PARAM(4);

    realtype ReactionFlux487 = PARAM(485) * PARAM(484) * SPVAR(236) * PARAM(5);

    realtype ReactionFlux488 = PARAM(254) * SPVAR(217) * PARAM(4);

    realtype ReactionFlux489 = PARAM(230) * PARAM(501) * (PARAM(16) + PARAM(17) + (PARAM(49) * SPVAR(54) + PARAM(61) * SPVAR(55) + PARAM(73) * SPVAR(56) + PARAM(85) * SPVAR(57) + PARAM(97) * SPVAR(58) + PARAM(109) * SPVAR(59) + PARAM(121) * SPVAR(60) + PARAM(133) * SPVAR(61) + PARAM(145) * SPVAR(62) + PARAM(157) * SPVAR(63) + PARAM(169) * SPVAR(64) + PARAM(181) * SPVAR(65) + PARAM(193) * SPVAR(66) + PARAM(205) * SPVAR(67) + PARAM(217) * SPVAR(68) + PARAM(229) * SPVAR(69) + PARAM(241) * SPVAR(70)) / (AUX_VAR_C_total + AUX_VAR_T_total + PARAM(9)) * (1.0 - AUX_VAR_H_PD1_C1)) * SPVAR(52) * AUX_VAR_V_T;

    realtype ReactionFlux490 = PARAM(495) * SPVAR(159) * PARAM(2);

    realtype ReactionFlux491 = PARAM(494) * SPVAR(141) * SPVAR(159) * PARAM(2);

    realtype ReactionFlux492 = PARAM(494) * PARAM(9) * SPVAR(159) * PARAM(3);

    realtype ReactionFlux493 = PARAM(496) * SPVAR(198) * PARAM(3);

    realtype ReactionFlux494 = PARAM(497) * SPVAR(199) * PARAM(3);

    realtype ReactionFlux495 = PARAM(498) * SPVAR(199) * SPVAR(200) * PARAM(4);

    realtype ReactionFlux496 = PARAM(499) * PARAM(498) * SPVAR(218) * PARAM(4);

    realtype ReactionFlux497 = PARAM(499) * PARAM(498) * SPVAR(237) * PARAM(5);

    realtype ReactionFlux498 = PARAM(254) * SPVAR(218) * PARAM(4);

    realtype ReactionFlux499 = PARAM(509) * (SPVAR(18) / PARAM(514) - SPVAR(43) / PARAM(515));

    realtype ReactionFlux500 = PARAM(510) * (SPVAR(18) / PARAM(514) - SPVAR(74) / PARAM(516));

    realtype ReactionFlux501 = PARAM(511) * (SPVAR(18) / PARAM(514) - SPVAR(160) / PARAM(517));

    realtype ReactionFlux502 = PARAM(512) * SPVAR(74) / PARAM(516) * AUX_VAR_V_T;

    realtype ReactionFlux503 = PARAM(512) * SPVAR(160) / PARAM(517) * PARAM(2);

    realtype ReactionFlux504 = PARAM(513) * SPVAR(18);

    realtype ReactionFlux505 = PARAM(518) * (SPVAR(19) / PARAM(523) - SPVAR(44) / PARAM(524));

    realtype ReactionFlux506 = PARAM(519) * (SPVAR(19) / PARAM(523) - SPVAR(75) / PARAM(525));

    realtype ReactionFlux507 = PARAM(520) * (SPVAR(19) / PARAM(523) - SPVAR(161) / PARAM(526));

    realtype ReactionFlux508 = PARAM(521) * SPVAR(75) / PARAM(525) * AUX_VAR_V_T;

    realtype ReactionFlux509 = PARAM(521) * SPVAR(161) / PARAM(526) * PARAM(2);

    realtype ReactionFlux510 = PARAM(522) * SPVAR(19);

    realtype ReactionFlux511 = PARAM(527) * (SPVAR(20) / PARAM(532) - SPVAR(45) / PARAM(533));

    realtype ReactionFlux512 = PARAM(528) * (SPVAR(20) / PARAM(532) - SPVAR(76) / PARAM(534));

    realtype ReactionFlux513 = PARAM(529) * (SPVAR(20) / PARAM(532) - SPVAR(162) / PARAM(535));

    realtype ReactionFlux514 = PARAM(530) * SPVAR(76) / PARAM(534) * AUX_VAR_V_T;

    realtype ReactionFlux515 = PARAM(530) * SPVAR(162) / PARAM(535) * PARAM(2);

    realtype ReactionFlux516 = PARAM(531) * SPVAR(20);

    realtype ReactionFlux517 = (PARAM(508) * SPVAR(240) * SPVAR(241) - PARAM(545) * SPVAR(238)) * PARAM(6);

    realtype ReactionFlux518 = (PARAM(536) * SPVAR(240) * SPVAR(242) - PARAM(546) * SPVAR(239)) * PARAM(6);

    realtype ReactionFlux519 = (2.0 * PARAM(537) * (SPVAR(240) * SPVAR(74) / PARAM(516)) - PARAM(547) * SPVAR(243)) * PARAM(6);

    realtype ReactionFlux520 = (PARAM(555) * PARAM(537) * SPVAR(240) * SPVAR(243) - 2.0 * PARAM(547) * SPVAR(244)) * PARAM(6);

    realtype ReactionFlux521 = (2.0 * PARAM(538) * (SPVAR(241) * SPVAR(75) / PARAM(525)) - PARAM(548) * SPVAR(245)) * PARAM(6);

    realtype ReactionFlux522 = (PARAM(556) * PARAM(538) * SPVAR(241) * SPVAR(245) - 2.0 * PARAM(548) * SPVAR(246)) * PARAM(6);

    realtype ReactionFlux523 = (2.0 * PARAM(539) * SPVAR(261) * SPVAR(263) - PARAM(549) * SPVAR(250)) * PARAM(6);

    realtype ReactionFlux524 = (PARAM(539) * SPVAR(261) * SPVAR(250) - 2.0 * PARAM(549) * SPVAR(251)) * PARAM(6);

    realtype ReactionFlux525 = (PARAM(540) * SPVAR(261) * SPVAR(264) - PARAM(550) * SPVAR(252)) * PARAM(6);

    realtype ReactionFlux526 = (4.0 * PARAM(541) * SPVAR(262) * SPVAR(263) - PARAM(551) * SPVAR(253)) * PARAM(6);

    realtype ReactionFlux527 = (PARAM(541) * SPVAR(262) * SPVAR(253) - 2.0 * PARAM(551) * SPVAR(255)) * PARAM(6);

    realtype ReactionFlux528 = (PARAM(541) * SPVAR(263) * SPVAR(255) - PARAM(551) * SPVAR(256)) * PARAM(6);

    realtype ReactionFlux529 = (PARAM(541) * SPVAR(253) * SPVAR(263) - 2.0 * PARAM(551) * SPVAR(254)) * PARAM(6);

    realtype ReactionFlux530 = (PARAM(541) * SPVAR(262) * SPVAR(254) - PARAM(551) * SPVAR(256)) * PARAM(6);

    realtype ReactionFlux531 = (2.0 * PARAM(542) * SPVAR(262) * SPVAR(264) - PARAM(552) * SPVAR(257)) * PARAM(6);

    realtype ReactionFlux532 = (PARAM(542) * SPVAR(257) * SPVAR(264) - 2.0 * PARAM(552) * SPVAR(258)) * PARAM(6);

    realtype ReactionFlux533 = (4.0 * PARAM(544) * (SPVAR(262) * SPVAR(76) / PARAM(534)) - PARAM(554) * SPVAR(265)) * PARAM(6);

    realtype ReactionFlux534 = (PARAM(557) * PARAM(544) * SPVAR(262) * SPVAR(265) - 2.0 * PARAM(554) * SPVAR(266)) * PARAM(6);

    realtype ReactionFlux535 = (2.0 * PARAM(543) * SPVAR(263) * SPVAR(247) - PARAM(553) * SPVAR(259)) * PARAM(6);

    realtype ReactionFlux536 = (PARAM(543) * SPVAR(259) * SPVAR(247) - 2.0 * PARAM(553) * SPVAR(260)) * PARAM(6);

    realtype ReactionFlux537 = (2.0 * PARAM(538) * (SPVAR(247) * SPVAR(75) / PARAM(525)) - PARAM(548) * SPVAR(248)) * PARAM(6);

    realtype ReactionFlux538 = (PARAM(556) * PARAM(538) * SPVAR(247) * SPVAR(248) - 2.0 * PARAM(548) * SPVAR(249)) * PARAM(6);

    realtype ReactionFlux539 = (PARAM(508) * SPVAR(269) * SPVAR(270) - PARAM(545) * SPVAR(267)) * PARAM(7);

    realtype ReactionFlux540 = (PARAM(536) * SPVAR(269) * SPVAR(271) - PARAM(546) * SPVAR(268)) * PARAM(7);

    realtype ReactionFlux541 = (2.0 * PARAM(537) * (SPVAR(269) * SPVAR(160) / PARAM(517)) - PARAM(547) * SPVAR(272)) * PARAM(7);

    realtype ReactionFlux542 = (PARAM(555) * PARAM(537) * SPVAR(269) * SPVAR(272) - 2.0 * PARAM(547) * SPVAR(273)) * PARAM(7);

    realtype ReactionFlux543 = (2.0 * PARAM(538) * (SPVAR(270) * SPVAR(161) / PARAM(526)) - PARAM(548) * SPVAR(274)) * PARAM(7);

    realtype ReactionFlux544 = (PARAM(556) * PARAM(538) * SPVAR(270) * SPVAR(274) - 2.0 * PARAM(548) * SPVAR(275)) * PARAM(7);

    realtype ReactionFlux545 = (2.0 * PARAM(539) * SPVAR(290) * SPVAR(292) - PARAM(549) * SPVAR(279)) * PARAM(7);

    realtype ReactionFlux546 = (PARAM(539) * SPVAR(290) * SPVAR(279) - 2.0 * PARAM(549) * SPVAR(280)) * PARAM(7);

    realtype ReactionFlux547 = (PARAM(540) * SPVAR(290) * SPVAR(293) - PARAM(550) * SPVAR(281)) * PARAM(7);

    realtype ReactionFlux548 = (4.0 * PARAM(541) * SPVAR(291) * SPVAR(292) - PARAM(551) * SPVAR(282)) * PARAM(7);

    realtype ReactionFlux549 = (PARAM(541) * SPVAR(291) * SPVAR(282) - 2.0 * PARAM(551) * SPVAR(284)) * PARAM(7);

    realtype ReactionFlux550 = (PARAM(541) * SPVAR(292) * SPVAR(284) - PARAM(551) * SPVAR(285)) * PARAM(7);

    realtype ReactionFlux551 = (PARAM(541) * SPVAR(282) * SPVAR(292) - 2.0 * PARAM(551) * SPVAR(283)) * PARAM(7);

    realtype ReactionFlux552 = (PARAM(541) * SPVAR(291) * SPVAR(283) - PARAM(551) * SPVAR(285)) * PARAM(7);

    realtype ReactionFlux553 = (2.0 * PARAM(542) * SPVAR(291) * SPVAR(293) - PARAM(552) * SPVAR(286)) * PARAM(7);

    realtype ReactionFlux554 = (PARAM(542) * SPVAR(286) * SPVAR(293) - 2.0 * PARAM(552) * SPVAR(287)) * PARAM(7);

    realtype ReactionFlux555 = (4.0 * PARAM(544) * (SPVAR(291) * SPVAR(162) / PARAM(535)) - PARAM(554) * SPVAR(294)) * PARAM(7);

    realtype ReactionFlux556 = (PARAM(557) * PARAM(544) * SPVAR(291) * SPVAR(294) - 2.0 * PARAM(554) * SPVAR(295)) * PARAM(7);

    realtype ReactionFlux557 = (2.0 * PARAM(543) * SPVAR(292) * SPVAR(276) - PARAM(553) * SPVAR(288)) * PARAM(7);

    realtype ReactionFlux558 = (PARAM(543) * SPVAR(288) * SPVAR(276) - 2.0 * PARAM(553) * SPVAR(289)) * PARAM(7);

    realtype ReactionFlux559 = (2.0 * PARAM(538) * (SPVAR(276) * SPVAR(161) / PARAM(526)) - PARAM(548) * SPVAR(277)) * PARAM(7);

    realtype ReactionFlux560 = (PARAM(556) * PARAM(538) * SPVAR(276) * SPVAR(277) - 2.0 * PARAM(548) * SPVAR(278)) * PARAM(7);

    realtype ReactionFlux561 = PARAM(544) * (SPVAR(77) * SPVAR(76) / PARAM(532)) - PARAM(554) * SPVAR(78);

    realtype ReactionFlux562 = PARAM(557) * PARAM(544) * SPVAR(77) * SPVAR(78) / PARAM(264) - PARAM(554) * SPVAR(79);

    realtype ReactionFlux563 = PARAM(544) * (SPVAR(46) * SPVAR(45) / PARAM(533)) - PARAM(554) * SPVAR(47);

    realtype ReactionFlux564 = PARAM(557) * PARAM(544) * SPVAR(46) * SPVAR(47) / PARAM(264) - PARAM(554) * SPVAR(48);

    realtype ReactionFlux565 = PARAM(577) * SPVAR(53) * AUX_VAR_H_Treg_T;

    realtype ReactionFlux566 = PARAM(577) * SPVAR(25) * AUX_VAR_H_Treg_P;

    realtype ReactionFlux567 = PARAM(578) * (PARAM(593) * AUX_VAR_V_T - SPVAR(80)) * (SPVAR(81) / (PARAM(591) + SPVAR(81)));

    realtype ReactionFlux568 = PARAM(596) * (PARAM(593) * AUX_VAR_V_T - SPVAR(80));

    realtype ReactionFlux569 = PARAM(579) * SPVAR(80);

    realtype ReactionFlux570 = PARAM(581) * SPVAR(81) * AUX_VAR_V_T;

    realtype ReactionFlux571 = PARAM(582) * SPVAR(82) * AUX_VAR_V_T;

    realtype ReactionFlux572 = PARAM(583) * SPVAR(83) * AUX_VAR_V_T;

    realtype ReactionFlux573 = PARAM(584) * SPVAR(52) * (1.0 - SPVAR(84) / (SPVAR(84) + PARAM(595)));

    realtype ReactionFlux574 = PARAM(585) * SPVAR(80) * (1.0 - SPVAR(84) / (SPVAR(84) + PARAM(587)));

    realtype ReactionFlux575 = PARAM(586) * SPVAR(80) * (1.0 - SPVAR(84) / (SPVAR(84) + PARAM(597)));

    realtype ReactionFlux576 = PARAM(588) * SPVAR(53) * (1.0 - SPVAR(53) / (PARAM(594) * AUX_VAR_V_T)) * SPVAR(83) / (PARAM(592) + SPVAR(83)) * (1.0 - SPVAR(84) / (SPVAR(84) + PARAM(597)));

    realtype ReactionFlux577 = PARAM(604) * SPVAR(24) * PARAM(0);

    realtype ReactionFlux578 = PARAM(598) * SPVAR(22) * PARAM(0);

    realtype ReactionFlux579 = PARAM(599) * SPVAR(23) * PARAM(0);

    realtype ReactionFlux580 = PARAM(600) * SPVAR(21) / (SPVAR(21) + PARAM(601)) * PARAM(0);

    realtype ReactionFlux581 = PARAM(605) * (SPVAR(21) / PARAM(610) - SPVAR(49) / PARAM(611));

    realtype ReactionFlux582 = PARAM(606) * (SPVAR(21) / PARAM(610) - SPVAR(84) / PARAM(612));

    realtype ReactionFlux583 = PARAM(607) * (SPVAR(21) / PARAM(610) - SPVAR(163) / PARAM(613));

    realtype ReactionFlux584 = PARAM(608) * SPVAR(84) / PARAM(612) * AUX_VAR_V_T;

    realtype ReactionFlux585 = PARAM(608) * SPVAR(163) / PARAM(613) * PARAM(2);

    realtype ReactionFlux586 = PARAM(609) * SPVAR(21) * PARAM(0);

    //dydt:

    //d(V_C.T0)/dt
    NV_DATA_S(ydot)[0] =  - ReactionFlux10 - ReactionFlux14 + ReactionFlux15 - ReactionFlux16 + ReactionFlux17;

    //d(V_C.T1)/dt
    NV_DATA_S(ydot)[1] =  - ReactionFlux26 - ReactionFlux32 + ReactionFlux33 - ReactionFlux34 + ReactionFlux35;

    //d(V_C.T2)/dt
    NV_DATA_S(ydot)[2] =  - ReactionFlux43 - ReactionFlux49 + ReactionFlux50 - ReactionFlux51 + ReactionFlux52;

    //d(V_C.T3)/dt
    NV_DATA_S(ydot)[3] =  - ReactionFlux60 - ReactionFlux66 + ReactionFlux67 - ReactionFlux68 + ReactionFlux69;

    //d(V_C.T4)/dt
    NV_DATA_S(ydot)[4] =  - ReactionFlux77 - ReactionFlux83 + ReactionFlux84 - ReactionFlux85 + ReactionFlux86;

    //d(V_C.T5)/dt
    NV_DATA_S(ydot)[5] =  - ReactionFlux94 - ReactionFlux100 + ReactionFlux101 - ReactionFlux102 + ReactionFlux103;

    //d(V_C.T6)/dt
    NV_DATA_S(ydot)[6] =  - ReactionFlux111 - ReactionFlux117 + ReactionFlux118 - ReactionFlux119 + ReactionFlux120;

    //d(V_C.T7)/dt
    NV_DATA_S(ydot)[7] =  - ReactionFlux128 - ReactionFlux134 + ReactionFlux135 - ReactionFlux136 + ReactionFlux137;

    //d(V_C.T8)/dt
    NV_DATA_S(ydot)[8] =  - ReactionFlux145 - ReactionFlux151 + ReactionFlux152 - ReactionFlux153 + ReactionFlux154;

    //d(V_C.T9)/dt
    NV_DATA_S(ydot)[9] =  - ReactionFlux162 - ReactionFlux168 + ReactionFlux169 - ReactionFlux170 + ReactionFlux171;

    //d(V_C.T10)/dt
    NV_DATA_S(ydot)[10] =  - ReactionFlux179 - ReactionFlux185 + ReactionFlux186 - ReactionFlux187 + ReactionFlux188;

    //d(V_C.T11)/dt
    NV_DATA_S(ydot)[11] =  - ReactionFlux196 - ReactionFlux202 + ReactionFlux203 - ReactionFlux204 + ReactionFlux205;

    //d(V_C.T12)/dt
    NV_DATA_S(ydot)[12] =  - ReactionFlux213 - ReactionFlux219 + ReactionFlux220 - ReactionFlux221 + ReactionFlux222;

    //d(V_C.T13)/dt
    NV_DATA_S(ydot)[13] =  - ReactionFlux230 - ReactionFlux236 + ReactionFlux237 - ReactionFlux238 + ReactionFlux239;

    //d(V_C.T14)/dt
    NV_DATA_S(ydot)[14] =  - ReactionFlux247 - ReactionFlux253 + ReactionFlux254 - ReactionFlux255 + ReactionFlux256;

    //d(V_C.T15)/dt
    NV_DATA_S(ydot)[15] =  - ReactionFlux264 - ReactionFlux270 + ReactionFlux271 - ReactionFlux272 + ReactionFlux273;

    //d(V_C.T16)/dt
    NV_DATA_S(ydot)[16] =  - ReactionFlux281 - ReactionFlux287 + ReactionFlux288 - ReactionFlux289 + ReactionFlux290;

    //d(V_C.T17)/dt
    NV_DATA_S(ydot)[17] =  - ReactionFlux298 - ReactionFlux304 + ReactionFlux305 - ReactionFlux306 + ReactionFlux307;

    //d(V_C.nivo)/dt
    NV_DATA_S(ydot)[18] = 1/PARAM(0)*( - ReactionFlux499 - ReactionFlux500 - ReactionFlux501 + ReactionFlux503 - ReactionFlux504);

    //d(V_C.durv)/dt
    NV_DATA_S(ydot)[19] = 1/PARAM(0)*( - ReactionFlux505 - ReactionFlux506 - ReactionFlux507 + ReactionFlux509 - ReactionFlux510);

    //d(V_C.ipi)/dt
    NV_DATA_S(ydot)[20] = 1/PARAM(0)*( - ReactionFlux511 - ReactionFlux512 - ReactionFlux513 + ReactionFlux515 - ReactionFlux516);

    //d(V_C.ENT)/dt
    NV_DATA_S(ydot)[21] = 1/PARAM(0)*(ReactionFlux578 + ReactionFlux579 - ReactionFlux580 - ReactionFlux581 - ReactionFlux582 - ReactionFlux583 + ReactionFlux585 - ReactionFlux586);

    //d(V_C.ENT_Buccal)/dt
    NV_DATA_S(ydot)[22] = 1/PARAM(0)*( - ReactionFlux578);

    //d(V_C.ENT_GI)/dt
    NV_DATA_S(ydot)[23] = 1/PARAM(0)*(ReactionFlux577 - ReactionFlux579);

    //d(V_C.Dose2)/dt
    NV_DATA_S(ydot)[24] = 1/PARAM(0)*( - ReactionFlux577);

    //d(V_P.T0)/dt
    NV_DATA_S(ydot)[25] =  - ReactionFlux11 + ReactionFlux14 - ReactionFlux15 - ReactionFlux566;

    //d(V_P.T1)/dt
    NV_DATA_S(ydot)[26] =  - ReactionFlux27 + ReactionFlux32 - ReactionFlux33;

    //d(V_P.T2)/dt
    NV_DATA_S(ydot)[27] =  - ReactionFlux44 + ReactionFlux49 - ReactionFlux50;

    //d(V_P.T3)/dt
    NV_DATA_S(ydot)[28] =  - ReactionFlux61 + ReactionFlux66 - ReactionFlux67;

    //d(V_P.T4)/dt
    NV_DATA_S(ydot)[29] =  - ReactionFlux78 + ReactionFlux83 - ReactionFlux84;

    //d(V_P.T5)/dt
    NV_DATA_S(ydot)[30] =  - ReactionFlux95 + ReactionFlux100 - ReactionFlux101;

    //d(V_P.T6)/dt
    NV_DATA_S(ydot)[31] =  - ReactionFlux112 + ReactionFlux117 - ReactionFlux118;

    //d(V_P.T7)/dt
    NV_DATA_S(ydot)[32] =  - ReactionFlux129 + ReactionFlux134 - ReactionFlux135;

    //d(V_P.T8)/dt
    NV_DATA_S(ydot)[33] =  - ReactionFlux146 + ReactionFlux151 - ReactionFlux152;

    //d(V_P.T9)/dt
    NV_DATA_S(ydot)[34] =  - ReactionFlux163 + ReactionFlux168 - ReactionFlux169;

    //d(V_P.T10)/dt
    NV_DATA_S(ydot)[35] =  - ReactionFlux180 + ReactionFlux185 - ReactionFlux186;

    //d(V_P.T11)/dt
    NV_DATA_S(ydot)[36] =  - ReactionFlux197 + ReactionFlux202 - ReactionFlux203;

    //d(V_P.T12)/dt
    NV_DATA_S(ydot)[37] =  - ReactionFlux214 + ReactionFlux219 - ReactionFlux220;

    //d(V_P.T13)/dt
    NV_DATA_S(ydot)[38] =  - ReactionFlux231 + ReactionFlux236 - ReactionFlux237;

    //d(V_P.T14)/dt
    NV_DATA_S(ydot)[39] =  - ReactionFlux248 + ReactionFlux253 - ReactionFlux254;

    //d(V_P.T15)/dt
    NV_DATA_S(ydot)[40] =  - ReactionFlux265 + ReactionFlux270 - ReactionFlux271;

    //d(V_P.T16)/dt
    NV_DATA_S(ydot)[41] =  - ReactionFlux282 + ReactionFlux287 - ReactionFlux288;

    //d(V_P.T17)/dt
    NV_DATA_S(ydot)[42] =  - ReactionFlux299 + ReactionFlux304 - ReactionFlux305;

    //d(V_P.nivo)/dt
    NV_DATA_S(ydot)[43] = 1/PARAM(1)*(ReactionFlux499);

    //d(V_P.durv)/dt
    NV_DATA_S(ydot)[44] = 1/PARAM(1)*(ReactionFlux505);

    //d(V_P.ipi)/dt
    NV_DATA_S(ydot)[45] = 1/PARAM(1)*(ReactionFlux511);

    //d(V_P.Treg_CTLA4)/dt
    NV_DATA_S(ydot)[46] =  - ReactionFlux563 - ReactionFlux564;

    //d(V_P.Treg_CTLA4_ipi)/dt
    NV_DATA_S(ydot)[47] = ReactionFlux563 - ReactionFlux564;

    //d(V_P.Treg_CTLA4_ipi_CTLA4)/dt
    NV_DATA_S(ydot)[48] = ReactionFlux564;

    //d(V_P.ENT)/dt
    NV_DATA_S(ydot)[49] = 1/PARAM(1)*(ReactionFlux581);

    //d(V_T.C_x)/dt
    NV_DATA_S(ydot)[50] =  - ReactionFlux1 + ReactionFlux4 + ReactionFlux37 + ReactionFlux54 + ReactionFlux71 + ReactionFlux88 + ReactionFlux105 + ReactionFlux122 + ReactionFlux139 + ReactionFlux156 + ReactionFlux173 + ReactionFlux190 + ReactionFlux207 + ReactionFlux224 + ReactionFlux241 + ReactionFlux258 + ReactionFlux275 + ReactionFlux292 + ReactionFlux309;

    //d(V_T.T_exh)/dt
    NV_DATA_S(ydot)[51] =  - ReactionFlux2 + ReactionFlux12 + ReactionFlux28 + ReactionFlux30 + ReactionFlux31 + ReactionFlux45 + ReactionFlux47 + ReactionFlux48 + ReactionFlux62 + ReactionFlux64 + ReactionFlux65 + ReactionFlux79 + ReactionFlux81 + ReactionFlux82 + ReactionFlux96 + ReactionFlux98 + ReactionFlux99 + ReactionFlux113 + ReactionFlux115 + ReactionFlux116 + ReactionFlux130 + ReactionFlux132 + ReactionFlux133 + ReactionFlux147 + ReactionFlux149 + ReactionFlux150 + ReactionFlux164 + ReactionFlux166 + ReactionFlux167 + ReactionFlux181 + ReactionFlux183 + ReactionFlux184 + ReactionFlux198 + ReactionFlux200 + ReactionFlux201 + ReactionFlux215 + ReactionFlux217 + ReactionFlux218 + ReactionFlux232 + ReactionFlux234 + ReactionFlux235 + ReactionFlux249 + ReactionFlux251 + ReactionFlux252 + ReactionFlux266 + ReactionFlux268 + ReactionFlux269 + ReactionFlux283 + ReactionFlux285 + ReactionFlux286 + ReactionFlux300 + ReactionFlux302 + ReactionFlux303;

    //d(V_T.C1)/dt
	NV_DATA_S(ydot)[52] = ReactionFlux3 - ReactionFlux4 - ReactionFlux37 - ReactionFlux54 - ReactionFlux71 - ReactionFlux88 - ReactionFlux105 - ReactionFlux122 - ReactionFlux139 - ReactionFlux156 - ReactionFlux173 - ReactionFlux190 - ReactionFlux207 - ReactionFlux224 - ReactionFlux241 - ReactionFlux258 - ReactionFlux275 - ReactionFlux292 - ReactionFlux309;

    //d(V_T.T0)/dt
    NV_DATA_S(ydot)[53] =  - ReactionFlux12 + ReactionFlux16 - ReactionFlux565 + ReactionFlux576;

    //d(V_T.T1)/dt
    NV_DATA_S(ydot)[54] =  - ReactionFlux28 - ReactionFlux30 - ReactionFlux31 + ReactionFlux34;

    //d(V_T.T2)/dt
    NV_DATA_S(ydot)[55] =  - ReactionFlux45 - ReactionFlux47 - ReactionFlux48 + ReactionFlux51;

    //d(V_T.T3)/dt
    NV_DATA_S(ydot)[56] =  - ReactionFlux62 - ReactionFlux64 - ReactionFlux65 + ReactionFlux68;

    //d(V_T.T4)/dt
    NV_DATA_S(ydot)[57] =  - ReactionFlux79 - ReactionFlux81 - ReactionFlux82 + ReactionFlux85;

    //d(V_T.T5)/dt
    NV_DATA_S(ydot)[58] =  - ReactionFlux96 - ReactionFlux98 - ReactionFlux99 + ReactionFlux102;

    //d(V_T.T6)/dt
    NV_DATA_S(ydot)[59] =  - ReactionFlux113 - ReactionFlux115 - ReactionFlux116 + ReactionFlux119;

    //d(V_T.T7)/dt
    NV_DATA_S(ydot)[60] =  - ReactionFlux130 - ReactionFlux132 - ReactionFlux133 + ReactionFlux136;

    //d(V_T.T8)/dt
    NV_DATA_S(ydot)[61] =  - ReactionFlux147 - ReactionFlux149 - ReactionFlux150 + ReactionFlux153;

    //d(V_T.T9)/dt
    NV_DATA_S(ydot)[62] =  - ReactionFlux164 - ReactionFlux166 - ReactionFlux167 + ReactionFlux170;

    //d(V_T.T10)/dt
    NV_DATA_S(ydot)[63] =  - ReactionFlux181 - ReactionFlux183 - ReactionFlux184 + ReactionFlux187;

    //d(V_T.T11)/dt
    NV_DATA_S(ydot)[64] =  - ReactionFlux198 - ReactionFlux200 - ReactionFlux201 + ReactionFlux204;

    //d(V_T.T12)/dt
    NV_DATA_S(ydot)[65] =  - ReactionFlux215 - ReactionFlux217 - ReactionFlux218 + ReactionFlux221;

    //d(V_T.T13)/dt
    NV_DATA_S(ydot)[66] =  - ReactionFlux232 - ReactionFlux234 - ReactionFlux235 + ReactionFlux238;

    //d(V_T.T14)/dt
    NV_DATA_S(ydot)[67] =  - ReactionFlux249 - ReactionFlux251 - ReactionFlux252 + ReactionFlux255;

    //d(V_T.T15)/dt
    NV_DATA_S(ydot)[68] =  - ReactionFlux266 - ReactionFlux268 - ReactionFlux269 + ReactionFlux272;

    //d(V_T.T16)/dt
    NV_DATA_S(ydot)[69] =  - ReactionFlux283 - ReactionFlux285 - ReactionFlux286 + ReactionFlux289;

    //d(V_T.T17)/dt
    NV_DATA_S(ydot)[70] =  - ReactionFlux300 - ReactionFlux302 - ReactionFlux303 + ReactionFlux306;

    //d(V_T.APC)/dt
    NV_DATA_S(ydot)[71] = ReactionFlux310 - ReactionFlux312;

    //d(V_T.mAPC)/dt
    NV_DATA_S(ydot)[72] = ReactionFlux312 - ReactionFlux313 - ReactionFlux314;

    //d(V_T.c)/dt
    NV_DATA_S(ydot)[73] = 1/AUX_VAR_V_T*(ReactionFlux316 + ReactionFlux317);

    //d(V_T.nivo)/dt
    NV_DATA_S(ydot)[74] = 1/AUX_VAR_V_T*(ReactionFlux500 - ReactionFlux502);

    //d(V_T.durv)/dt
    NV_DATA_S(ydot)[75] = 1/AUX_VAR_V_T*(ReactionFlux506 - ReactionFlux508);

    //d(V_T.ipi)/dt
    NV_DATA_S(ydot)[76] = 1/AUX_VAR_V_T*(ReactionFlux512 - ReactionFlux514);

    //d(V_T.Treg_CTLA4)/dt
    NV_DATA_S(ydot)[77] =  - ReactionFlux561 - ReactionFlux562;

    //d(V_T.Treg_CTLA4_ipi)/dt
    NV_DATA_S(ydot)[78] = ReactionFlux561 - ReactionFlux562;

    //d(V_T.Treg_CTLA4_ipi_CTLA4)/dt
    NV_DATA_S(ydot)[79] = ReactionFlux562;

    //d(V_T.MDSC)/dt
    NV_DATA_S(ydot)[80] = ReactionFlux567 + ReactionFlux568 - ReactionFlux569;

    //d(V_T.CCL2)/dt
    NV_DATA_S(ydot)[81] = 1/AUX_VAR_V_T*( - ReactionFlux570 + ReactionFlux573);

    //d(V_T.NO)/dt
    NV_DATA_S(ydot)[82] = 1/AUX_VAR_V_T*( - ReactionFlux571 + ReactionFlux574);

    //d(V_T.ArgI)/dt
    NV_DATA_S(ydot)[83] = 1/AUX_VAR_V_T*( - ReactionFlux572 + ReactionFlux575);

    //d(V_T.ENT)/dt
    NV_DATA_S(ydot)[84] = 1/AUX_VAR_V_T*(ReactionFlux582 - ReactionFlux584);

    //d(V_LN.nT0)/dt
    NV_DATA_S(ydot)[85] = ReactionFlux5 - ReactionFlux6 - ReactionFlux7;

    //d(V_LN.aT0)/dt
    NV_DATA_S(ydot)[86] = ReactionFlux7 - ReactionFlux8;

    //d(V_LN.T0)/dt
    NV_DATA_S(ydot)[87] = ReactionFlux9 - ReactionFlux13 - ReactionFlux17;

    //d(V_LN.IL2)/dt
    NV_DATA_S(ydot)[88] = 1/PARAM(2)*( - ReactionFlux18 - ReactionFlux19 - ReactionFlux20 + ReactionFlux36 + ReactionFlux53 + ReactionFlux70 + ReactionFlux87 + ReactionFlux104 + ReactionFlux121 + ReactionFlux138 + ReactionFlux155 + ReactionFlux172 + ReactionFlux189 + ReactionFlux206 + ReactionFlux223 + ReactionFlux240 + ReactionFlux257 + ReactionFlux274 + ReactionFlux291 + ReactionFlux308);

    //d(V_LN.nT1)/dt
    NV_DATA_S(ydot)[89] = ReactionFlux21 - ReactionFlux22 - ReactionFlux23;

    //d(V_LN.aT1)/dt
    NV_DATA_S(ydot)[90] = ReactionFlux23 - ReactionFlux24;

    //d(V_LN.T1)/dt
    NV_DATA_S(ydot)[91] = ReactionFlux25 - ReactionFlux29 - ReactionFlux35;

    //d(V_LN.nT2)/dt
    NV_DATA_S(ydot)[92] = ReactionFlux38 - ReactionFlux39 - ReactionFlux40;

    //d(V_LN.aT2)/dt
    NV_DATA_S(ydot)[93] = ReactionFlux40 - ReactionFlux41;

    //d(V_LN.T2)/dt
    NV_DATA_S(ydot)[94] = ReactionFlux42 - ReactionFlux46 - ReactionFlux52;

    //d(V_LN.nT3)/dt
    NV_DATA_S(ydot)[95] = ReactionFlux55 - ReactionFlux56 - ReactionFlux57;

    //d(V_LN.aT3)/dt
    NV_DATA_S(ydot)[96] = ReactionFlux57 - ReactionFlux58;

    //d(V_LN.T3)/dt
    NV_DATA_S(ydot)[97] = ReactionFlux59 - ReactionFlux63 - ReactionFlux69;

    //d(V_LN.nT4)/dt
    NV_DATA_S(ydot)[98] = ReactionFlux72 - ReactionFlux73 - ReactionFlux74;

    //d(V_LN.aT4)/dt
    NV_DATA_S(ydot)[99] = ReactionFlux74 - ReactionFlux75;

    //d(V_LN.T4)/dt
    NV_DATA_S(ydot)[100] = ReactionFlux76 - ReactionFlux80 - ReactionFlux86;

    //d(V_LN.nT5)/dt
    NV_DATA_S(ydot)[101] = ReactionFlux89 - ReactionFlux90 - ReactionFlux91;

    //d(V_LN.aT5)/dt
    NV_DATA_S(ydot)[102] = ReactionFlux91 - ReactionFlux92;

    //d(V_LN.T5)/dt
    NV_DATA_S(ydot)[103] = ReactionFlux93 - ReactionFlux97 - ReactionFlux103;

    //d(V_LN.nT6)/dt
    NV_DATA_S(ydot)[104] = ReactionFlux106 - ReactionFlux107 - ReactionFlux108;

    //d(V_LN.aT6)/dt
    NV_DATA_S(ydot)[105] = ReactionFlux108 - ReactionFlux109;

    //d(V_LN.T6)/dt
    NV_DATA_S(ydot)[106] = ReactionFlux110 - ReactionFlux114 - ReactionFlux120;

    //d(V_LN.nT7)/dt
    NV_DATA_S(ydot)[107] = ReactionFlux123 - ReactionFlux124 - ReactionFlux125;

    //d(V_LN.aT7)/dt
    NV_DATA_S(ydot)[108] = ReactionFlux125 - ReactionFlux126;

    //d(V_LN.T7)/dt
    NV_DATA_S(ydot)[109] = ReactionFlux127 - ReactionFlux131 - ReactionFlux137;

    //d(V_LN.nT8)/dt
    NV_DATA_S(ydot)[110] = ReactionFlux140 - ReactionFlux141 - ReactionFlux142;

    //d(V_LN.aT8)/dt
    NV_DATA_S(ydot)[111] = ReactionFlux142 - ReactionFlux143;

    //d(V_LN.T8)/dt
    NV_DATA_S(ydot)[112] = ReactionFlux144 - ReactionFlux148 - ReactionFlux154;

    //d(V_LN.nT9)/dt
    NV_DATA_S(ydot)[113] = ReactionFlux157 - ReactionFlux158 - ReactionFlux159;

    //d(V_LN.aT9)/dt
    NV_DATA_S(ydot)[114] = ReactionFlux159 - ReactionFlux160;

    //d(V_LN.T9)/dt
    NV_DATA_S(ydot)[115] = ReactionFlux161 - ReactionFlux165 - ReactionFlux171;

    //d(V_LN.nT10)/dt
    NV_DATA_S(ydot)[116] = ReactionFlux174 - ReactionFlux175 - ReactionFlux176;

    //d(V_LN.aT10)/dt
    NV_DATA_S(ydot)[117] = ReactionFlux176 - ReactionFlux177;

    //d(V_LN.T10)/dt
    NV_DATA_S(ydot)[118] = ReactionFlux178 - ReactionFlux182 - ReactionFlux188;

    //d(V_LN.nT11)/dt
    NV_DATA_S(ydot)[119] = ReactionFlux191 - ReactionFlux192 - ReactionFlux193;

    //d(V_LN.aT11)/dt
    NV_DATA_S(ydot)[120] = ReactionFlux193 - ReactionFlux194;

    //d(V_LN.T11)/dt
    NV_DATA_S(ydot)[121] = ReactionFlux195 - ReactionFlux199 - ReactionFlux205;

    //d(V_LN.nT12)/dt
    NV_DATA_S(ydot)[122] = ReactionFlux208 - ReactionFlux209 - ReactionFlux210;

    //d(V_LN.aT12)/dt
    NV_DATA_S(ydot)[123] = ReactionFlux210 - ReactionFlux211;

    //d(V_LN.T12)/dt
    NV_DATA_S(ydot)[124] = ReactionFlux212 - ReactionFlux216 - ReactionFlux222;

    //d(V_LN.nT13)/dt
    NV_DATA_S(ydot)[125] = ReactionFlux225 - ReactionFlux226 - ReactionFlux227;

    //d(V_LN.aT13)/dt
    NV_DATA_S(ydot)[126] = ReactionFlux227 - ReactionFlux228;

    //d(V_LN.T13)/dt
    NV_DATA_S(ydot)[127] = ReactionFlux229 - ReactionFlux233 - ReactionFlux239;

    //d(V_LN.nT14)/dt
    NV_DATA_S(ydot)[128] = ReactionFlux242 - ReactionFlux243 - ReactionFlux244;

    //d(V_LN.aT14)/dt
    NV_DATA_S(ydot)[129] = ReactionFlux244 - ReactionFlux245;

    //d(V_LN.T14)/dt
    NV_DATA_S(ydot)[130] = ReactionFlux246 - ReactionFlux250 - ReactionFlux256;

    //d(V_LN.nT15)/dt
    NV_DATA_S(ydot)[131] = ReactionFlux259 - ReactionFlux260 - ReactionFlux261;

    //d(V_LN.aT15)/dt
    NV_DATA_S(ydot)[132] = ReactionFlux261 - ReactionFlux262;

    //d(V_LN.T15)/dt
    NV_DATA_S(ydot)[133] = ReactionFlux263 - ReactionFlux267 - ReactionFlux273;

    //d(V_LN.nT16)/dt
    NV_DATA_S(ydot)[134] = ReactionFlux276 - ReactionFlux277 - ReactionFlux278;

    //d(V_LN.aT16)/dt
    NV_DATA_S(ydot)[135] = ReactionFlux278 - ReactionFlux279;

    //d(V_LN.T16)/dt
    NV_DATA_S(ydot)[136] = ReactionFlux280 - ReactionFlux284 - ReactionFlux290;

    //d(V_LN.nT17)/dt
    NV_DATA_S(ydot)[137] = ReactionFlux293 - ReactionFlux294 - ReactionFlux295;

    //d(V_LN.aT17)/dt
    NV_DATA_S(ydot)[138] = ReactionFlux295 - ReactionFlux296;

    //d(V_LN.T17)/dt
    NV_DATA_S(ydot)[139] = ReactionFlux297 - ReactionFlux301 - ReactionFlux307;

    //d(V_LN.APC)/dt
    NV_DATA_S(ydot)[140] = ReactionFlux311;

    //d(V_LN.mAPC)/dt
    NV_DATA_S(ydot)[141] = ReactionFlux313 - ReactionFlux315;

    //d(V_LN.P0)/dt
    NV_DATA_S(ydot)[142] = 1/PARAM(2)*(ReactionFlux319 - ReactionFlux320 - ReactionFlux321);

    //d(V_LN.P1)/dt
    NV_DATA_S(ydot)[143] = 1/PARAM(2)*(ReactionFlux329 - ReactionFlux330 - ReactionFlux331);

    //d(V_LN.P2)/dt
    NV_DATA_S(ydot)[144] = 1/PARAM(2)*(ReactionFlux339 - ReactionFlux340 - ReactionFlux341);

    //d(V_LN.P3)/dt
    NV_DATA_S(ydot)[145] = 1/PARAM(2)*(ReactionFlux349 - ReactionFlux350 - ReactionFlux351);

    //d(V_LN.P4)/dt
    NV_DATA_S(ydot)[146] = 1/PARAM(2)*(ReactionFlux359 - ReactionFlux360 - ReactionFlux361);

    //d(V_LN.P5)/dt
    NV_DATA_S(ydot)[147] = 1/PARAM(2)*(ReactionFlux369 - ReactionFlux370 - ReactionFlux371);

    //d(V_LN.P6)/dt
    NV_DATA_S(ydot)[148] = 1/PARAM(2)*(ReactionFlux379 - ReactionFlux380 - ReactionFlux381);

    //d(V_LN.P7)/dt
    NV_DATA_S(ydot)[149] = 1/PARAM(2)*(ReactionFlux389 - ReactionFlux390 - ReactionFlux391);

    //d(V_LN.P8)/dt
    NV_DATA_S(ydot)[150] = 1/PARAM(2)*(ReactionFlux399 - ReactionFlux400 - ReactionFlux401);

    //d(V_LN.P9)/dt
    NV_DATA_S(ydot)[151] = 1/PARAM(2)*(ReactionFlux409 - ReactionFlux410 - ReactionFlux411);

    //d(V_LN.P10)/dt
    NV_DATA_S(ydot)[152] = 1/PARAM(2)*(ReactionFlux419 - ReactionFlux420 - ReactionFlux421);

    //d(V_LN.P11)/dt
    NV_DATA_S(ydot)[153] = 1/PARAM(2)*(ReactionFlux429 - ReactionFlux430 - ReactionFlux431);

    //d(V_LN.P12)/dt
    NV_DATA_S(ydot)[154] = 1/PARAM(2)*(ReactionFlux439 - ReactionFlux440 - ReactionFlux441);

    //d(V_LN.P13)/dt
    NV_DATA_S(ydot)[155] = 1/PARAM(2)*(ReactionFlux449 - ReactionFlux450 - ReactionFlux451);

    //d(V_LN.P14)/dt
    NV_DATA_S(ydot)[156] = 1/PARAM(2)*(ReactionFlux459 - ReactionFlux460 - ReactionFlux461);

    //d(V_LN.P15)/dt
    NV_DATA_S(ydot)[157] = 1/PARAM(2)*(ReactionFlux469 - ReactionFlux470 - ReactionFlux471);

    //d(V_LN.P16)/dt
    NV_DATA_S(ydot)[158] = 1/PARAM(2)*(ReactionFlux479 - ReactionFlux480 - ReactionFlux481);

    //d(V_LN.P17)/dt
    NV_DATA_S(ydot)[159] = 1/PARAM(2)*(ReactionFlux489 - ReactionFlux490 - ReactionFlux491);

    //d(V_LN.nivo)/dt
    NV_DATA_S(ydot)[160] = 1/PARAM(2)*(ReactionFlux501 + ReactionFlux502 - ReactionFlux503);

    //d(V_LN.durv)/dt
    NV_DATA_S(ydot)[161] = 1/PARAM(2)*(ReactionFlux507 + ReactionFlux508 - ReactionFlux509);

    //d(V_LN.ipi)/dt
    NV_DATA_S(ydot)[162] = 1/PARAM(2)*(ReactionFlux513 + ReactionFlux514 - ReactionFlux515);

    //d(V_LN.ENT)/dt
    NV_DATA_S(ydot)[163] = 1/PARAM(2)*(ReactionFlux583 + ReactionFlux584 - ReactionFlux585);

    //d(V_e.P0)/dt
    NV_DATA_S(ydot)[164] = 1/PARAM(3)*(ReactionFlux322 - ReactionFlux323);

    //d(V_e.p0)/dt
    NV_DATA_S(ydot)[165] = 1/PARAM(3)*(ReactionFlux323 - ReactionFlux324 - ReactionFlux325 + ReactionFlux326);

    //d(V_e.P1)/dt
    NV_DATA_S(ydot)[166] = 1/PARAM(3)*(ReactionFlux332 - ReactionFlux333);

    //d(V_e.p1)/dt
    NV_DATA_S(ydot)[167] = 1/PARAM(3)*(ReactionFlux333 - ReactionFlux334 - ReactionFlux335 + ReactionFlux336);

    //d(V_e.P2)/dt
    NV_DATA_S(ydot)[168] = 1/PARAM(3)*(ReactionFlux342 - ReactionFlux343);

    //d(V_e.p2)/dt
    NV_DATA_S(ydot)[169] = 1/PARAM(3)*(ReactionFlux343 - ReactionFlux344 - ReactionFlux345 + ReactionFlux346);

    //d(V_e.P3)/dt
    NV_DATA_S(ydot)[170] = 1/PARAM(3)*(ReactionFlux352 - ReactionFlux353);

    //d(V_e.p3)/dt
    NV_DATA_S(ydot)[171] = 1/PARAM(3)*(ReactionFlux353 - ReactionFlux354 - ReactionFlux355 + ReactionFlux356);

    //d(V_e.P4)/dt
    NV_DATA_S(ydot)[172] = 1/PARAM(3)*(ReactionFlux362 - ReactionFlux363);

    //d(V_e.p4)/dt
    NV_DATA_S(ydot)[173] = 1/PARAM(3)*(ReactionFlux363 - ReactionFlux364 - ReactionFlux365 + ReactionFlux366);

    //d(V_e.P5)/dt
    NV_DATA_S(ydot)[174] = 1/PARAM(3)*(ReactionFlux372 - ReactionFlux373);

    //d(V_e.p5)/dt
    NV_DATA_S(ydot)[175] = 1/PARAM(3)*(ReactionFlux373 - ReactionFlux374 - ReactionFlux375 + ReactionFlux376);

    //d(V_e.P6)/dt
    NV_DATA_S(ydot)[176] = 1/PARAM(3)*(ReactionFlux382 - ReactionFlux383);

    //d(V_e.p6)/dt
    NV_DATA_S(ydot)[177] = 1/PARAM(3)*(ReactionFlux383 - ReactionFlux384 - ReactionFlux385 + ReactionFlux386);

    //d(V_e.P7)/dt
    NV_DATA_S(ydot)[178] = 1/PARAM(3)*(ReactionFlux392 - ReactionFlux393);

    //d(V_e.p7)/dt
    NV_DATA_S(ydot)[179] = 1/PARAM(3)*(ReactionFlux393 - ReactionFlux394 - ReactionFlux395 + ReactionFlux396);

    //d(V_e.P8)/dt
    NV_DATA_S(ydot)[180] = 1/PARAM(3)*(ReactionFlux402 - ReactionFlux403);

    //d(V_e.p8)/dt
    NV_DATA_S(ydot)[181] = 1/PARAM(3)*(ReactionFlux403 - ReactionFlux404 - ReactionFlux405 + ReactionFlux406);

    //d(V_e.P9)/dt
    NV_DATA_S(ydot)[182] = 1/PARAM(3)*(ReactionFlux412 - ReactionFlux413);

    //d(V_e.p9)/dt
    NV_DATA_S(ydot)[183] = 1/PARAM(3)*(ReactionFlux413 - ReactionFlux414 - ReactionFlux415 + ReactionFlux416);

    //d(V_e.P10)/dt
    NV_DATA_S(ydot)[184] = 1/PARAM(3)*(ReactionFlux422 - ReactionFlux423);

    //d(V_e.p10)/dt
    NV_DATA_S(ydot)[185] = 1/PARAM(3)*(ReactionFlux423 - ReactionFlux424 - ReactionFlux425 + ReactionFlux426);

    //d(V_e.P11)/dt
    NV_DATA_S(ydot)[186] = 1/PARAM(3)*(ReactionFlux432 - ReactionFlux433);

    //d(V_e.p11)/dt
    NV_DATA_S(ydot)[187] = 1/PARAM(3)*(ReactionFlux433 - ReactionFlux434 - ReactionFlux435 + ReactionFlux436);

    //d(V_e.P12)/dt
    NV_DATA_S(ydot)[188] = 1/PARAM(3)*(ReactionFlux442 - ReactionFlux443);

    //d(V_e.p12)/dt
    NV_DATA_S(ydot)[189] = 1/PARAM(3)*(ReactionFlux443 - ReactionFlux444 - ReactionFlux445 + ReactionFlux446);

    //d(V_e.P13)/dt
    NV_DATA_S(ydot)[190] = 1/PARAM(3)*(ReactionFlux452 - ReactionFlux453);

    //d(V_e.p13)/dt
    NV_DATA_S(ydot)[191] = 1/PARAM(3)*(ReactionFlux453 - ReactionFlux454 - ReactionFlux455 + ReactionFlux456);

    //d(V_e.P14)/dt
    NV_DATA_S(ydot)[192] = 1/PARAM(3)*(ReactionFlux462 - ReactionFlux463);

    //d(V_e.p14)/dt
    NV_DATA_S(ydot)[193] = 1/PARAM(3)*(ReactionFlux463 - ReactionFlux464 - ReactionFlux465 + ReactionFlux466);

    //d(V_e.P15)/dt
    NV_DATA_S(ydot)[194] = 1/PARAM(3)*(ReactionFlux472 - ReactionFlux473);

    //d(V_e.p15)/dt
    NV_DATA_S(ydot)[195] = 1/PARAM(3)*(ReactionFlux473 - ReactionFlux474 - ReactionFlux475 + ReactionFlux476);

    //d(V_e.P16)/dt
    NV_DATA_S(ydot)[196] = 1/PARAM(3)*(ReactionFlux482 - ReactionFlux483);

    //d(V_e.p16)/dt
    NV_DATA_S(ydot)[197] = 1/PARAM(3)*(ReactionFlux483 - ReactionFlux484 - ReactionFlux485 + ReactionFlux486);

    //d(V_e.P17)/dt
    NV_DATA_S(ydot)[198] = 1/PARAM(3)*(ReactionFlux492 - ReactionFlux493);

    //d(V_e.p17)/dt
    NV_DATA_S(ydot)[199] = 1/PARAM(3)*(ReactionFlux493 - ReactionFlux494 - ReactionFlux495 + ReactionFlux496);

    //d(A_e.M1)/dt
    NV_DATA_S(ydot)[200] = 1/PARAM(4)*( - ReactionFlux318 - ReactionFlux325 + ReactionFlux326 - ReactionFlux335 + ReactionFlux336 - ReactionFlux345 + ReactionFlux346 - ReactionFlux355 + ReactionFlux356 - ReactionFlux365 + ReactionFlux366 - ReactionFlux375 + ReactionFlux376 - ReactionFlux385 + ReactionFlux386 - ReactionFlux395 + ReactionFlux396 - ReactionFlux405 + ReactionFlux406 - ReactionFlux415 + ReactionFlux416 - ReactionFlux425 + ReactionFlux426 - ReactionFlux435 + ReactionFlux436 - ReactionFlux445 + ReactionFlux446 - ReactionFlux455 + ReactionFlux456 - ReactionFlux465 + ReactionFlux466 - ReactionFlux475 + ReactionFlux476 - ReactionFlux485 + ReactionFlux486 - ReactionFlux495 + ReactionFlux496);

    //d(A_e.M1p0)/dt
    NV_DATA_S(ydot)[201] = 1/PARAM(4)*(ReactionFlux325 - ReactionFlux326 - ReactionFlux328);

    //d(A_e.M1p1)/dt
    NV_DATA_S(ydot)[202] = 1/PARAM(4)*(ReactionFlux335 - ReactionFlux336 - ReactionFlux338);

    //d(A_e.M1p2)/dt
    NV_DATA_S(ydot)[203] = 1/PARAM(4)*(ReactionFlux345 - ReactionFlux346 - ReactionFlux348);

    //d(A_e.M1p3)/dt
    NV_DATA_S(ydot)[204] = 1/PARAM(4)*(ReactionFlux355 - ReactionFlux356 - ReactionFlux358);

    //d(A_e.M1p4)/dt
    NV_DATA_S(ydot)[205] = 1/PARAM(4)*(ReactionFlux365 - ReactionFlux366 - ReactionFlux368);

    //d(A_e.M1p5)/dt
    NV_DATA_S(ydot)[206] = 1/PARAM(4)*(ReactionFlux375 - ReactionFlux376 - ReactionFlux378);

    //d(A_e.M1p6)/dt
    NV_DATA_S(ydot)[207] = 1/PARAM(4)*(ReactionFlux385 - ReactionFlux386 - ReactionFlux388);

    //d(A_e.M1p7)/dt
    NV_DATA_S(ydot)[208] = 1/PARAM(4)*(ReactionFlux395 - ReactionFlux396 - ReactionFlux398);

    //d(A_e.M1p8)/dt
    NV_DATA_S(ydot)[209] = 1/PARAM(4)*(ReactionFlux405 - ReactionFlux406 - ReactionFlux408);

    //d(A_e.M1p9)/dt
    NV_DATA_S(ydot)[210] = 1/PARAM(4)*(ReactionFlux415 - ReactionFlux416 - ReactionFlux418);

    //d(A_e.M1p10)/dt
    NV_DATA_S(ydot)[211] = 1/PARAM(4)*(ReactionFlux425 - ReactionFlux426 - ReactionFlux428);

    //d(A_e.M1p11)/dt
    NV_DATA_S(ydot)[212] = 1/PARAM(4)*(ReactionFlux435 - ReactionFlux436 - ReactionFlux438);

    //d(A_e.M1p12)/dt
    NV_DATA_S(ydot)[213] = 1/PARAM(4)*(ReactionFlux445 - ReactionFlux446 - ReactionFlux448);

    //d(A_e.M1p13)/dt
    NV_DATA_S(ydot)[214] = 1/PARAM(4)*(ReactionFlux455 - ReactionFlux456 - ReactionFlux458);

    //d(A_e.M1p14)/dt
    NV_DATA_S(ydot)[215] = 1/PARAM(4)*(ReactionFlux465 - ReactionFlux466 - ReactionFlux468);

    //d(A_e.M1p15)/dt
    NV_DATA_S(ydot)[216] = 1/PARAM(4)*(ReactionFlux475 - ReactionFlux476 - ReactionFlux478);

    //d(A_e.M1p16)/dt
    NV_DATA_S(ydot)[217] = 1/PARAM(4)*(ReactionFlux485 - ReactionFlux486 - ReactionFlux488);

    //d(A_e.M1p17)/dt
    NV_DATA_S(ydot)[218] = 1/PARAM(4)*(ReactionFlux495 - ReactionFlux496 - ReactionFlux498);

    //d(A_s.M1)/dt
    NV_DATA_S(ydot)[219] = 1/PARAM(5)*(ReactionFlux318 + ReactionFlux327 + ReactionFlux337 + ReactionFlux347 + ReactionFlux357 + ReactionFlux367 + ReactionFlux377 + ReactionFlux387 + ReactionFlux397 + ReactionFlux407 + ReactionFlux417 + ReactionFlux427 + ReactionFlux437 + ReactionFlux447 + ReactionFlux457 + ReactionFlux467 + ReactionFlux477 + ReactionFlux487 + ReactionFlux497);

    //d(A_s.M1p0)/dt
    NV_DATA_S(ydot)[220] = 1/PARAM(5)*( - ReactionFlux327 + ReactionFlux328);

    //d(A_s.M1p1)/dt
    NV_DATA_S(ydot)[221] = 1/PARAM(5)*( - ReactionFlux337 + ReactionFlux338);

    //d(A_s.M1p2)/dt
    NV_DATA_S(ydot)[222] = 1/PARAM(5)*( - ReactionFlux347 + ReactionFlux348);

    //d(A_s.M1p3)/dt
    NV_DATA_S(ydot)[223] = 1/PARAM(5)*( - ReactionFlux357 + ReactionFlux358);

    //d(A_s.M1p4)/dt
    NV_DATA_S(ydot)[224] = 1/PARAM(5)*( - ReactionFlux367 + ReactionFlux368);

    //d(A_s.M1p5)/dt
    NV_DATA_S(ydot)[225] = 1/PARAM(5)*( - ReactionFlux377 + ReactionFlux378);

    //d(A_s.M1p6)/dt
    NV_DATA_S(ydot)[226] = 1/PARAM(5)*( - ReactionFlux387 + ReactionFlux388);

    //d(A_s.M1p7)/dt
    NV_DATA_S(ydot)[227] = 1/PARAM(5)*( - ReactionFlux397 + ReactionFlux398);

    //d(A_s.M1p8)/dt
    NV_DATA_S(ydot)[228] = 1/PARAM(5)*( - ReactionFlux407 + ReactionFlux408);

    //d(A_s.M1p9)/dt
    NV_DATA_S(ydot)[229] = 1/PARAM(5)*( - ReactionFlux417 + ReactionFlux418);

    //d(A_s.M1p10)/dt
    NV_DATA_S(ydot)[230] = 1/PARAM(5)*( - ReactionFlux427 + ReactionFlux428);

    //d(A_s.M1p11)/dt
    NV_DATA_S(ydot)[231] = 1/PARAM(5)*( - ReactionFlux437 + ReactionFlux438);

    //d(A_s.M1p12)/dt
    NV_DATA_S(ydot)[232] = 1/PARAM(5)*( - ReactionFlux447 + ReactionFlux448);

    //d(A_s.M1p13)/dt
    NV_DATA_S(ydot)[233] = 1/PARAM(5)*( - ReactionFlux457 + ReactionFlux458);

    //d(A_s.M1p14)/dt
    NV_DATA_S(ydot)[234] = 1/PARAM(5)*( - ReactionFlux467 + ReactionFlux468);

    //d(A_s.M1p15)/dt
    NV_DATA_S(ydot)[235] = 1/PARAM(5)*( - ReactionFlux477 + ReactionFlux478);

    //d(A_s.M1p16)/dt
    NV_DATA_S(ydot)[236] = 1/PARAM(5)*( - ReactionFlux487 + ReactionFlux488);

    //d(A_s.M1p17)/dt
    NV_DATA_S(ydot)[237] = 1/PARAM(5)*( - ReactionFlux497 + ReactionFlux498);

    //d(syn_T_C1.PD1_PDL1)/dt
    NV_DATA_S(ydot)[238] = 1/PARAM(6)*(ReactionFlux517);

    //d(syn_T_C1.PD1_PDL2)/dt
    NV_DATA_S(ydot)[239] = 1/PARAM(6)*(ReactionFlux518);

    //d(syn_T_C1.PD1)/dt
    NV_DATA_S(ydot)[240] = 1/PARAM(6)*( - ReactionFlux517 - ReactionFlux518 - ReactionFlux519 - ReactionFlux520);

    //d(syn_T_C1.PDL1)/dt
    NV_DATA_S(ydot)[241] = 1/PARAM(6)*( - ReactionFlux517 - ReactionFlux521 - ReactionFlux522);

    //d(syn_T_C1.PDL2)/dt
    NV_DATA_S(ydot)[242] = 1/PARAM(6)*( - ReactionFlux518);

    //d(syn_T_C1.PD1_nivo)/dt
    NV_DATA_S(ydot)[243] = 1/PARAM(6)*(ReactionFlux519 - ReactionFlux520);

    //d(syn_T_C1.PD1_nivo_PD1)/dt
    NV_DATA_S(ydot)[244] = 1/PARAM(6)*(ReactionFlux520);

    //d(syn_T_C1.PDL1_durv)/dt
    NV_DATA_S(ydot)[245] = 1/PARAM(6)*(ReactionFlux521 - ReactionFlux522);

    //d(syn_T_C1.PDL1_durv_PDL1)/dt
    NV_DATA_S(ydot)[246] = 1/PARAM(6)*(ReactionFlux522);

    //d(syn_T_C1.TPDL1)/dt
    NV_DATA_S(ydot)[247] = 1/PARAM(6)*( - ReactionFlux535 - ReactionFlux536 - ReactionFlux537 - ReactionFlux538);

    //d(syn_T_C1.TPDL1_durv)/dt
    NV_DATA_S(ydot)[248] = 1/PARAM(6)*(ReactionFlux537 - ReactionFlux538);

    //d(syn_T_C1.TPDL1_durv_TPDL1)/dt
    NV_DATA_S(ydot)[249] = 1/PARAM(6)*(ReactionFlux538);

    //d(syn_T_C1.CD28_CD80)/dt
    NV_DATA_S(ydot)[250] = 1/PARAM(6)*(ReactionFlux523 - ReactionFlux524);

    //d(syn_T_C1.CD28_CD80_CD28)/dt
    NV_DATA_S(ydot)[251] = 1/PARAM(6)*(ReactionFlux524);

    //d(syn_T_C1.CD28_CD86)/dt
    NV_DATA_S(ydot)[252] = 1/PARAM(6)*(ReactionFlux525);

    //d(syn_T_C1.CD80_CTLA4)/dt
    NV_DATA_S(ydot)[253] = 1/PARAM(6)*(ReactionFlux526 - ReactionFlux527 - ReactionFlux529);

    //d(syn_T_C1.CD80_CTLA4_CD80)/dt
    NV_DATA_S(ydot)[254] = 1/PARAM(6)*(ReactionFlux529 - ReactionFlux530);

    //d(syn_T_C1.CTLA4_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[255] = 1/PARAM(6)*(ReactionFlux527 - ReactionFlux528);

    //d(syn_T_C1.CD80_CTLA4_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[256] = 1/PARAM(6)*(ReactionFlux528 + ReactionFlux530);

    //d(syn_T_C1.CD86_CTLA4)/dt
    NV_DATA_S(ydot)[257] = 1/PARAM(6)*(ReactionFlux531 - ReactionFlux532);

    //d(syn_T_C1.CD86_CTLA4_CD86)/dt
    NV_DATA_S(ydot)[258] = 1/PARAM(6)*(ReactionFlux532);

    //d(syn_T_C1.TPDL1_CD80)/dt
    NV_DATA_S(ydot)[259] = 1/PARAM(6)*(ReactionFlux535 - ReactionFlux536);

    //d(syn_T_C1.TPDL1_CD80_TPDL1)/dt
    NV_DATA_S(ydot)[260] = 1/PARAM(6)*(ReactionFlux536);

    //d(syn_T_C1.CD28)/dt
    NV_DATA_S(ydot)[261] = 1/PARAM(6)*( - ReactionFlux523 - ReactionFlux524 - ReactionFlux525);

    //d(syn_T_C1.CTLA4)/dt
    NV_DATA_S(ydot)[262] = 1/PARAM(6)*( - ReactionFlux526 - ReactionFlux527 - ReactionFlux530 - ReactionFlux531 - ReactionFlux533 - ReactionFlux534);

    //d(syn_T_C1.CD80)/dt
    NV_DATA_S(ydot)[263] = 1/PARAM(6)*( - ReactionFlux523 - ReactionFlux526 - ReactionFlux528 - ReactionFlux529 - ReactionFlux535);

    //d(syn_T_C1.CD86)/dt
    NV_DATA_S(ydot)[264] = 1/PARAM(6)*( - ReactionFlux525 - ReactionFlux531 - ReactionFlux532);

    //d(syn_T_C1.CTLA4_ipi)/dt
    NV_DATA_S(ydot)[265] = 1/PARAM(6)*(ReactionFlux533 - ReactionFlux534);

    //d(syn_T_C1.CTLA4_ipi_CTLA4)/dt
    NV_DATA_S(ydot)[266] = 1/PARAM(6)*(ReactionFlux534);

    //d(syn_T_APC.PD1_PDL1)/dt
    NV_DATA_S(ydot)[267] = 1/PARAM(7)*(ReactionFlux539);

    //d(syn_T_APC.PD1_PDL2)/dt
    NV_DATA_S(ydot)[268] = 1/PARAM(7)*(ReactionFlux540);

    //d(syn_T_APC.PD1)/dt
    NV_DATA_S(ydot)[269] = 1/PARAM(7)*( - ReactionFlux539 - ReactionFlux540 - ReactionFlux541 - ReactionFlux542);

    //d(syn_T_APC.PDL1)/dt
    NV_DATA_S(ydot)[270] = 1/PARAM(7)*( - ReactionFlux539 - ReactionFlux543 - ReactionFlux544);

    //d(syn_T_APC.PDL2)/dt
    NV_DATA_S(ydot)[271] = 1/PARAM(7)*( - ReactionFlux540);

    //d(syn_T_APC.PD1_nivo)/dt
    NV_DATA_S(ydot)[272] = 1/PARAM(7)*(ReactionFlux541 - ReactionFlux542);

    //d(syn_T_APC.PD1_nivo_PD1)/dt
    NV_DATA_S(ydot)[273] = 1/PARAM(7)*(ReactionFlux542);

    //d(syn_T_APC.PDL1_durv)/dt
    NV_DATA_S(ydot)[274] = 1/PARAM(7)*(ReactionFlux543 - ReactionFlux544);

    //d(syn_T_APC.PDL1_durv_PDL1)/dt
    NV_DATA_S(ydot)[275] = 1/PARAM(7)*(ReactionFlux544);

    //d(syn_T_APC.TPDL1)/dt
    NV_DATA_S(ydot)[276] = 1/PARAM(7)*( - ReactionFlux557 - ReactionFlux558 - ReactionFlux559 - ReactionFlux560);

    //d(syn_T_APC.TPDL1_durv)/dt
    NV_DATA_S(ydot)[277] = 1/PARAM(7)*(ReactionFlux559 - ReactionFlux560);

    //d(syn_T_APC.TPDL1_durv_TPDL1)/dt
    NV_DATA_S(ydot)[278] = 1/PARAM(7)*(ReactionFlux560);

    //d(syn_T_APC.CD28_CD80)/dt
    NV_DATA_S(ydot)[279] = 1/PARAM(7)*(ReactionFlux545 - ReactionFlux546);

    //d(syn_T_APC.CD28_CD80_CD28)/dt
    NV_DATA_S(ydot)[280] = 1/PARAM(7)*(ReactionFlux546);

    //d(syn_T_APC.CD28_CD86)/dt
    NV_DATA_S(ydot)[281] = 1/PARAM(7)*(ReactionFlux547);

    //d(syn_T_APC.CD80_CTLA4)/dt
    NV_DATA_S(ydot)[282] = 1/PARAM(7)*(ReactionFlux548 - ReactionFlux549 - ReactionFlux551);

    //d(syn_T_APC.CD80_CTLA4_CD80)/dt
    NV_DATA_S(ydot)[283] = 1/PARAM(7)*(ReactionFlux551 - ReactionFlux552);

    //d(syn_T_APC.CTLA4_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[284] = 1/PARAM(7)*(ReactionFlux549 - ReactionFlux550);

    //d(syn_T_APC.CD80_CTLA4_CD80_CTLA4)/dt
    NV_DATA_S(ydot)[285] = 1/PARAM(7)*(ReactionFlux550 + ReactionFlux552);

    //d(syn_T_APC.CD86_CTLA4)/dt
    NV_DATA_S(ydot)[286] = 1/PARAM(7)*(ReactionFlux553 - ReactionFlux554);

    //d(syn_T_APC.CD86_CTLA4_CD86)/dt
    NV_DATA_S(ydot)[287] = 1/PARAM(7)*(ReactionFlux554);

    //d(syn_T_APC.TPDL1_CD80)/dt
    NV_DATA_S(ydot)[288] = 1/PARAM(7)*(ReactionFlux557 - ReactionFlux558);

    //d(syn_T_APC.TPDL1_CD80_TPDL1)/dt
    NV_DATA_S(ydot)[289] = 1/PARAM(7)*(ReactionFlux558);

    //d(syn_T_APC.CD28)/dt
    NV_DATA_S(ydot)[290] = 1/PARAM(7)*( - ReactionFlux545 - ReactionFlux546 - ReactionFlux547);

    //d(syn_T_APC.CTLA4)/dt
    NV_DATA_S(ydot)[291] = 1/PARAM(7)*( - ReactionFlux548 - ReactionFlux549 - ReactionFlux552 - ReactionFlux553 - ReactionFlux555 - ReactionFlux556);

    //d(syn_T_APC.CD80)/dt
    NV_DATA_S(ydot)[292] = 1/PARAM(7)*( - ReactionFlux545 - ReactionFlux548 - ReactionFlux550 - ReactionFlux551 - ReactionFlux557);

    //d(syn_T_APC.CD86)/dt
    NV_DATA_S(ydot)[293] = 1/PARAM(7)*( - ReactionFlux547 - ReactionFlux553 - ReactionFlux554);

    //d(syn_T_APC.CTLA4_ipi)/dt
    NV_DATA_S(ydot)[294] = 1/PARAM(7)*(ReactionFlux555 - ReactionFlux556);

    //d(syn_T_APC.CTLA4_ipi_CTLA4)/dt
    NV_DATA_S(ydot)[295] = 1/PARAM(7)*(ReactionFlux556);

    return(0);
}
int ODE_system::g(realtype t, N_Vector y, realtype *gout, void *user_data){

    ODE_system* ptrOde = static_cast<ODE_system*>(user_data);

    //Assignment rules:

    //V_T.C1 < (0.9 * cell)
    gout[0] = 0.9 * PARAM(9) - (SPVAR(52));

    return(0);
}

bool ODE_system::triggerComponentEvaluate(int i, realtype t, bool curr) {

    bool discrete = false;
    realtype diff = 0;
    bool eval = false;
    //Assignment rules:

    switch(i)
    {
    case 0:
        //V_T.C1 < (0.9 * cell)
        diff = 0.9 * _class_parameter[9] - (NV_DATA_S(_y)[52]);
        break;
    default:
        break;
    }
    if (!discrete){
        eval = diff == 0 ? curr : (diff > 0);
    }
    return eval;
}

bool ODE_system::eventEvaluate(int i) {
    bool eval = false;
    switch(i)
    {
    case 0:
        eval = getSatisfied(0);
        break;
    default:
        break;
    }
    return eval;
}

bool ODE_system::eventExecution(int i, bool delayed, realtype& dt){

    bool setDelay = false;

    //Assignment rules:

    switch(i)
    {
    case 0:
        NV_DATA_S(_y)[52] = 0.0 * _class_parameter[9];
        break;
    default:
        break;
    }
    return setDelay;
}
void ODE_system::update_y_other(void){

    return;
}
std::string ODE_system::getHeader(){

    std::string s = "";
    s += ",V_C.T0";
    s += ",V_C.T1";
    s += ",V_C.T2";
    s += ",V_C.T3";
    s += ",V_C.T4";
    s += ",V_C.T5";
    s += ",V_C.T6";
    s += ",V_C.T7";
    s += ",V_C.T8";
    s += ",V_C.T9";
    s += ",V_C.T10";
    s += ",V_C.T11";
    s += ",V_C.T12";
    s += ",V_C.T13";
    s += ",V_C.T14";
    s += ",V_C.T15";
    s += ",V_C.T16";
    s += ",V_C.T17";
    s += ",V_C.nivo";
    s += ",V_C.durv";
    s += ",V_C.ipi";
    s += ",V_C.ENT";
    s += ",V_C.ENT_Buccal";
    s += ",V_C.ENT_GI";
    s += ",V_C.Dose2";
    s += ",V_P.T0";
    s += ",V_P.T1";
    s += ",V_P.T2";
    s += ",V_P.T3";
    s += ",V_P.T4";
    s += ",V_P.T5";
    s += ",V_P.T6";
    s += ",V_P.T7";
    s += ",V_P.T8";
    s += ",V_P.T9";
    s += ",V_P.T10";
    s += ",V_P.T11";
    s += ",V_P.T12";
    s += ",V_P.T13";
    s += ",V_P.T14";
    s += ",V_P.T15";
    s += ",V_P.T16";
    s += ",V_P.T17";
    s += ",V_P.nivo";
    s += ",V_P.durv";
    s += ",V_P.ipi";
    s += ",V_P.Treg_CTLA4";
    s += ",V_P.Treg_CTLA4_ipi";
    s += ",V_P.Treg_CTLA4_ipi_CTLA4";
    s += ",V_P.ENT";
    s += ",V_T.C_x";
    s += ",V_T.T_exh";
    s += ",V_T.C1";
    s += ",V_T.T0";
    s += ",V_T.T1";
    s += ",V_T.T2";
    s += ",V_T.T3";
    s += ",V_T.T4";
    s += ",V_T.T5";
    s += ",V_T.T6";
    s += ",V_T.T7";
    s += ",V_T.T8";
    s += ",V_T.T9";
    s += ",V_T.T10";
    s += ",V_T.T11";
    s += ",V_T.T12";
    s += ",V_T.T13";
    s += ",V_T.T14";
    s += ",V_T.T15";
    s += ",V_T.T16";
    s += ",V_T.T17";
    s += ",V_T.APC";
    s += ",V_T.mAPC";
    s += ",V_T.c";
    s += ",V_T.nivo";
    s += ",V_T.durv";
    s += ",V_T.ipi";
    s += ",V_T.Treg_CTLA4";
    s += ",V_T.Treg_CTLA4_ipi";
    s += ",V_T.Treg_CTLA4_ipi_CTLA4";
    s += ",V_T.MDSC";
    s += ",V_T.CCL2";
    s += ",V_T.NO";
    s += ",V_T.ArgI";
    s += ",V_T.ENT";
    s += ",V_LN.nT0";
    s += ",V_LN.aT0";
    s += ",V_LN.T0";
    s += ",V_LN.IL2";
    s += ",V_LN.nT1";
    s += ",V_LN.aT1";
    s += ",V_LN.T1";
    s += ",V_LN.nT2";
    s += ",V_LN.aT2";
    s += ",V_LN.T2";
    s += ",V_LN.nT3";
    s += ",V_LN.aT3";
    s += ",V_LN.T3";
    s += ",V_LN.nT4";
    s += ",V_LN.aT4";
    s += ",V_LN.T4";
    s += ",V_LN.nT5";
    s += ",V_LN.aT5";
    s += ",V_LN.T5";
    s += ",V_LN.nT6";
    s += ",V_LN.aT6";
    s += ",V_LN.T6";
    s += ",V_LN.nT7";
    s += ",V_LN.aT7";
    s += ",V_LN.T7";
    s += ",V_LN.nT8";
    s += ",V_LN.aT8";
    s += ",V_LN.T8";
    s += ",V_LN.nT9";
    s += ",V_LN.aT9";
    s += ",V_LN.T9";
    s += ",V_LN.nT10";
    s += ",V_LN.aT10";
    s += ",V_LN.T10";
    s += ",V_LN.nT11";
    s += ",V_LN.aT11";
    s += ",V_LN.T11";
    s += ",V_LN.nT12";
    s += ",V_LN.aT12";
    s += ",V_LN.T12";
    s += ",V_LN.nT13";
    s += ",V_LN.aT13";
    s += ",V_LN.T13";
    s += ",V_LN.nT14";
    s += ",V_LN.aT14";
    s += ",V_LN.T14";
    s += ",V_LN.nT15";
    s += ",V_LN.aT15";
    s += ",V_LN.T15";
    s += ",V_LN.nT16";
    s += ",V_LN.aT16";
    s += ",V_LN.T16";
    s += ",V_LN.nT17";
    s += ",V_LN.aT17";
    s += ",V_LN.T17";
    s += ",V_LN.APC";
    s += ",V_LN.mAPC";
    s += ",V_LN.P0";
    s += ",V_LN.P1";
    s += ",V_LN.P2";
    s += ",V_LN.P3";
    s += ",V_LN.P4";
    s += ",V_LN.P5";
    s += ",V_LN.P6";
    s += ",V_LN.P7";
    s += ",V_LN.P8";
    s += ",V_LN.P9";
    s += ",V_LN.P10";
    s += ",V_LN.P11";
    s += ",V_LN.P12";
    s += ",V_LN.P13";
    s += ",V_LN.P14";
    s += ",V_LN.P15";
    s += ",V_LN.P16";
    s += ",V_LN.P17";
    s += ",V_LN.nivo";
    s += ",V_LN.durv";
    s += ",V_LN.ipi";
    s += ",V_LN.ENT";
    s += ",V_e.P0";
    s += ",V_e.p0";
    s += ",V_e.P1";
    s += ",V_e.p1";
    s += ",V_e.P2";
    s += ",V_e.p2";
    s += ",V_e.P3";
    s += ",V_e.p3";
    s += ",V_e.P4";
    s += ",V_e.p4";
    s += ",V_e.P5";
    s += ",V_e.p5";
    s += ",V_e.P6";
    s += ",V_e.p6";
    s += ",V_e.P7";
    s += ",V_e.p7";
    s += ",V_e.P8";
    s += ",V_e.p8";
    s += ",V_e.P9";
    s += ",V_e.p9";
    s += ",V_e.P10";
    s += ",V_e.p10";
    s += ",V_e.P11";
    s += ",V_e.p11";
    s += ",V_e.P12";
    s += ",V_e.p12";
    s += ",V_e.P13";
    s += ",V_e.p13";
    s += ",V_e.P14";
    s += ",V_e.p14";
    s += ",V_e.P15";
    s += ",V_e.p15";
    s += ",V_e.P16";
    s += ",V_e.p16";
    s += ",V_e.P17";
    s += ",V_e.p17";
    s += ",A_e.M1";
    s += ",A_e.M1p0";
    s += ",A_e.M1p1";
    s += ",A_e.M1p2";
    s += ",A_e.M1p3";
    s += ",A_e.M1p4";
    s += ",A_e.M1p5";
    s += ",A_e.M1p6";
    s += ",A_e.M1p7";
    s += ",A_e.M1p8";
    s += ",A_e.M1p9";
    s += ",A_e.M1p10";
    s += ",A_e.M1p11";
    s += ",A_e.M1p12";
    s += ",A_e.M1p13";
    s += ",A_e.M1p14";
    s += ",A_e.M1p15";
    s += ",A_e.M1p16";
    s += ",A_e.M1p17";
    s += ",A_s.M1";
    s += ",A_s.M1p0";
    s += ",A_s.M1p1";
    s += ",A_s.M1p2";
    s += ",A_s.M1p3";
    s += ",A_s.M1p4";
    s += ",A_s.M1p5";
    s += ",A_s.M1p6";
    s += ",A_s.M1p7";
    s += ",A_s.M1p8";
    s += ",A_s.M1p9";
    s += ",A_s.M1p10";
    s += ",A_s.M1p11";
    s += ",A_s.M1p12";
    s += ",A_s.M1p13";
    s += ",A_s.M1p14";
    s += ",A_s.M1p15";
    s += ",A_s.M1p16";
    s += ",A_s.M1p17";
    s += ",syn_T_C1.PD1_PDL1";
    s += ",syn_T_C1.PD1_PDL2";
    s += ",syn_T_C1.PD1";
    s += ",syn_T_C1.PDL1";
    s += ",syn_T_C1.PDL2";
    s += ",syn_T_C1.PD1_nivo";
    s += ",syn_T_C1.PD1_nivo_PD1";
    s += ",syn_T_C1.PDL1_durv";
    s += ",syn_T_C1.PDL1_durv_PDL1";
    s += ",syn_T_C1.TPDL1";
    s += ",syn_T_C1.TPDL1_durv";
    s += ",syn_T_C1.TPDL1_durv_TPDL1";
    s += ",syn_T_C1.CD28_CD80";
    s += ",syn_T_C1.CD28_CD80_CD28";
    s += ",syn_T_C1.CD28_CD86";
    s += ",syn_T_C1.CD80_CTLA4";
    s += ",syn_T_C1.CD80_CTLA4_CD80";
    s += ",syn_T_C1.CTLA4_CD80_CTLA4";
    s += ",syn_T_C1.CD80_CTLA4_CD80_CTLA4";
    s += ",syn_T_C1.CD86_CTLA4";
    s += ",syn_T_C1.CD86_CTLA4_CD86";
    s += ",syn_T_C1.TPDL1_CD80";
    s += ",syn_T_C1.TPDL1_CD80_TPDL1";
    s += ",syn_T_C1.CD28";
    s += ",syn_T_C1.CTLA4";
    s += ",syn_T_C1.CD80";
    s += ",syn_T_C1.CD86";
    s += ",syn_T_C1.CTLA4_ipi";
    s += ",syn_T_C1.CTLA4_ipi_CTLA4";
    s += ",syn_T_APC.PD1_PDL1";
    s += ",syn_T_APC.PD1_PDL2";
    s += ",syn_T_APC.PD1";
    s += ",syn_T_APC.PDL1";
    s += ",syn_T_APC.PDL2";
    s += ",syn_T_APC.PD1_nivo";
    s += ",syn_T_APC.PD1_nivo_PD1";
    s += ",syn_T_APC.PDL1_durv";
    s += ",syn_T_APC.PDL1_durv_PDL1";
    s += ",syn_T_APC.TPDL1";
    s += ",syn_T_APC.TPDL1_durv";
    s += ",syn_T_APC.TPDL1_durv_TPDL1";
    s += ",syn_T_APC.CD28_CD80";
    s += ",syn_T_APC.CD28_CD80_CD28";
    s += ",syn_T_APC.CD28_CD86";
    s += ",syn_T_APC.CD80_CTLA4";
    s += ",syn_T_APC.CD80_CTLA4_CD80";
    s += ",syn_T_APC.CTLA4_CD80_CTLA4";
    s += ",syn_T_APC.CD80_CTLA4_CD80_CTLA4";
    s += ",syn_T_APC.CD86_CTLA4";
    s += ",syn_T_APC.CD86_CTLA4_CD86";
    s += ",syn_T_APC.TPDL1_CD80";
    s += ",syn_T_APC.TPDL1_CD80_TPDL1";
    s += ",syn_T_APC.CD28";
    s += ",syn_T_APC.CTLA4";
    s += ",syn_T_APC.CD80";
    s += ",syn_T_APC.CD86";
    s += ",syn_T_APC.CTLA4_ipi";
    s += ",syn_T_APC.CTLA4_ipi_CTLA4";
    return s;
}
realtype ODE_system::get_unit_conversion_species(int i) const{

    static std::vector<realtype> scalor = {
        //sp_var
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        999.9999999999999,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1000000.0,
        1000000.0,
        1000000.0,
        1000000.0,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1000000.0,
        1000000.0,
        1000000.0,
        1000000.0,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        999999999.9999999,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        1.66053872801495e-24,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999999999.9999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        999.9999999999999,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        1.66053872801495e-12,
        //sp_other
    };
    return scalor[i];
}
realtype ODE_system::get_unit_conversion_nspvar(int i) const{

    static std::vector<realtype> scalor = {
    };
    return scalor[i];
}
};
