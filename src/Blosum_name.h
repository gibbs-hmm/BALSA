/****************************************************************/
/* BALSA - Bayesian Algorithm for Local Sequence Alignment      */
/*                                                              */
/* Please acknowledge the program authors on any publication of */
/* scientific results based in part on use of the program and   */
/* cite the following articles in which the program was         */
/* described.                                                   */
/*                                                              */
/* Webb-Robertson, B. J., L. A. McCue, et al. (2008). "Measuring*/
/* Global Credibility with Application to Local Sequence        */
/* Alignment." PLoS CompBiol. 4(5): e1000077.                   */
/* Webb, B. J., J. S. Liu, et al. (2002). "BALSA: Bayesian      */
/* algorithm for local sequence alignment." Nucleic Acids Res   */
/* 30(5): 1268-77.                                              */
/*                                                              */
/*                                                              */
/* Copyright (C) 2006   Health Research Inc.                    */
/* HEALTH RESEARCH INCORPORATED (HRI),                          */
/* ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455.             */
/* Email:  gibbsamp@wadsworth.org                              */
/*                                                              */
/*                                                              */
/* Copyright (C) 2009   Brown University                        */
/* Brown University                                             */
/* Providence, RI 02912                                         */
/* Email:  gibbs@brown.edu                                     */
/****************************************************************/
/*                                                              */
/* This program is free software; you can redistribute it       */
/* and/or modify it under the terms of the GNU General Public   */
/* License as published by the Free Software Foundation;        */
/* either version 2 of the License, or (at your option)         */
/* any later version.                                           */
/*                                                              */
/* This program is distributed in the hope that it will be      */
/* useful, but WITHOUT ANY WARRANTY; without even the implied   */
/* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR      */
/* PURPOSE. See the GNU General Public License for more         */
/* details.                                                     */
/*                                                              */
/* You should have received a copy of the GNU General Public    */
/* License along with this program; if not, write to the        */
/* Free Software Foundation, Inc., 51 Franklin Street,          */
/* Fifth Floor, Boston, MA  02110-1301, USA.                    */
/****************************************************************/

//#ifndef LETTER_H
//#define LETTER_H


static char const *name_related[]={"BLOSUM_30","BLOSUM_40", "BLOSUM_45", "BLOSUM_50","BLOSUM_55","BLOSUM_62", 
			     "BLOSUM_70","BLOSUM_80", "PAM250",
"PAM1_DNA", "PAM2_DNA", "PAM3_DNA", "PAM4_DNA", "PAM5_DNA", "PAM6_DNA", "PAM7_DNA", "PAM8_DNA", "PAM9_DNA", "PAM10_DNA", 
"PAM11_DNA", "PAM12_DNA", "PAM13_DNA", "PAM14_DNA", "PAM15_DNA", "PAM16_DNA", "PAM17_DNA", "PAM18_DNA", "PAM19_DNA", "PAM20_DNA", 
"PAM21_DNA", "PAM22_DNA", "PAM23_DNA", "PAM24_DNA", "PAM25_DNA", "PAM26_DNA", "PAM27_DNA", "PAM28_DNA", "PAM29_DNA", "PAM30_DNA", 
"PAM31_DNA", "PAM32_DNA", "PAM33_DNA", "PAM34_DNA", "PAM35_DNA", "PAM36_DNA", "PAM37_DNA", "PAM38_DNA", "PAM39_DNA", "PAM40_DNA", 
"PAM41_DNA", "PAM42_DNA", "PAM43_DNA", "PAM44_DNA", "PAM45_DNA", "PAM46_DNA", "PAM47_DNA", "PAM48_DNA", "PAM49_DNA", "PAM50_DNA", 
"PAM51_DNA", "PAM52_DNA", "PAM53_DNA", "PAM54_DNA", "PAM55_DNA", "PAM56_DNA", "PAM57_DNA", "PAM58_DNA", "PAM59_DNA", "PAM60_DNA", 
"PAM61_DNA", "PAM62_DNA", "PAM63_DNA", "PAM64_DNA", "PAM65_DNA", "PAM66_DNA", "PAM67_DNA", "PAM68_DNA", "PAM69_DNA", "PAM70_DNA", 
"PAM71_DNA", "PAM72_DNA", "PAM73_DNA", "PAM74_DNA", "PAM75_DNA", "PAM76_DNA", "PAM77_DNA", "PAM78_DNA", "PAM79_DNA", "PAM80_DNA", 
"PAM81_DNA", "PAM82_DNA", "PAM83_DNA", "PAM84_DNA", "PAM85_DNA", "PAM86_DNA", "PAM87_DNA", "PAM88_DNA", "PAM89_DNA", "PAM90_DNA", 
"PAM91_DNA", "PAM92_DNA", "PAM93_DNA", "PAM94_DNA", "PAM95_DNA", "PAM96_DNA", "PAM97_DNA", "PAM98_DNA", "PAM99_DNA", "PAM100_DNA", 
"PAM101_DNA", "PAM102_DNA", "PAM103_DNA", "PAM104_DNA", "PAM105_DNA", "PAM106_DNA", "PAM107_DNA", "PAM108_DNA", "PAM109_DNA", "PAM110_DNA", 
"PAM111_DNA", "PAM112_DNA", "PAM113_DNA", "PAM114_DNA", "PAM115_DNA", "PAM116_DNA", "PAM117_DNA", "PAM118_DNA", "PAM119_DNA", "PAM120_DNA", 
"PAM121_DNA", "PAM122_DNA", "PAM123_DNA", "PAM124_DNA", "PAM125_DNA", "PAM126_DNA", "PAM127_DNA", "PAM128_DNA", "PAM129_DNA", "PAM130_DNA", 
"PAM131_DNA", "PAM132_DNA", "PAM133_DNA", "PAM134_DNA", "PAM135_DNA", "PAM136_DNA", "PAM137_DNA", "PAM138_DNA", "PAM139_DNA", "PAM140_DNA", 
"PAM141_DNA", "PAM142_DNA", "PAM143_DNA", "PAM144_DNA", "PAM145_DNA", "PAM146_DNA", "PAM147_DNA", "PAM148_DNA", "PAM149_DNA", "PAM150_DNA", 
"PAM151_DNA", "PAM152_DNA", "PAM153_DNA", "PAM154_DNA", "PAM155_DNA", "PAM156_DNA", "PAM157_DNA", "PAM158_DNA", "PAM159_DNA", "PAM160_DNA", 
"PAM161_DNA", "PAM162_DNA", "PAM163_DNA", "PAM164_DNA", "PAM165_DNA", "PAM166_DNA", "PAM167_DNA", "PAM168_DNA", "PAM169_DNA", "PAM170_DNA", 
"PAM171_DNA", "PAM172_DNA", "PAM173_DNA", "PAM174_DNA", "PAM175_DNA", "PAM176_DNA", "PAM177_DNA", "PAM178_DNA", "PAM179_DNA", "PAM180_DNA", 
"PAM181_DNA", "PAM182_DNA", "PAM183_DNA", "PAM184_DNA", "PAM185_DNA", "PAM186_DNA", "PAM187_DNA", "PAM188_DNA", "PAM189_DNA", "PAM190_DNA", 
"PAM191_DNA", "PAM192_DNA", "PAM193_DNA", "PAM194_DNA", "PAM195_DNA", "PAM196_DNA", "PAM197_DNA", "PAM198_DNA", "PAM199_DNA", "PAM200_DNA", 
"PAM201_DNA", "PAM202_DNA", "PAM203_DNA", "PAM204_DNA", "PAM205_DNA", "PAM206_DNA", "PAM207_DNA", "PAM208_DNA", "PAM209_DNA", "PAM210_DNA", 
"PAM211_DNA", "PAM212_DNA", "PAM213_DNA", "PAM214_DNA", "PAM215_DNA", "PAM216_DNA", "PAM217_DNA", "PAM218_DNA", "PAM219_DNA", "PAM220_DNA", 
"PAM221_DNA", "PAM222_DNA", "PAM223_DNA", "PAM224_DNA", "PAM225_DNA", "PAM226_DNA", "PAM227_DNA", "PAM228_DNA", "PAM229_DNA", "PAM230_DNA", 
"PAM231_DNA", "PAM232_DNA", "PAM233_DNA", "PAM234_DNA", "PAM235_DNA", "PAM236_DNA", "PAM237_DNA", "PAM238_DNA", "PAM239_DNA", "PAM240_DNA", 
"PAM241_DNA", "PAM242_DNA", "PAM243_DNA", "PAM244_DNA", "PAM245_DNA", "PAM246_DNA", "PAM247_DNA", "PAM248_DNA", "PAM249_DNA", "PAM250_DNA", 
"PAM251_DNA", "PAM252_DNA", "PAM253_DNA", "PAM254_DNA", "PAM255_DNA", "PAM256_DNA", "PAM257_DNA", "PAM258_DNA", "PAM259_DNA", "PAM260_DNA", 
"PAM261_DNA", "PAM262_DNA", "PAM263_DNA", "PAM264_DNA", "PAM265_DNA", "PAM266_DNA", "PAM267_DNA", "PAM268_DNA", "PAM269_DNA", "PAM270_DNA", 
"PAM271_DNA", "PAM272_DNA", "PAM273_DNA", "PAM274_DNA", "PAM275_DNA", "PAM276_DNA", "PAM277_DNA", "PAM278_DNA", "PAM279_DNA", "PAM280_DNA", 
"PAM281_DNA", "PAM282_DNA", "PAM283_DNA", "PAM284_DNA", "PAM285_DNA", "PAM286_DNA", "PAM287_DNA", "PAM288_DNA", "PAM289_DNA", "PAM290_DNA", 
"PAM291_DNA", "PAM292_DNA", "PAM293_DNA", "PAM294_DNA", "PAM295_DNA", "PAM296_DNA", "PAM297_DNA", "PAM298_DNA", "PAM299_DNA", "PAM300_DNA", 
"PAM301_DNA", "PAM302_DNA", "PAM303_DNA", "PAM304_DNA", "PAM305_DNA", "PAM306_DNA", "PAM307_DNA", "PAM308_DNA", "PAM309_DNA", "PAM310_DNA", 
"PAM311_DNA", "PAM312_DNA", "PAM313_DNA", "PAM314_DNA", "PAM315_DNA", "PAM316_DNA", "PAM317_DNA", "PAM318_DNA", "PAM319_DNA", "PAM320_DNA", 
"PAM321_DNA", "PAM322_DNA", "PAM323_DNA", "PAM324_DNA", "PAM325_DNA", "PAM326_DNA", "PAM327_DNA", "PAM328_DNA", "PAM329_DNA", "PAM330_DNA", 
"PAM331_DNA", "PAM332_DNA", "PAM333_DNA", "PAM334_DNA", "PAM335_DNA", "PAM336_DNA", "PAM337_DNA", "PAM338_DNA", "PAM339_DNA", "PAM340_DNA", 
"PAM341_DNA", "PAM342_DNA", "PAM343_DNA", "PAM344_DNA", "PAM345_DNA", "PAM346_DNA", "PAM347_DNA", "PAM348_DNA", "PAM349_DNA", "PAM350_DNA", 
"PAM351_DNA", "PAM352_DNA", "PAM353_DNA", "PAM354_DNA", "PAM355_DNA", "PAM356_DNA", "PAM357_DNA", "PAM358_DNA", "PAM359_DNA", "PAM360_DNA", 
"PAM361_DNA", "PAM362_DNA", "PAM363_DNA", "PAM364_DNA", "PAM365_DNA", "PAM366_DNA", "PAM367_DNA", "PAM368_DNA", "PAM369_DNA", "PAM370_DNA", 
"PAM371_DNA", "PAM372_DNA", "PAM373_DNA", "PAM374_DNA", "PAM375_DNA", "PAM376_DNA", "PAM377_DNA", "PAM378_DNA", "PAM379_DNA", "PAM380_DNA", 
"PAM381_DNA", "PAM382_DNA", "PAM383_DNA", "PAM384_DNA", "PAM385_DNA", "PAM386_DNA", "PAM387_DNA", "PAM388_DNA", "PAM389_DNA", "PAM390_DNA", 
"PAM391_DNA", "PAM392_DNA", "PAM393_DNA", "PAM394_DNA", "PAM395_DNA", "PAM396_DNA", "PAM397_DNA", "PAM398_DNA", "PAM399_DNA", "PAM400_DNA", 
"PAM401_DNA", "PAM402_DNA", "PAM403_DNA", "PAM404_DNA", "PAM405_DNA", "PAM406_DNA", "PAM407_DNA", "PAM408_DNA", "PAM409_DNA", "PAM410_DNA", 
"PAM411_DNA", "PAM412_DNA", "PAM413_DNA", "PAM414_DNA", "PAM415_DNA", "PAM416_DNA", "PAM417_DNA", "PAM418_DNA", "PAM419_DNA", "PAM420_DNA", 
"PAM421_DNA", "PAM422_DNA", "PAM423_DNA", "PAM424_DNA", "PAM425_DNA", "PAM426_DNA", "PAM427_DNA", "PAM428_DNA", "PAM429_DNA", "PAM430_DNA", 
"PAM431_DNA", "PAM432_DNA", "PAM433_DNA", "PAM434_DNA", "PAM435_DNA", "PAM436_DNA", "PAM437_DNA", "PAM438_DNA", "PAM439_DNA", "PAM440_DNA", 
"PAM441_DNA", "PAM442_DNA", "PAM443_DNA", "PAM444_DNA", "PAM445_DNA", "PAM446_DNA", "PAM447_DNA", "PAM448_DNA", "PAM449_DNA", "PAM450_DNA", 
"PAM451_DNA", "PAM452_DNA", "PAM453_DNA", "PAM454_DNA", "PAM455_DNA", "PAM456_DNA", "PAM457_DNA", "PAM458_DNA", "PAM459_DNA", "PAM460_DNA", 
"PAM461_DNA", "PAM462_DNA", "PAM463_DNA", "PAM464_DNA", "PAM465_DNA", "PAM466_DNA", "PAM467_DNA", "PAM468_DNA", "PAM469_DNA", "PAM470_DNA", 
"PAM471_DNA", "PAM472_DNA", "PAM473_DNA", "PAM474_DNA", "PAM475_DNA", "PAM476_DNA", "PAM477_DNA", "PAM478_DNA", "PAM479_DNA", "PAM480_DNA", 
"PAM481_DNA", "PAM482_DNA", "PAM483_DNA", "PAM484_DNA", "PAM485_DNA", "PAM486_DNA", "PAM487_DNA", "PAM488_DNA", "PAM489_DNA", "PAM490_DNA", 
"PAM491_DNA", "PAM492_DNA", "PAM493_DNA", "PAM494_DNA", "PAM495_DNA", "PAM496_DNA", "PAM497_DNA", "PAM498_DNA", "PAM499_DNA", "PAM500_DNA", 
"PAM_DNA"};


//#endif

