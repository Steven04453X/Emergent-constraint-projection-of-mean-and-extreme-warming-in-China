#!/bin/bash
# this code control the model input of 01 to 04
# so that we can run all the models sequentially

s_Model=("ACCESS-CM2" "ACCESS-ESM1-5" "AWI-CM-1-1-MR" \
          "BCC-CSM2-MR" "CAMS-CSM1-0" "CanESM5" \
          "CanESM5-CanOE" "CNRM-CM6-1" "CNRM-ESM2-1" \
          "EC-Earth3" "EC-Earth3-Veg" "FGOALS-f3-L" \
          "FGOALS-g3" "GFDL-CM4" "GFDL-ESM4" \
          "GISS-E2-1-G" "HadGEM3-GC31-LL" "INM-CM4-8" \
          "IPSL-CM6A-LR" "MIROC6" "MPI-ESM1-2-LR" \
          "MRI-ESM2-0" "UKESM1-0-LL")

s_Exp=("ssp126" "ssp245" "ssp370" "ssp585") #

Dir="/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Code/06.PerfectModelTest/01.ForMeanWarming/"
for iModel in ${s_Model[*]}
do
	echo $iModel
	for iExp in ${s_Exp[*]}
	do
		echo $iExp
		# 01.ConstrainGSAT_AllReal_ForPMT.ncl
		FileRun="01.ConstrainGSAT_AllReal_ForPMT.ncl"
		echo $FileRun
		ncl s_PseudoModel=\"${iModel}\" s_Exp=\"${iExp}\" ${Dir}$FileRun
		echo ""
		#
		# 02.ScaledWarmingInChina_withDif_RegMMEForScaling_ForPMT.ncl
		FileRun="02.ScaledWarmingInChina_withDif_RegMMEForScaling_ForPMT.ncl"
		echo $FileRun
		ncl s_PseudoModel=\"${iModel}\" s_Exp=\"${iExp}\" ${Dir}$FileRun
		echo ""
		#
		# 03.LocalTrendAndProjectedWarming_ForPMT.ncl
		FileRun="03.LocalTrendAndProjectedWarming_ForPMT.ncl"
		echo $FileRun
		ncl s_PseudoModel=\"${iModel}\" s_Exp=\"${iExp}\" ${Dir}$FileRun
		echo ""
		#
		FileRun="04.FurtherCorrectConstrainedWarming.ncl"
		echo $FileRun
		ncl s_PseudoModel=\"${iModel}\" s_Exp=\"${iExp}\" ${Dir}$FileRun
		echo ""
	done
	echo ""
done