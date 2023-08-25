#!/bin/bash
# this code control the model input of 01 to 04
# so that we can run all the models sequentially

s_Model=("ACCESS-CM2" "ACCESS-ESM1-5" "AWI-CM-1-1-MR" \
         "BCC-CSM2-MR" "CAMS-CSM1-0" "CNRM-CM6-1" \
         "CNRM-ESM2-1" "CanESM5" "EC-Earth3" \
         "EC-Earth3-Veg" "FGOALS-g3" "GFDL-CM4" \
         "GFDL-ESM4" "IPSL-CM6A-LR" "INM-CM4-8" \
         "MIROC6" "MRI-ESM2-0" "UKESM1-0-LL")

s_Exp=("ssp126" "ssp245" "ssp370" "ssp585") #

Dir="/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Code/06.PerfectModelTest/02.ForExtremeWarming/"
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
		FileRun="02.CorTXxWarmingByScaling_withDif_RegMMEForScaling_ForPMT.ncl"
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

rm -rf  FurtherCorrectConstrainedWarming_PseModel_ACCESS-ESM1-5_ssp126_17Models.nc Regression_LocalTrendAndProjectedWarming_PseModel_ACCESS-ESM1-5_ssp126_17Models.nc