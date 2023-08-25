#!/bin/bash
# this code control the model input of 01 to 04
# so that we can run all the models sequentially
# usage: bash 0.ControlBash.sh ssp245


# s_Exp=("ssp370" "ssp585") #
s_Exp=$1
# 
echo $s_Exp

Dir="/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Code/07.UncerQualification/02.ForExtremeWarming/"
for iExp in ${s_Exp[*]}
do
		echo $iExp
		# 01.Dis_of_pc_and_constrainedGSAT
		FileRun="01.Dis_of_pc_and_constrainedGSAT.ncl"
		echo $FileRun
		ncl s_Exp=\"${iExp}\" ${Dir}$FileRun
		echo ""
		#
		# 02.Dis_of_pc1_and_constrainedLocalWarming_Yc1.ncl
		FileRun="02.Dis_of_pc1_and_constrainedLocalWarming_Yc1.ncl"
		echo $FileRun
		nohup ncl s_Exp=\"${iExp}\" ${Dir}$FileRun &> 02.${iExp}.out 
		tail -100 02.${iExp}.out 
		rm -rf 02.${iExp}.out 
		echo ""
		#
		# 03.Dis_of_pc2_and_constrainedLocalWarmingByLocalTrend_Yc2.ncl
		FileRun="03.Dis_of_pc2_and_constrainedLocalWarmingByLocalTrend_Yc2.ncl"
		echo $FileRun
		ncl s_Exp=\"${iExp}\" ${Dir}$FileRun
		echo ""
		#
		# 01.CorGSAT_AndLocWarmingByScaling_ForEachReal.ncl
		FileRun="04.ConEachRealForInterVarEst/01.CorGSAT_AndLocWarmingByScaling_ForEachReal.ncl"
		echo $FileRun
		ncl s_Exp=\"${iExp}\" ${Dir}$FileRun
		#
		# 02.FurtherCorrectConstrainedWarming_ForEachReal.ncl
		FileRun="04.ConEachRealForInterVarEst/02.FurtherCorrectConstrainedWarming_ForEachReal.ncl"
		echo $FileRun
		ncl s_Exp=\"${iExp}\" ${Dir}$FileRun
		#
		# 03.EstInterVar_fromEachLENS_1000times.ncl
		FileRun="04.ConEachRealForInterVarEst/03.EstInterVar_fromEachLENS_1000times.ncl"
		echo $FileRun
		ncl s_Exp=\"${iExp}\" ${Dir}$FileRun
	echo ""
done














