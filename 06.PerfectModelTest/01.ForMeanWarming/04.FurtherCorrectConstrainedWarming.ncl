; ====================================================================
; Author: Ziming Chen 
; This code is for:
;   - further correcting the constrained warming over the region 
; where there is a significant correlation between the residual warming trend 
; and the constrained warming 
;   1. indicate the model, exp, projection terms 
;   2. read the regression coef and the sig level, and then take the 
; coef which is significant at the 10% level 
;   3. read the constrained warming of MME and each model, and read 
; the present-day obs warming trend for further correction 
;   4. correct the MME warming at first, and then for each model 
; Attention: This code is for perfect mdoel test. So we should indicate the pseudo-obs model 
; ====================================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code/my_FunctionInPNNL.ncl"
;
undef("my_SelChina")
function my_SelChina(Var:list)
; 本程序挑选中国地区，并mask中国以外的区域
; 注意Var中的矩阵必须包含整个中国，至少是大陆地区，注意最后两位必须是lat x lon
;   1. 读取中国区域的mask场，并将其插值到与Var矩阵同样的维度上
;   2. 挑选包含整个中国的范围，然后基于Var中矩阵的维度，调整mask场的维度
;   3. mask中国以外的区域
local r_Var 
begin
    r_Var             = Var[0]
;
    s_DirRead         = "$WACCEM/my_data/Obs/CN05.1/"
    s_FileRead        = s_DirRead + "CN05.1_Pre_1961_2015_month_025x025.nc"
    FileMask          = addfile(s_FileRead, "r")
    r_Mask            = FileMask->pre(0, :, :)
    
    i_Mask            = where(ismissing(r_Mask) .eq. False, 1, r_Mask)
    copy_VarCoords(r_Mask, i_Mask)
    r_latM            = i_Mask&lat 
    r_lonM            = i_Mask&lon 
    
    i_dims            = dimsizes(r_Var)
    if (dimsizes(i_dims) .eq. 2) then 
        r_VarChina    = r_Var({min(r_latM):max(r_latM)}, \
                              {min(r_lonM):max(r_lonM)})
    else if (dimsizes(i_dims) .eq. 3) then
        r_VarChina    = r_Var(:, {min(r_latM):max(r_latM)}, \
                              {min(r_lonM):max(r_lonM)})
    else if (dimsizes(i_dims) .eq. 4) then
        r_VarChina    = r_Var(:, :, {min(r_latM):max(r_latM)}, \
                              {min(r_lonM):max(r_lonM)})
    end if
    end if
    end if
    ;
    r_latInter       := r_VarChina&lat
    r_lonInter       := r_VarChina&lon
    i_Mask2          := i_Mask({r_latInter}, {r_lonInter}) 
    ;
    if (dimsizes(i_dims) .eq. 2) then 
        i_Mask2_      = i_Mask2
    else if (dimsizes(i_dims) .eq. 3) then
        i_Mask2_      = conform(r_VarChina , i_Mask2, (/1, 2/))
    else if (dimsizes(i_dims) .eq. 4) then
        i_Mask2_      = conform(r_VarChina , i_Mask2, (/2, 3/))
    end if
    end if
    end if
    copy_VarCoords(r_VarChina, i_Mask2_)

    r_VarChina_       = where(i_Mask2_ .ne. 0, r_VarChina, \
                                      r_VarChina@_FillValue)
    copy_VarCoords(r_VarChina, r_VarChina_)
    
    ;; 台湾经纬度范围
    r_lonTaiwan       = (/120, 125/)
    r_latTaiwan       = (/20, 25/)
    r_Temp           := r_VarChina
    r_Temp            = my_mask([/r_Temp/], 0)
    if (dimsizes(i_dims) .eq. 2) then 
        r_VarChina_({r_latTaiwan(0):r_latTaiwan(1)}, \
                        {r_lonTaiwan(0):r_lonTaiwan(1)}) = r_Temp(\
                                            {r_latTaiwan(0):r_latTaiwan(1)}, \
                                            {r_lonTaiwan(0):r_lonTaiwan(1)})
    else if (dimsizes(i_dims) .eq. 3) then
        r_VarChina_(:, {r_latTaiwan(0):r_latTaiwan(1)}, \
                       {r_lonTaiwan(0):r_lonTaiwan(1)}) = r_Temp(:, \
                                            {r_latTaiwan(0):r_latTaiwan(1)}, \
                                            {r_lonTaiwan(0):r_lonTaiwan(1)})
    else if (dimsizes(i_dims) .eq. 4) then
        r_VarChina_(:, :, {r_latTaiwan(0):r_latTaiwan(1)}, \
                          {r_lonTaiwan(0):r_lonTaiwan(1)}) = r_Temp(:, :, \
                                            {r_latTaiwan(0):r_latTaiwan(1)}, \
                                            {r_lonTaiwan(0):r_lonTaiwan(1)})
    end if
    end if
    end if
    ;; 台湾经纬度范围
return(r_VarChina_)
;
end
;
begin
;
; ====================================================================
;   1. indicate the model, exp, projection terms 
;
if (isdefined("s_PseudoModel") .eq. False) then
    s_PseudoModel               = "CanESM5" ; "EC-Earth3"
end if
;
s_4Models                       = "_23Models" ; "_5Models" ; "_30Models" ; "" ; 
; s_ModelSelected                 = (/"ACCESS-CM2", "CNRM-CM6-1", "EC-Earth3", \
;                             "GFDL-CM4", "MPI-ESM1-2-LR"/)
s_ModelSelected                 = (/ "ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                         "BCC-CSM2-MR", "CAMS-CSM1-0", "CanESM5", \
                         "CanESM5-CanOE", "CNRM-CM6-1", "CNRM-ESM2-1", \
                         "EC-Earth3", "EC-Earth3-Veg", "FGOALS-f3-L", \
                         "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
                         "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8", \
                         "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR", \
                         "MRI-ESM2-0", "UKESM1-0-LL"/)
s_Month                         = "AnMean"
if (isdefined("s_Exp") .eq. False) then
    s_Exp                       = "ssp245" ;
end if
s_Terms                         = (/"Near", "Mid", "Long"/)
;
i_yrConstrained                 = (/1981, 2014/)
; ====================================================================
;   2. read the regression coef and the sig level, and then take the 
; coef which is significant at the 10% level 
s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/06.PerfectModelTest/01.ForMeanWarming/"
s_FileRead                      = s_DirRead + \
                            "Regression_LocalTrendAndProjectedWarming_PseModel_"\
                            + s_PseudoModel + "_" + s_Exp + s_4Models + ".nc"
FileReg                         = addfile(s_FileRead, "r")
;
s_VarName2                      = "RegInTimDim_" + s_Exp
r_RegSigLev                     = FileReg->$s_VarName2$(:, :, {-90:90}, :)
r_RegSigLev                     = my_LonTest([/r_RegSigLev/], 360)
;
r_Reg                           = r_RegSigLev(0, :, :, :)
r_Reg                           = where(r_RegSigLev(1, :, :, :) .le. 0.1, r_Reg,\
                                        r_Reg@_FillValue)
;
s_DirRead = "$WACCEM/my_data/CMIP6/tas/POST/ProjectionOfAllReal/"
s_FileRead                      = s_DirRead + "ProjectionInDifTerms_AllReal" + \
                                s_4Models + ".nc"
FileProjChanges                 = addfile(s_FileRead, "r")
;
s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/06.PerfectModelTest/01.ForMeanWarming/"
s_FileRead                      = s_DirRead + \
                "ConstrainedWarmingInChina_RegMMEForScaling_PseModel_" + \
                s_PseudoModel + s_4Models + ".nc"
FileConWarming                  = addfile(s_FileRead, "r")
s_VarName2                      = "ConstrainedWarmingChina_MME_" + s_Exp
r_ConWarmingMME                 = FileConWarming->$s_VarName2$
r_ConWarmingMME_Cor             = r_ConWarmingMME
;
s_DirRead = "$WACCEM/my_data/CMIP6/tas/POST/GSAT_AndRegionalTrend_AllReal/"
s_FileRead                      = s_DirRead + "GSAT_AndRegionalTrend_" + \
            i_yrConstrained(0) + "-" + i_yrConstrained(1) + s_4Models + ".nc"
FileReTrend                     = addfile(s_FileRead, "r")
s_VarName2                      = "ResidualTrend_historical"
r_HistReTrend                   = FileReTrend->$s_VarName2$(1, :, :)
printMinMax(r_HistReTrend, False)
;
l_list                         := my_CMIP6Model_Exp_Path_ForAllReal(\
                                s_PseudoModel, s_Exp, "Amon")
; s_ExpModel_Path                := l_list[0]
s_Real                         := l_list[1]
l_list_                        := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                s_PseudoModel, s_Exp, "Amon")
; s_ExpModel_Path_               := l_list_[0]
s_Real_                        := l_list_[1]
if (all(ismissing(s_Real)) .or. dimsizes(s_Real_) .gt. dimsizes(s_Real) ) then
        s_Real                 := s_Real_
        ; s_ExpModel_Path        := s_ExpModel_Path_
end if
do iReal = 0, dimsizes(s_Real)-1
    s_VarName2                  = "ResidualTrend_" + s_PseudoModel + \
                                "_historical_" + s_Real(iReal)
    if (my_VarNameExist(FileReTrend, s_VarName2) .eq. False) then
        print(s_VarName2 + " Missing!")
        continue
    end if
    r_Temp                     := FileReTrend->$s_VarName2$(1, {-90:90}, :)
    r_Temp                     := my_LonTest([/r_Temp/], 360)
    ;
    s_VarName2                  = "tas_inDifTerms_" + s_PseudoModel + "_" + \
                                s_Exp + "_" + s_Real(iReal)
    if (my_VarNameExist(FileProjChanges, s_VarName2) .eq. False) then
        print(s_VarName2 + " Missing!")
        continue
    end if
    r_PseWarming_Temp          := FileProjChanges->$s_VarName2$(:, {-90:90}, :)
    r_PseWarming_Temp          := my_LonTest([/r_PseWarming_Temp/], 360)
    ;
    if (isdefined("r_TrendObsMean") .eq. False) then
        i_dims                 := array_append_record(dimsizes(s_Real), \
                                                      dimsizes(r_Temp), 0)
        r_TrendObsMean          = new(i_dims, float)
        copy_VarCoords(r_Temp, r_TrendObsMean(0, :, :))
        r_TrendObsMean!0        = "Real"
        r_TrendObsMean&Real     = s_Real
        ;
        i_dims                 := array_append_record(dimsizes(s_Real), \
                                                dimsizes(r_PseWarming_Temp), 0)
        r_PseWarming           := new(i_dims, float)
        copy_VarCoords(r_PseWarming_Temp, r_PseWarming(0, :, :, :))
        copy_VarCoords(r_TrendObsMean(:, 0, 0), r_PseWarming(:, 0, 0, 0))
    end if
    r_TrendObsMean(iReal, :, :) = r_Temp
    r_PseWarming(iReal, :, :, :)= r_PseWarming_Temp
end do
r_TrendObsMean                 := dim_avg_n_Wrap(r_TrendObsMean, 0)
r_PseWarming                   := dim_avg_n_Wrap(r_PseWarming, 0)
print("Read the pseudo-obs residual trend!")
print("")
;
; ====================================================================
;   4. correct the MME warming at first, and then for each model 
;; constrain MME 
do iPeriod = 0, dimsizes(s_Terms)-1
    r_HistReTrend_Temp         := where(ismissing(r_Reg(iPeriod, :, :)), \
                                        r_HistReTrend@_FillValue, r_HistReTrend)
    copy_VarCoords(r_HistReTrend, r_HistReTrend_Temp)
    r_TrendObsMean_Temp        := where(ismissing(r_HistReTrend_Temp), \
                                    r_TrendObsMean@_FillValue, r_TrendObsMean)
    copy_VarCoords(r_TrendObsMean, r_TrendObsMean_Temp)
    ;
    dx                         := r_HistReTrend_Temp
    dx                          = r_HistReTrend_Temp - r_TrendObsMean_Temp
    ;
    dy                          = dx * r_Reg(iPeriod, :, :)
    ;
    r_Temp                     := r_ConWarmingMME(iPeriod, :, :) 
    r_Temp                      = where(ismissing(r_HistReTrend_Temp), \
                                        r_Temp@_FillValue, r_Temp)
    r_Temp                      = r_Temp - dy 
    ;
    r_ConWarmingMME_Cor(iPeriod, :, :) = where(ismissing(r_Temp) .eq. False, \
                                    r_Temp, r_ConWarmingMME_Cor(iPeriod, :, :))
end do
;
do iModel = 0, dimsizes(s_ModelSelected)-1
    if (s_ModelSelected(iModel) .eq. s_PseudoModel) then
        print(s_ModelSelected(iModel) + " is pseudo-obs model!")
        print("Skip!")
        continue
    end if
    ;
    l_list                     := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp, "Amon")
    s_ExpModel_Path            := l_list[0]
    s_Real                     := l_list[1]
    l_list_                    := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp, "Amon")
    s_ExpModel_Path_           := l_list_[0]
    s_Real_                    := l_list_[1]
    if (all(ismissing(s_Real)) .or. dimsizes(s_Real_).gt.dimsizes(s_Real)) then
        s_Real                 := s_Real_
        s_ExpModel_Path        := s_ExpModel_Path_
    end if
    ;
    do iReal = 0, dimsizes(s_Real)-1
            ;
        s_VarNameProj          := "tas_inDifTerms_" + s_ModelSelected(iModel) + \
                                "_" + s_Exp + "_" + s_Real(iReal)
        if (my_VarNameExist(FileProjChanges, s_VarNameProj) .eq. False) then
            print(s_VarNameProj + "  Missing!")
            continue
        end if
        r_ProjTerm             := FileProjChanges->$s_VarNameProj$(:,{-90:90},:)
        r_ProjTerm             := my_LonTest([/r_ProjTerm/], 360)
        ;
        ;; store up 
        if (isdefined("r_ProjTerm_AllReal").eq.False .or. iReal.eq.0) then
            i_dims                     := array_append_record(\
                                    dimsizes(s_Real), dimsizes(r_ProjTerm), 0)
            r_ProjTerm_AllReal         := new(i_dims, float)
            copy_VarCoords(r_ProjTerm, r_ProjTerm_AllReal(0, :, :, :))
            r_ProjTerm_AllReal!0        = "Real"
            r_ProjTerm_AllReal&Real     = s_Real
        end if
        r_ProjTerm_AllReal(iReal,:,:,:) = r_ProjTerm
    end do
    ;
    if (isdefined("r_ProjTerm_AllModel") .eq. False) then
        i_dims                         := array_append_record(\
                                    dimsizes(s_ModelSelected), \
                                    dimsizes(r_ProjTerm_AllReal(0, :, :, :)), 0)
        r_ProjTerm_AllModel            := new(i_dims, float)
        copy_VarCoords(r_ProjTerm_AllReal(0,:,:,:), r_ProjTerm_AllModel(0,:,:,:))
        r_ProjTerm_AllModel!0           = "Model"
        r_ProjTerm_AllModel&Model       = s_ModelSelected
    end if
    r_ProjTerm_AllModel(iModel, :, :, :)= dim_avg_n_Wrap(r_ProjTerm_AllReal, 0)
end do
r_ProjTerm_AllModel    := dim_avg_n_Wrap(r_ProjTerm_AllModel, 0)
;
;;; Area-Mean over China 
r_China                := my_SelChina([/r_ConWarmingMME_Cor/])
;
r_China_Raw            := my_SelChina([/r_ProjTerm_AllModel/])
;
r_China_Pseudo         := my_SelChina([/r_PseWarming/])
;
;; PCC & RMSE
r_PCC                  := new(3, float)
copy_VarCoords(r_ProjTerm_AllModel(:, 0, 0), r_PCC)
r_PCC_Raw               = r_PCC
r_RMSE                  = r_PCC
r_RMSE_Raw              = r_PCC
do iTerm = 0, dimsizes(r_PCC)-1
    l_list             := my_PCC_RMSE_for2D0(r_China_Pseudo(iTerm, :, :), \
                                             r_China(iTerm, :, :), 0)
    r_PCC(iTerm)        = l_list[0]
    r_RMSE(iTerm)       = l_list[1]
    l_list             := my_PCC_RMSE_for2D0(r_China_Pseudo(iTerm, :, :), \
                                             r_China_Raw(iTerm, :, :), 0)
    r_PCC_Raw(iTerm)    = l_list[0]
    r_RMSE_Raw(iTerm)   = l_list[1]
end do
print("RMSE in 3 Term: Constrained, Raw")
print(r_RMSE + ", " + r_RMSE_Raw)
print("")
;
print("3 Term Warming: Constrained, Raw, Pseudo")
r_China                := my_AreaMean2([/r_China/], r_China&lat, r_China&lon, 0)
r_China_Raw            := my_AreaMean2([/r_China_Raw/], r_China_Raw&lat, \
                                       r_China_Raw&lon, 0)
r_China_Pseudo         := my_AreaMean2([/r_China_Pseudo/], r_China_Pseudo&lat, \
                                       r_China_Pseudo&lon, 0)
print(r_China + " K, " + r_China_Raw + " K, " + r_China_Pseudo)
print("")
; exit()
; CanESM5
; (0) RMSE in 3 Term: Constrained, Raw
; (0) 0.354402, 0.656835
; (1) 0.451973, 0.988818
; (2) 0.594326, 1.36425
; (0) 
; (0) 3 Term Warming: Constrained, Raw, Pseudo
; (0) 1.59462 K, 1.07743 K, 1.62685
; (1) 2.67818 K, 1.87176 K, 2.75225
; (2) 3.97032 K, 2.94043 K, 4.21124
;
; EC-Earth3
; (0) RMSE in 3 Term: Constrained, Raw
; (0) 0.260411, 0.146651
; (1) 0.612657, 0.196653
; (2) 1.10971, 0.439666
; (0) 
; (0) 3 Term Warming: Constrained, Raw, Pseudo
; (0) 0.807184 K, 1.10588 K, 1.02941
; (1) 1.40217 K, 1.90907 K, 1.96869
; (2) 2.2849 K, 2.9814 K, 3.351
;
s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/06.PerfectModelTest/01.ForMeanWarming/"
system("mkdir -p " + s_DirOut)
s_FileOut               = s_DirOut + "FurtherCorrectConstrainedWarming_PseModel_" + \
                        s_PseudoModel + "_" + s_Exp + s_4Models + ".nc"
ls                      = systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileOut             = addfile(s_FileOut, "c")
    GAttr               = True
    GAttr@history       = systemfunc("date") + "  " + systemfunc("pwd") + \
                        "/04.FurtherCorrectConstrainedWarming.ncl"
    fileattdef(FileOut, GAttr)
    ;
else 
    FileOut             = addfile(s_FileOut, "w")
end if
s_VarName2              = "ConstrainedWarming_" + s_Exp
r_ConWarmingMME_Cor@long_name = s_VarName2
r_ConWarmingMME_Cor@CSAT= r_China
FileOut->$s_VarName2$   = r_ConWarmingMME_Cor
;
s_VarName2              = "PseudoWarming_" + s_Exp
r_PseWarming@long_name  = s_VarName2 + " of " + s_PseudoModel
r_PseWarming@CSAT       = r_China_Pseudo
FileOut->$s_VarName2$   = r_PseWarming
;
s_VarName2              = "RawProjectedWarming_" + s_Exp
r_ProjTerm_AllModel@long_name = s_VarName2
r_ProjTerm_AllModel@CSAT= r_China_Raw
FileOut->$s_VarName2$   = r_ProjTerm_AllModel
;
s_VarName2              = "RMSE_China_ConstrainedWarming_" + s_Exp
r_RMSE@long_name        = s_VarName2
r_RMSE@Value            = r_RMSE
FileOut->$s_VarName2$   = r_RMSE
;
s_VarName2              = "RMSE_China_RawProjectedWarming_" + s_Exp
r_RMSE_Raw@long_name    = s_VarName2
r_RMSE_Raw@Value        = r_RMSE_Raw
FileOut->$s_VarName2$   = r_RMSE_Raw
;
;; for each CMIP6
do iModel = 0, dimsizes(s_ModelSelected)-1
    if (s_ModelSelected(iModel) .eq. s_PseudoModel) then
        print(s_ModelSelected(iModel) + " is pseudo-obs model!")
        print("Skip!")
        continue
    end if
    s_VarName2          = "ConstrainedWarmingChina_" + s_ModelSelected(iModel)+\
                        s_Exp
    if (my_VarNameExist(FileConWarming, s_VarName2) .eq. False) then
        print(s_VarName2 + " Missing!")
        continue
    end if
    r_ConWarming       := FileConWarming->$s_VarName2$
    r_ConWarming_Cor   := r_ConWarming
    ;
    s_VarName2          = "ResidualTrend_" + s_ModelSelected(iModel) + \
                        "_historical"
    r_HistReTrend      := FileReTrend->$s_VarName2$(1, :, :)
    ;
    do iPeriod = 0, dimsizes(s_Terms)-1
        r_HistReTrend_Temp := where(ismissing(r_Reg(iPeriod, :, :)), \
                                        r_HistReTrend@_FillValue, r_HistReTrend)
        copy_VarCoords(r_HistReTrend, r_HistReTrend_Temp)
        r_TrendObsMean_Temp:= where(ismissing(r_HistReTrend_Temp), \
                                    r_TrendObsMean@_FillValue, r_TrendObsMean)
        copy_VarCoords(r_TrendObsMean, r_TrendObsMean_Temp)
    ;
        dx                 := r_HistReTrend_Temp
        dx                  = r_HistReTrend_Temp - r_TrendObsMean_Temp
    ;
        dy                  = dx * r_Reg(iPeriod, :, :)
    ;
        r_Temp             := r_ConWarming(iPeriod, :, :) 
        r_Temp              = where(ismissing(r_HistReTrend_Temp), \
                                        r_Temp@_FillValue, r_Temp)
        r_Temp              = r_Temp - dy 
    ;
        r_ConWarming_Cor(iPeriod, :, :) = where(ismissing(r_Temp) .eq. False, \
                                    r_Temp, r_ConWarming_Cor(iPeriod, :, :))
    end do
    ;
    s_VarName2             := "ConstrainedWarming_" + s_ModelSelected(iModel) + \
                            "_" + s_Exp
    FileOut->$s_VarName2$   = r_ConWarming_Cor
    print(s_VarName2 + " Done!")
    ;
    ;;; Area-Mean over China 
    r_China                := my_SelChina([/r_ConWarming_Cor/])
    r_China                := my_AreaMean2([/r_China/], r_China&lat, r_China&lon, 0)
    ; r_China_Raw            := my_SelChina([/r_ConWarming/])
    ; r_China_Raw            := my_AreaMean2([/r_China_Raw/], r_China_Raw&lat, \
    ;                                    r_China_Raw&lon, 0)
    ; print("3 Term Warming: Constrained, Raw, Pseudo")
    ; print(r_China + " K, " + r_China_Raw + " K, " + r_China_Pseudo)
    print("")
end do
print(s_FileOut + "")

end





















