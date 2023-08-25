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
s_Exp                           = "ssp245" ;"ssp126" ; "ssp370" ; "ssp585" ; 
s_Terms                         = (/"Near", "Mid", "Long"/)
;
i_yrConstrained                 = (/1981, 2014/)
; ====================================================================
;   2. read the regression coef and the sig level, and then take the 
; coef which is significant at the 10% level 
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/03.ConstraintByUsingLocalTrend/"
s_FileRead                      = s_DirRead + \
                            "Regression_LocalTrendAndProjectedWarming_" + s_Exp\
                            + s_4Models + ".nc"
FileReg                         = addfile(s_FileRead, "r")
;
s_VarName2                      = "RegInModelDim_" + s_Exp
r_RegSigLev                     = FileReg->$s_VarName2$(:, :, {-90:90}, :)
r_RegSigLev                     = my_LonTest([/r_RegSigLev/], 360)
;
r_Reg                           = r_RegSigLev(0, :, :, :)
r_Reg                           = where(r_RegSigLev(1, :, :, :) .le. 0.1, r_Reg,\
                                        r_Reg@_FillValue)
;
; ====================================================================
;   3. read the constrained warming of MME and each model, and read 
; the present-day obs warming trend for further correction 
;; Obs Residual Warming Trend 
s_DirOut2 = "/global/cfs/cdirs/m1867/zmchen/my_data/Obs/TEMP/POST/YearMon/"
s_FileRead                      = s_DirOut2 + \
                "SurfaceAirTemp_Post_YearMon_1900-2014.nc"
FileObs                         = addfile(s_FileRead, "r")

s_VarNameObs                    = (/ "Cowtan_WayV2", "BEST", \
                                     "GISTEMP", "NOAAGlobalTemp" /)
;
opt                             = True
opt@DstGridLon                  = r_Reg&lon
opt@DstGridLat                  = r_Reg&lat
opt@DstReginal                  = True
opt@InterpMethod                = "conserve"
opt@ForceOverwrite              = True
;
do iObs = 0, dimsizes(s_VarNameObs)-1
    s_VarName2                  = "tas_" + s_VarNameObs(iObs)
    r_HadCRUT4                 := FileObs->$s_VarName2$(\
                                {i_yrConstrained(0):i_yrConstrained(1)}, :,:,:)
    r_HadCRUT4                  = where(abs(r_HadCRUT4) .lt. 100, r_HadCRUT4, \
                                        r_HadCRUT4@_FillValue)
    ;; Interpolation 
    opt@SrcGridLon             := r_HadCRUT4&lon
    opt@SrcGridLat             := r_HadCRUT4&lat
    r_HadCRUT4                 := ESMF_regrid(r_HadCRUT4, opt)
    system("rm -rf ./destination_grid_file.nc ./PET0.RegridWeightGen.Log ./source_grid_file.nc ./weights_file.nc")
    ;
    r_HadCRUT4                 := dim_avg_n_Wrap(r_HadCRUT4, 1)
    r_HadCRUT4!1                = "lat"
    r_HadCRUT4!2                = "lon"
    r_HadCRUT4_Glo             := my_AreaMean2([/r_HadCRUT4/], \
                                            r_HadCRUT4&lat, r_HadCRUT4&lon, 0)
    r_HadCRUT4_Glo             := conform(r_HadCRUT4, r_HadCRUT4_Glo, (/0/))
    copy_VarCoords(r_HadCRUT4, r_HadCRUT4_Glo)
    r_HadCRUT4                  = r_HadCRUT4 - r_HadCRUT4_Glo
    ;
    r_Temp_Trend                = trend_manken(r_HadCRUT4, False, 0)
    ;
    r_Temp_Trend                = where(abs(r_Temp_Trend) .lt. 100, r_Temp_Trend, \
                                        r_Temp_Trend@_FillValue)
    if (isdefined("r_TrendObs") .eq. False .or. iObs .eq. 0) then
        i_dims                 := array_append_record(dimsizes(s_VarNameObs), \
                                            dimsizes(r_Temp_Trend(0, :, :)), 0)
        r_TrendObs              = new(i_dims, float)
        copy_VarCoords(r_Temp_Trend(0, :, :), r_TrendObs(0, :, :))
        r_TrendObs!0            = "ObsData"
        r_TrendObs&ObsData      = ispan(1, dimsizes(s_VarNameObs), 1)
        r_TrendObs@ObsData      = s_VarNameObs
    end if
    r_TrendObs(iObs, :, :)      = r_Temp_Trend(1, :, :) * 10.
end do 
printMinMax(r_TrendObs, False)
r_TrendObsMean                  = dim_avg_n_Wrap(r_TrendObs, 0)
;
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/"
s_FileRead                      = s_DirRead + \
                "ConstrainedWarmingInChina_RegMMEForScaling" + s_4Models + ".nc"
FileConWarming                  = addfile(s_FileRead, "r")
s_VarName2                      = "ConstrainedWarmingChina_MME_" + s_Exp
r_ConWarmingMME                 = FileConWarming->$s_VarName2$
r_ConWarmingMME_Cor             = r_ConWarmingMME
;
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/GSAT_AndRegionalTrend_AllReal/"
s_FileRead                      = s_DirRead + "GSAT_AndRegionalTrend_" + \
            i_yrConstrained(0) + "-" + i_yrConstrained(1) + s_4Models + ".nc"
FileReTrend                     = addfile(s_FileRead, "r")
s_VarName2                      = "ResidualTrend_historical"
r_HistReTrend                   = FileReTrend->$s_VarName2$(1, :, :)
printMinMax(r_HistReTrend, False)
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
;;; Area-Mean over China 
r_China                := my_SelChina([/r_ConWarmingMME_Cor/])
r_China                := my_AreaMean2([/r_China/], r_China&lat, r_China&lon, 0)
r_China_Raw            := my_SelChina([/r_ConWarmingMME/])
r_China_Raw            := my_AreaMean2([/r_China_Raw/], r_China_Raw&lat, \
                                       r_China_Raw&lon, 0)
print("3 Term Warming:")
print(r_China + " K, " + r_China_Raw + " K")
print("")
;
s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/03.ConstraintByUsingLocalTrend/"
system("mkdir -p " + s_DirOut)
s_FileOut               = s_DirOut + "FurtherCorrectConstrainedWarming_" + \
                        s_Exp + s_4Models + ".nc"
ls                      = systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileOut             = addfile(s_FileOut, "c")
    GAttr               = True
    GAttr@history       = systemfunc("date") + "  " + systemfunc("pwd") + \
                        "/06.FurtherCorrectConstrainedWarming.ncl"
    fileattdef(FileOut, GAttr)
    ;
else 
    FileOut             = addfile(s_FileOut, "w")
end if
s_VarName2              = "ConstrainedWarming_" + s_Exp
r_ConWarmingMME_Cor@long_name = s_VarName2
FileOut->$s_VarName2$   = r_ConWarmingMME_Cor
;
;; for each CMIP6
do iModel = 0, dimsizes(s_ModelSelected)-1
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
    r_China_Raw            := my_SelChina([/r_ConWarming/])
    r_China_Raw            := my_AreaMean2([/r_China_Raw/], r_China_Raw&lat, \
                                       r_China_Raw&lon, 0)
    print("3 Term Warming:")
    print(r_China + " K, " + r_China_Raw + " K")
    print("")
end do
print(s_FileOut + "")

end





















