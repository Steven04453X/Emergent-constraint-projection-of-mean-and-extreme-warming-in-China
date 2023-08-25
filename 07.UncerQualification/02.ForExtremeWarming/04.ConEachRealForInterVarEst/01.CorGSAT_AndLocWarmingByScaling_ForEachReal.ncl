; ==================================================
; 本程序基于当前GSAT趋势与GSAT未来升温、GSAT与局地升温的回归系数，对每一个样本成员的升温进行约束
;   1. 指定模式、试验、年份范围
;   2. 读取回归系数，然后读取每一个样本成员的GSAT升温趋势，约束每一个样本的GSAT升温
;   3. 随后读取每一个样本的局地未来升温，基于GSAT和局地未来升温的回归系数，约束每一个样本的局地未来升温
; 注意本程序是针对每一个样本成员进行约束，以估计当前内部变率的贡献
; 注意本程序采用三个预估时期，GSAT与局地升温的回归系数，作为scaling的系数
; 注意本程序是通过订正的方式得到约束预估的升温，即先订正每一个模式的升温，然后
; 进行集合平均，这样可以避开截距系数的影响
; ==================================================
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
undef("my_3TermsMean")
function my_3TermsMean(Array, yrTerms)
local r_Array 
begin
;
i_yrProj                              = yrTerms
r_Array                               = Array
;
i_dims                                = dimsizes(r_Array)
if (dimsizes(i_dims) .eq. 1) then
    r_TempAllReal                     = r_Array(:2)
else if (dimsizes(i_dims) .eq. 2) then
    r_TempAllReal                     = r_Array(:2, :)
else if (dimsizes(i_dims) .eq. 3) then
    r_TempAllReal                     = r_Array(:2, :, :)
end if
end if
end if
r_TempAllReal!0                       = "Terms"
r_TempAllReal&Terms                  := ispan(0, 2, 1)
do iTerm = 0, dimsizes(i_yrProj(:, 0))-1
    i_TermYr                         := i_yrProj(iTerm, :)
    if (dimsizes(i_dims) .eq. 1) then
        r_TempAllReal(iTerm)          = dim_avg_n_Wrap(\
                                    r_Array({i_TermYr(0):i_TermYr(1)}), 0)
    else if (dimsizes(i_dims) .eq. 2) then
        r_TempAllReal(iTerm, :)       = dim_avg_n_Wrap(\
                                    r_Array({i_TermYr(0):i_TermYr(1)}, :), 0)
    else if (dimsizes(i_dims) .eq. 3) then
        r_TempAllReal(iTerm, :, :)    = dim_avg_n_Wrap(\
                                    r_Array({i_TermYr(0):i_TermYr(1)}, :, :), 0)
    end if 
    end if
    end if
end do 
;
return(r_TempAllReal)
;
end
;
begin
;
; ==================================================
;   1. 指定模式、试验、年份范围
s_4Models                 = "_18Models" ; "_5Models" ; "_30Models" ; "" ; 
s_4ModelsGSAT_Con         = "_23Models" ; Keep it as a largest ensemble
s_ModelSelected           = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
; s_ModelSel= (/ "ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR", \
;                 "CAMS-CSM1-0", "CanESM5", "CanESM5-CanOE", "CESM2-WACCM", \
;                 "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3", "EC-Earth3-Veg", \
;                 "FGOALS-f3-L", "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
;                 "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8", "INM-CM5-0", \
;                 "IPSL-CM6A-LR", "MCM-UA-1-0", "MIROC6", "MIROC-ES2L", \
;                 "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NESM3", \
;                 "NorESM2-LM", "UKESM1-0-LL"/)

s_Month                   = "AnMean"
if (isdefined("s_Exp") .eq. False) then
    s_Exp                 = "ssp245" ; "ssp126" ; "ssp585" ; "ssp370" ; 
end if
i_yrProj                  = (/ (/2021, 2040/), (/2041, 2060/), \
                                   (/2080, 2099/) /)
i_yrConstrained           = (/1981, 2014/)
;
;; For Output 
s_DirOut = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/02.ForExtremeWarming/04.ConEachRealForInterVarEst/"
system("mkdir -p " + s_DirOut)
s_FileOut                        = s_DirOut + \
                    "ConstrainedWarmingInChina_RegMMEForScaling_ForEachReal_"+ \
                    s_Exp + s_4Models + ".nc"
ls                               = systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileOut                      = addfile(s_FileOut, "c")
    GAttr                        = True
    GAttr@history                = systemfunc("date")+";  "+systemfunc("pwd") +\
                        "/01.CorGSAT_AndLocWarmingByScaling_ForEachReal.ncl"
    GAttr@Modles                 = s_ModelSelected
    fileattdef(FileOut, GAttr)
else 
    FileOut                      = addfile(s_FileOut, "w")
end if
;
s_FileOut                        = s_DirOut + "InterVarInProj_" + s_Exp + \
                                    s_4Models + ".nc"
ls                               = systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileOut_ProjInterVar         = addfile(s_FileOut, "c")
    GAttr                        = True
    GAttr@history                = systemfunc("date")+";  "+systemfunc("pwd") +\
                        "/01.CorGSAT_AndLocWarmingByScaling_ForEachReal.ncl"
    GAttr@Modles                 = s_ModelSelected
    fileattdef(FileOut_ProjInterVar, GAttr)
else 
    FileOut_ProjInterVar         = addfile(s_FileOut, "w")
end if
;
;; Obs Warming Trend 
s_DirOut2 = "/global/cfs/cdirs/m1867/zmchen/my_data/Obs/TEMP/POST/YearMon/"
s_FileRead                      = s_DirOut2 + \
                "SurfaceAirTemp_Post_YearMon_1900-2014.nc"
; s_FileRead                      = s_DirOut2 + \
;                 "SurfaceAirTemp_Post_YearMon_5x5_1900-2014.nc"
FileObs                         = addfile(s_FileRead, "r")

; s_VarNameObs                    = (/ "Cowtan_WayV2", "GISTEMP" /)
s_VarNameObs                    = (/ "Cowtan_WayV2", "BEST", \
                                     "GISTEMP", "NOAAGlobalTemp" /)
r_TrendObs                      = new(dimsizes(s_VarNameObs), float)
do iObs = 0, dimsizes(s_VarNameObs)-1
    s_VarName2                  = "tas_" + s_VarNameObs(iObs)
    r_HadCRUT4                 := FileObs->$s_VarName2$(\
                                {i_yrConstrained(0):i_yrConstrained(1)}, :,:,:)
    r_HadCRUT4                 := dim_avg_n_Wrap(r_HadCRUT4, 1)
    r_HadCRUT4!1                = "lat"
    r_HadCRUT4!2                = "lon"
    r_HadCRUT4                 := my_AreaMean2([/r_HadCRUT4/], \
                                            r_HadCRUT4&lat, r_HadCRUT4&lon, 0)
    ;
    r_Temp_Trend                = trend_manken(r_HadCRUT4, False, 0)
    r_TrendObs(iObs)            = r_Temp_Trend(1)
end do 
r_TrendObs                      = r_TrendObs * 10.
print(r_TrendObs + "")
r_TrendSTD                      = stddev(r_TrendObs)
print("Obs STD: " + r_TrendSTD)
; exit()
r_TrendObs                     := avg(r_TrendObs)
r_TrendObs@units                = "K per decade"
;
;; Blend Effect of surface temperature
r_TrendObs                      = r_TrendObs + 0.013
r_TrendSTD                      = r_TrendSTD + 0.005
;
; =========================================================
;   2. 读取回归系数，然后读取每一个样本成员的GSAT升温趋势，约束每一个样本的GSAT升温
;   3. 随后读取每一个样本的局地未来升温，基于GSAT和局地未来升温的回归系数，约束每一个样本的局地未来升温
;
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/GSAT_AndRegionalTrend_AllReal/"
s_FileRead                = s_DirRead + "ConstrainedGSAT_inThreeTerms" + \
                            s_4ModelsGSAT_Con + ".nc"
FileConstraint            = addfile(s_FileRead, "r")
s_VarName2                = "RegYintercept_" + s_Exp
r_RegGSAT_Trend           = FileConstraint->$s_VarName2$(0, :)
;
s_FileRead                = s_DirRead + "GSAT_AndRegionalTrend_1981-2014" + \
                            s_4ModelsGSAT_Con + ".nc"
FileGSAT_Trend            = addfile(s_FileRead, "r")
;
;; Raw projection of GSAT 
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/GMSAT/MultiMembers/"
s_FileRead                = s_DirRead + "GMSAT_Series_MultiMember_" +s_Exp+".nc"
FileGMSAT                 = addfile(s_FileRead, "r")
s_FileRead                = s_DirRead + "GMSAT_Series_MultiMember_historical.nc"
FileGMSAT_Hist            = addfile(s_FileRead, "r")
;
s_DirRead   = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tasmax/POST/TXxRegGSAT/InTimeDim/"
s_FileRead              = s_DirRead + "GSAT_reg_TXxWarmingInTimeDim_AnMean_" + \
                            s_Exp + "_20yr" + s_4Models + ".nc"
FileRead                = addfile(s_FileRead, "r")
s_VarName2              = "Reg_MME_AcrossModel_" + s_Exp
r_Reg                   = FileRead->$s_VarName2$
;
do iModel = 0, dimsizes(s_ModelSelected)-1
    l_list               := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp, "day")
    s_Real               := l_list[1]
    l_list_              := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp, "day")
    s_Real_              := l_list_[1]
    if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real) ) then
        s_Real           := s_Real_
    end if
    do iReal = 0, dimsizes(s_Real)-1
        s_VarName2        = "GSAT_Trend_" + s_ModelSelected(iModel) + \
                            "_historical_" + s_Real(iReal)
        if (my_VarNameExist(FileGSAT_Trend, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue
        end if
        r_PD_GSAT_Trend  := FileGSAT_Trend->$s_VarName2$(1)
        ;
        r_dx             := r_PD_GSAT_Trend - r_TrendObs
        r_dy             := r_RegGSAT_Trend * r_dx
        ;
        ;; GMSAT 
        s_VarName2        = "GMSAT_" + s_ModelSelected(iModel) + "_" + s_Exp + \
                            "_" + s_Real(iReal)
        if (my_VarNameExist(FileGMSAT, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue
        end if
        r_TempGMSAT      := FileGMSAT->$s_VarName2$
        r_TempGMSAT      := dim_avg_n_Wrap(r_TempGMSAT, 1)
        s_VarName2        = "GMSAT_" + s_ModelSelected(iModel) + \
                                    "_historical_" + s_Real(iReal)
        r_TempGMSAT_Hist := FileGMSAT_Hist->$s_VarName2$({1995:2014}, :)
        r_TempGMSAT_Hist := dim_avg_n_Wrap(r_TempGMSAT_Hist, 0)
        r_TempGMSAT_Hist := dim_avg_n_Wrap(r_TempGMSAT_Hist, 0)
        r_TempGMSAT       = r_TempGMSAT - r_TempGMSAT_Hist
        ;
        r_TempGMSAT      := my_3TermsMean(r_TempGMSAT, i_yrProj)
        ; print(r_TempGMSAT + "")
        ;
        ;; GSAT constraint 
        r_TempGMSAT_Con  := r_TempGMSAT
        r_TempGMSAT_Con   = r_TempGMSAT - r_dy
        ;
        ;; local warming 
        s_VarName         = "TXx" ; "pr" ; 
        s_DirRead = "$WACCEM/my_data/CMIP6/tasmax/POST/TXx/AllReal/"
        s_FileRead        = s_DirRead + "TXx_AnMean_" + s_Exp + "_" + \
                            s_ModelSelected(iModel) + "_CMIP6_AllReal.nc"
        FileRead_Raw     := addfile(s_FileRead, "r")
        s_VarName2       := s_VarName + "_" + s_ModelSelected(iModel) + "_" + \
                            s_Exp + "_" + s_Real(iReal)
        if (my_VarNameExist(FileRead_Raw, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue 
        end if
        r_Temp           := FileRead_Raw->$s_VarName2$(:, {-90:90}, {0:360})
        r_Temp           := my_3TermsMean(r_Temp, i_yrProj)
        ;
        s_FileRead        = s_DirRead + "TXx_AnMean_historical_" + \
                            s_ModelSelected(iModel) + "_CMIP6_AllReal.nc"
        FileRead_RawHist := addfile(s_FileRead, "r")
        s_VarName2       := s_VarName + "_" + s_ModelSelected(iModel) + \
                            "_historical_" + s_Real(iReal)
        if (my_VarNameExist(FileRead_RawHist, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue 
        end if
        r_TempHist       := FileRead_RawHist->$s_VarName2$({1995:2014}, \
                                                           {-90:90}, {0:360})
        r_TempHist       := dim_avg_n_Wrap(r_TempHist, 0)
        r_TempHist       := conform(r_Temp, r_TempHist, (/1, 2/))
        copy_VarCoords(r_Temp, r_TempHist)
        r_Temp            = r_Temp - r_TempHist
        ;
        ;; local warming constraint 
        r_dx             := r_TempGMSAT - r_TempGMSAT_Con
        r_dx             := conform(r_Temp, r_dx, (/0/))
        copy_VarCoords(r_Temp, r_dx)
        r_dy             := r_Reg * r_dx 
        r_VarCon         := r_Temp 
        r_VarCon          = r_Temp - r_dy
        ;
        s_VarName2        = "ConstrainedWarmingChina_" + \
                        s_ModelSelected(iModel) + "_" + s_Exp + "_" + \
                        s_Real(iReal)
        r_VarCon@long_name= s_VarName2
        FileOut->$s_VarName2$ = r_VarCon
        if (isdefined("r_RawAll") .eq. False .or. iReal .eq. 0) then
            i_dims       := array_append_record(dimsizes(s_Real), \
                                                dimsizes(r_Temp), 0)
            r_RawAll     := new(i_dims, float)
            copy_VarCoords(r_VarCon, r_RawAll(0, :, :, :))
            r_RawAll!0    = "Real"

        end if
        r_RawAll(iReal, :, :, :) = r_Temp + 0
        print(s_Real(iReal) + " Done!")
    end do
    if (isdefined("r_RawAll")) then
        if (dimsizes(s_Real) .ge. 3) then
            r_RawAll     := dim_stddev_n_Wrap(r_RawAll, 0)
            s_VarName2   := "InterVarInProj_" + s_ModelSelected(iModel) + "_" +\
                            s_Exp
            r_RawAll@long_name    = s_VarName2
            FileOut_ProjInterVar->$s_VarName2$ = r_RawAll
        end if
    ;
        delete(r_RawAll)
    end if
    print(s_ModelSelected(iModel) + "Done!")
    print("")
end do
print(s_FileOut + "")
exit()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
end



































