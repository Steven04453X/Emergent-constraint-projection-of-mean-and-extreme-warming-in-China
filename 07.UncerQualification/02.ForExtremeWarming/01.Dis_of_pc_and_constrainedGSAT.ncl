; ============================================================
; Author: Ziming Chen 
; This code is for:
;   - cal the variance explained by the emergent constraint 
;   - this is the first step 
;   - providing the PDF for the regression coef which is used 
;   - to constrain the GSAT 
;   1. indicate the models, exp, projection terms, constrained 
; time range 
;   2. read the present-day GSAT trend, projected GSAT 
;   3. produce a random number array with the shape of 1000 x N 
;   4. for each random group, cal the regression coef, and then 
; constrained the GSAT projection 
;   5. output the results 
; ============================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code/my_FunctionInPNNL.ncl"
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
if (dimsizes(i_dims) .eq. 3) then
    r_TempAllReal                     = r_Array(:2, :, :)
else if (dimsizes(i_dims) .eq. 2) then
    r_TempAllReal                     = r_Array(:2, :)
else if (dimsizes(i_dims) .eq. 1) then
    r_TempAllReal                     = r_Array(:2)
end if
end if
end if
r_TempAllReal!0                       = "Terms"
r_TempAllReal&Terms                  := ispan(0, 2, 1)
do iTerm = 0, dimsizes(i_yrProj(:, 0))-1
    i_TermYr                         := i_yrProj(iTerm, :)
    if (dimsizes(i_dims) .eq. 1) then
        r_TempAllReal(iTerm)          = avg(r_Array({i_TermYr(0):i_TermYr(1)}))
    else if (dimsizes(i_dims) .eq. 2) then
        r_TempAllReal(iTerm, :, :)    = dim_avg_n_Wrap(\
                                    r_Array({i_TermYr(0):i_TermYr(1)}, :), 0)
    else if (dimsizes(i_dims) .eq. 3) then
        r_TempAllReal(iTerm, :, :)    = dim_avg_n_Wrap(\
                                    r_Array({i_TermYr(0):i_TermYr(1)}, :, :), 0)
    end if
    end if
    end if
    ;
end do 
;
return(r_TempAllReal)
;
end
;
begin
;
; ============================================================
;   1. indicate the models, exp, projection terms, constrained 
; time range 
s_ModelSelect           = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
s_4Models               = "_18Models"
s_4ModelsGSAT_Con       = "_23Models" ; Keep it as a largest ensemble
if (isdefined("s_Exp") .eq. False) then
    s_Exp               = "ssp245"
end if
if (s_Exp .eq. "ssp126") then
    s_ModelSelected    := (/"ACCESS-CM2", "AWI-CM-1-1-MR", \
                        "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                        "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
end if
i_yrBaseState           = (/1995, 2014/)
i_yrConstrained         = (/1981, 2014/)
i_yrProj                = (/ (/2021, 2040/), (/2041, 2060/), (/2080, 2099/) /)
s_yrProj                = (/"near", "mid", "long"/)
;
s_VarName               = "tas"
;
; ============================================================
;   2. read the present-day GSAT trend, projected GSAT 
;; the present-day GSAT trend in the pseduo-obs model 
s_DirRead = "$WACCEM/my_data/CMIP6/tas/POST/GSAT_AndRegionalTrend_AllReal/"
s_FileRead                = s_DirRead + "GSAT_AndRegionalTrend_1981-2014" + \
                            s_4ModelsGSAT_Con + ".nc"
FilePD_Trend              = addfile(s_FileRead, "r")
;
s_DirRead  = "$WACCEM/my_data/CMIP6/tas/POST/GMSAT/MultiMembers/"
s_FileRead                = s_DirRead + "GMSAT_Series_MultiMember_"+s_Exp + ".nc"
FileProjGSAT              = addfile(s_FileRead, "r")
s_FileRead                = s_DirRead + "GMSAT_Series_MultiMember_historical.nc"
FileHistGSAT              = addfile(s_FileRead, "r")
;
s_DirRead  = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/ProjectionOfAllReal/"
s_FileRead                = s_DirRead + "ProjectionInDifTerms_AllReal" + \
                            s_4ModelsGSAT_Con + ".nc"
FileTasProjChanges        = addfile(s_FileRead, "r")
;
do iModel = 0, dimsizes(s_ModelSelect)-1
    ;
    l_list               := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp, "Amon")
    s_Real               := l_list[1]
    l_list_              := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp, "Amon")
    s_ExpModel_Path_     := l_list_[0]
    s_Real_              := l_list_[1]
    if (all(ismissing(s_Real)) .or. dimsizes(s_Real_).gt.dimsizes(s_Real)) then
        s_Real           := s_Real_
    end if
    ;
    do iReal = 0, dimsizes(s_Real)-1    
        ;; PD Trend 
        s_VarNamePD_GSAT_Trend  = "GSAT_Trend_" + s_ModelSelect(iModel)\
                                    + "_historical_" + s_Real(iReal)
        if (my_VarNameExist(FilePD_Trend, s_VarNamePD_GSAT_Trend) .eq. False) then
            print(s_VarNamePD_GSAT_Trend + "  Missing!")
            continue
        end if;
        r_PD_GSAT_Trend        := FilePD_Trend->$s_VarNamePD_GSAT_Trend$(1)
    ;
        ;; GSAT 
        s_VarName2              = "GMSAT_" + s_ModelSelect(iModel) + "_" + \
                                    s_Exp + "_" + s_Real(iReal)
        if (my_VarNameExist(FileProjGSAT, s_VarName2) .eq. False) then
            print(s_VarName2 + "  Missing!")
            continue
        end if;
        r_GSAT                 := FileProjGSAT->$s_VarName2$
        r_GSAT                 := dim_avg_n_Wrap(r_GSAT, 1)
        ;
        s_VarName2              = "GMSAT_" + s_ModelSelect(iModel) + \
                                "_historical_" + s_Real(iReal)
        r_GSAT_Hist            := FileHistGSAT->$s_VarName2$(\
                                        {i_yrBaseState(0):i_yrBaseState(1)}, :)
        r_GSAT_Hist            := dim_avg_n_Wrap(r_GSAT_Hist, 1)
        r_GSAT_Hist            := dim_avg_n_Wrap(r_GSAT_Hist, 0)
        ;
        r_GSAT                  = r_GSAT - r_GSAT_Hist
        ;
        ; ;; tas projection 
        ; s_VarNameProj          := "tas_inDifTerms_" + s_ModelSelect(iModel) + \
        ;                             "_" + s_Exp(iExp) + "_" + s_Real(iReal)
        ; if (my_VarNameExist(FileProjChanges, s_VarNameProj) .eq. False) then
        ;     print(s_VarNameProj + "  Missing!")
        ;     continue
        ; end if
        ; r_ProjTerm             := FileTasProjChanges->$s_VarNameProj$
        ;
        if (isdefined("r_PD_GSAT_Trend_AllReal") .eq. False .or. iReal .eq. 0) then
            i_dims             := dimsizes(s_Real)
            r_PD_GSAT_Trend_AllReal    := new(i_dims, float)
            r_PD_GSAT_Trend_AllReal!0   = "Real"
            r_PD_GSAT_Trend_AllReal&Real= s_Real
            ;
            i_dims             := array_append_record(dimsizes(s_Real), \
                                                      dimsizes(r_GSAT), 0)
            r_GSAT_AllReal     := new(i_dims, float)
            copy_VarCoords(r_PD_GSAT_Trend_AllReal, r_GSAT_AllReal(:, 0))
            copy_VarCoords(r_GSAT, r_GSAT_AllReal(0, :))
            ;
            ; i_dims                     := array_append_record(\
            ;                         dimsizes(s_Real), dimsizes(r_ProjTerm), 0)
            ; r_ProjTerm_AllReal         := new(i_dims, float)
            ; copy_VarCoords(r_ProjTerm, r_ProjTerm_AllReal(0, :, :, :))
            ; copy_VarCoords(r_PD_GSAT_Trend_AllReal, r_ProjTerm_AllReal(:,0,0,0))
        end if
        r_PD_GSAT_Trend_AllReal(iReal)  = r_PD_GSAT_Trend
        r_GSAT_AllReal(iReal, :)        = r_GSAT
        ; r_ProjTerm_AllReal(iReal,:,:,:) = r_ProjTerm
    end do 
    if (isdefined("r_PD_GSAT_Trend_AllReal") .eq. False) then
        print("r_PD_GSAT_Trend_All Missing!")
        continue
    end if
    if (isdefined("r_GSAT_AllReal") .eq. False) then
        print("r_GSAT_AllReal Missing!")
        continue
    end if
    r_PD_GSAT_Trend_AllReal    := dim_avg_n_Wrap(r_PD_GSAT_Trend_AllReal, 0)
    r_GSAT_AllReal             := dim_avg_n_Wrap(r_GSAT_AllReal, 0)
    ; r_ProjTerm_GSAT_AllReal    := dim_avg_n_Wrap(r_ProjTerm_GSAT_AllReal, 0)
    ;
    if (isdefined("r_PD_GSAT_Trend_All") .eq. False .or. iModel .eq. 0) then
        i_dims                 := dimsizes(s_ModelSelect)
        r_PD_GSAT_Trend_All     = new(i_dims, float)
        r_PD_GSAT_Trend_All!0   = "Models"
        r_PD_GSAT_Trend_All&Models = s_ModelSelect
        ;
        i_dims                 := array_append_record(dimsizes(s_ModelSelect), \
                                                      3, 0)
        r_GSAT_All             := new(i_dims, float)
        copy_VarCoords(r_PD_GSAT_Trend_All, r_GSAT_All(:, 0))
        r_GSAT_All!1            = "Terms"
        r_GSAT_All&Terms        = ispan(0, 2, 1)
        ;
    end if
    ;
    r_PD_GSAT_Trend_All(iModel) = r_PD_GSAT_Trend_AllReal
    r_GSAT_All(iModel, :)       = my_3TermsMean(r_GSAT_AllReal, i_yrProj)
    delete([/r_GSAT_AllReal, r_PD_GSAT_Trend_AllReal/])
end do
;
; ============================================================
;   3. produce a random number array with the shape of 1000 x N
nBoot                           = 100
nDim                            = 0
opt                             = True
opt@rseed3                      = "clock"
;
i_dims                         := array_append_record(dimsizes(s_yrProj), \
                                                      nBoot, 0)
r_rcBoot                       := new(i_dims, float)
copy_VarCoords(r_GSAT_All(0, :), r_rcBoot(:, 0))
r_rcBoot!1                      = "boot"
r_rcBoot&boot                   = ispan(0, nBoot-1, 1)
r_rcBootAvg                     = r_rcBoot(:, 0)
r_rcBootStd                     = r_rcBoot(:, 0)
do iTerm = 0, dimsizes(s_yrProj)-1
    BootStrap                  := bootstrap_regcoef(r_PD_GSAT_Trend_All, \
                                                    r_GSAT_All(:, iTerm), nBoot, \
                                                    nDim, opt)
    rcBoot                     := BootStrap[0] ; ; Bootstrapped regression coefficients in ascending order
    ; printVarSummary(rcBoot)
    ; printMinMax(rcBoot, False)
;
    rcBootAvg                  := BootStrap[1] ; Average of the boot strapped regression coefficients
    rcBootStd                  := BootStrap[2] ; ; Standard deviation of bootstrapped regression coefficients
    ;
    r_rcBoot(iTerm, :)          = rcBoot
    r_rcBootAvg(iTerm)          = rcBootAvg
    r_rcBootStd(iTerm)          = rcBootStd
end do
;
; ============================================================ 
;   4. for each random group, cal the regression coef, and then 
; constrained the GSAT projection 
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
r_ConstraintModels             := new((/nBoot, dimsizes(s_yrProj), \
                                       dimsizes(s_ModelSelect)/), float)
copy_VarCoords(r_rcBoot(0, :), r_ConstraintModels(:, 0, 0))
copy_VarCoords(r_rcBoot(:, 0), r_ConstraintModels(0, :, 0))
r_ConstraintModels!2            = "Models"
r_ConstraintModels&Models       = s_ModelSelect
do iTerm = 0, dimsizes(s_yrProj)-1
    do iModel = 0, dimsizes(s_ModelSelect)-1
        r_delta_x              := r_PD_GSAT_Trend_All(iModel) - r_TrendObs
        r_delta_y              := r_rcBoot(iTerm, :) * r_delta_x
            ;
        r_ConstraintModels(:, iTerm, iModel) = r_GSAT_All(iModel, iTerm) - \
                                                r_delta_y
    end do
end do
printVarSummary(r_ConstraintModels)
;
; ============================================================
;   5. output the results 
s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/02.ForExtremeWarming/"
system("mkdir -p " + s_DirOut)
s_FileOut                       = s_DirOut + "Dis_of_pc_and_constrainedGSAT_" +\
                                s_Exp + s_4Models + ".nc"
system("rm -rf " + s_FileOut)
setfileoption("nc", "Format", "LargeFile")
FileOut                         = addfile(s_FileOut, "c")
GAttr                           = True
GAttr@history                   = systemfunc("date") + "  " + systemfunc("pwd")\
                                + "01.Dis_of_pc_and_constrainedGSAT.ncl"
fileattdef(FileOut, GAttr)
s_VarName2                      = "rcBoot_pc_" + s_Exp
FileOut->$s_VarName2$           = r_rcBoot
;
do iModel = 0, dimsizes(s_ModelSelect)-1
    s_VarName2                  = "rcBoot_ConstrainedGSAT_" + s_Exp + "_" + \
                                    s_ModelSelect(iModel)
    r_Temp                     := r_ConstraintModels(:, :, iModel)
    FileOut->$s_VarName2$       = r_Temp
end do
;
print(s_FileOut + "")
;
end











