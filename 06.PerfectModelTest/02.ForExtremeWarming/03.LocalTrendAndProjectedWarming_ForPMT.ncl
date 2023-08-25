; ===================================================================
; Author: Ziming Chen 
; This code is for: 
;   - cal the regression and correlation between the residual warming and present-day local trend 
;   - residual warming is cal by two methods: 
;       1) constrained local warming minus the constrained GSAT 
;       2) raw projection minus the constrained local warming
;   - before cal the present-day local trend, deduct the GSAT part 
;   1. indicate the models, exp, 
;   2. read the constrained local warming (it seems that we should cal it), 
; and then read the raw projection warming, and read the constrained GSAT 
; for each model finally
;   3. cal two residual wamring in the projection 
;   4. read the present-day (PD) residual trend
;   5. cal the inter-model correlation and regression btw the PD residual trend 
; and the two residual warming 
;   6. cal the inter-model regression btw the the PD residual trend 
; and the constrained local warming
;   6 output the results
; Attention: This code is for perfect mdoel test. So we should indicate the pseudo-obs model 
; ===================================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code//my_FunctionInPNNL.ncl"
;
undef("my_Reg_SigLev")
function my_Reg_SigLev(GSAT[*]:numeric, SSP[*]:numeric)
local r_GSAT2, r_SSP2
begin
;
r_GSAT2                       = GSAT 
r_SSP2                        = SSP 
;
        r_Reg                := regCoef_n(r_GSAT2, r_SSP2, 0, 0)
        ; copy_VarCoords(r_SSP2(0, :, :), r_Reg)
    ;
        r_tval               := r_Reg@tval 
        ; r_tval               := reshape(r_tval, dimsizes(r_Reg))
        ; copy_VarCoords(r_SSP2(0, :, :), r_tval)
        r_df                 := r_Reg@nptxy
        ; r_df                 := reshape(r_df, dimsizes(r_Reg))
        ; copy_VarCoords(r_tval, r_df)
    ;
        r_df                  = r_df - 2
        r_Temp               := r_tval
        r_Temp                = 0.5
        r_Sig                := betainc(r_df/(r_df + r_tval^2), r_df/2., r_Temp)
        ; copy_VarCoords(r_Temp, r_Sig)
    ;
        r_Temp               := r_Reg@yintercept
        ; r_Temp               := reshape(r_Temp, dimsizes(r_Reg))
        ; copy_VarCoords(r_Sig, r_Temp)
        ;
; i_dims                        = array_append_record(2, dimsizes(r_Reg), 0)
r_Reg_All                     = new(2, float)
; copy_VarCoords(r_Reg, r_Reg_All(0, :, :))
r_Reg_All!0                   = "Reg_SigLev"
r_Reg_All&Reg_SigLev          = (/0, 1/)
r_Reg_All(0)                  = r_Reg + 0
r_Reg_All(1)                  = r_Sig + 0
;
r_Reg_All@long_name           = "regression btw GSAT and local warming"
r_Temp@long_name              = "y intercept of regression"
return([/r_Reg_All, r_Temp/])
;
end
;
begin
;
; ===================================================================
;   1. indicate the models, exp, 
if (isdefined("s_PseudoModel") .eq. False) then
    s_PseudoModel               = "CanESM5" ; "EC-Earth3"
end if
;
s_4Models                       = "_18Models" ; "_5Models" ; "_30Models" ; 
s_4ModelsGSAT_Con               = "_23Models" ; Keep it as a largest ensemble
; s_ModelSelected                 = (/"ACCESS-CM2", "CNRM-CM6-1", "EC-Earth3", \
;                             "GFDL-CM4", "MPI-ESM1-2-LR"/)
s_ModelSelected                 = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
s_Month                         = "AnMean"
if (isdefined("s_Exp") .eq. False) then
    s_Exp                       = "ssp245" ; "ssp126" ;
end if
;
if (s_Exp .eq. "ssp126") then
    s_ModelSelected            := (/"ACCESS-CM2", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
end if
;
; ===================================================================
;   2. read the constrained local warming (it seems that we should cal it), 
; and then read the raw projection warming, and read the constrained GSAT 
; for each model finally
;; constrained local warming
s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/06.PerfectModelTest/02.ForExtremeWarming/"
s_FileRead                      = s_DirRead + \
                    "ConstrainedTXxWarmingInChina_RegMMEForScaling_PseModel_" + \
                    s_PseudoModel + s_4Models + ".nc"
FileConLocalWarming             = addfile(s_FileRead, "r")
;
; ;; raw projection warming 
; s_DirRead = "$WACCEM/my_data/CMIP6/tas/POST/ProjectionOfAllReal/"
; s_FileRead                      = s_DirRead + "ProjectionInDifTerms_AllReal" + \
;                                 s_4Models + ".nc"
; FileRawProj                     = addfile(s_FileRead, "r")
;
;; constrained GSAT 
s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/06.PerfectModelTest/02.ForExtremeWarming//"
s_FileRead                      = s_DirRead + \
                                "ConstrainedGSAT_inThreeTerms_PseModel_" + \
                                s_PseudoModel + s_4Models + ".nc"
FileConGSAT                     = addfile(s_FileRead, "r")
;
;; residual trend 
;;;;;;; Attention: Here is only 22 Models now!
s_DirRead = "$WACCEM/my_data/CMIP6/tas/POST/GSAT_AndRegionalTrend_AllReal/"
s_FileRead                      = s_DirRead + "GSAT_AndRegionalTrend_1981-2014" + \
                                s_4ModelsGSAT_Con + ".nc"
FileReTrend                     = addfile(s_FileRead, "r")
;;;;;;; Attention: Here is only 22 Models now!
; ===================================================================
;   3. cal two residual wamring in the projection 
do iModel = 0, dimsizes(s_ModelSelected)-1
    if (s_ModelSelected(iModel) .eq. s_PseudoModel) then
        print(s_ModelSelected(iModel) + " is pseudo-obs model")
        print("Skip!")
        continue
    end if
    ;       1) constrained local warming minus the constrained GSAT 
    s_VarName2                  = "ConstrainedWarmingChina_" + \
                                s_ModelSelected(iModel) + s_Exp
    if (my_VarNameExist(FileConLocalWarming, s_VarName2) .eq. False) then
        print(s_VarName2 + "  Missing!")
        continue
    end if
    r_ConLocWarming            := FileConLocalWarming->$s_VarName2$
    ;
    ; ==========================================================
    ;   4. read the present-day (PD) residual trend
    l_list                     := my_CMIP6Model_Exp_Path_ForAllReal(\
                                s_ModelSelected(iModel), s_Exp, "day")
    ; s_ExpModel_Path            := l_list[0]
    s_Real                     := l_list[1]
    l_list_                    := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                s_ModelSelected(iModel), s_Exp, "day")
    ; s_ExpModel_Path_           := l_list_[0]
    s_Real_                    := l_list_[1]
    if (all(ismissing(s_Real)) .or. dimsizes(s_Real_) .gt. dimsizes(s_Real) ) then
        s_Real             := s_Real_
            ; s_ExpModel_Path    := s_ExpModel_Path_
    end if
    do iReal = 0, dimsizes(s_Real)-1
        s_VarName2              = "ResidualTrend_" + s_ModelSelected(iModel) + \
                                "_historical_" + s_Real(iReal)
        if (my_VarNameExist(FileReTrend, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue
        end if
        r_Temp                 := FileReTrend->$s_VarName2$(1, :, :)
        ;
        if (isdefined("r_ReTrendAllReal") .eq. False .or. iReal .eq. 0) then
            i_dims             := array_append_record(dimsizes(s_Real), \
                                                      dimsizes(r_Temp), 0)
            r_ReTrendAllReal   := new(i_dims, float)
            copy_VarCoords(r_Temp, r_ReTrendAllReal(0, :, :))
            r_ReTrendAllReal!0  = "Real"
            r_ReTrendAllReal&Real = s_Real
        end if
        r_ReTrendAllReal(iReal, :, :) = r_Temp
    end do
    if (isdefined("r_ReTrendAllReal") .eq. False) then
        print("r_ReTrendAllReal Missing! " + s_ModelSelected(iModel))
        continue
    end if
    r_ReTrendAllReal           := dim_avg_n_Wrap(r_ReTrendAllReal, 0)
    ;
    if (isdefined("r_ReTrendAllModel") .eq. False .or. iModel .eq. 0) then
        i_dims                 := array_append_record(dimsizes(s_ModelSelected), \
                                                      dimsizes(r_ReTrendAllReal), 0)
        r_ReTrendAllModel      := new(i_dims, float)
        copy_VarCoords(r_ReTrendAllReal, r_ReTrendAllModel(0, :, :))
        r_ReTrendAllModel!0     = "Models"
        r_ReTrendAllModel&Models= s_ModelSelected
        ;
        i_dims                 := array_append_record(dimsizes(s_ModelSelected), \
                                                      dimsizes(r_ConLocWarming), 0)
        r_ConLocWarmingAllModel = new(i_dims, float)
        copy_VarCoords(r_ConLocWarming, r_ConLocWarmingAllModel(0, :, :, :))
        copy_VarCoords(r_ReTrendAllModel(:, 0, 0), \
                       r_ConLocWarmingAllModel(:, 0, 0, 0))
    end if
    r_ReTrendAllModel(iModel, :, :) = r_ReTrendAllReal
    r_ConLocWarmingAllModel(iModel, :, :, :) = r_ConLocWarming
    ;
    delete([/r_ReTrendAllReal/])
end do
; ====================================================================
;   5. cal the inter-model correlation and regression btw the PD residual trend 
; and the two residual warming 
;; Regression for each model 
i_dims                         := array_append_record(2, \
                            dimsizes(r_ConLocWarmingAllModel(0, :, :, :)), 0)
r_Reg_All                       = new(i_dims, float)
        copy_VarCoords(r_ConLocWarmingAllModel(0,:,:,:), r_Reg_All(0, :, :, :))
r_Reg_All!0                     = "RegAndSigLevel"
r_Reg_All&RegAndSigLevel       := (/0, 1/)
r_yintercept                   := r_Reg_All(0, :, :, :)
        ;
i_dimsPeriod                   := 3
do iPeriod = 0, i_dimsPeriod-1
    do ilat = 0, dimsizes(r_Reg_All&lat)-1
        do ilon = 0, dimsizes(r_Reg_All&lon)-1
            l_list             := my_Reg_SigLev(\
                                r_ReTrendAllModel(:, ilat, ilon),\
                                r_ConLocWarmingAllModel(:, iPeriod, ilat, ilon))
            r_Reg_All(:, iPeriod, ilat, ilon) = l_list[0]
            r_yintercept(iPeriod, ilat, ilon) = l_list[1]
        end do 
    end do 
end do


;

;   6. cal the inter-model regression btw the the PD residual trend 
; and the constrained local warming
;
; ===============================================================
;   6 output the results
s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/06.PerfectModelTest/02.ForExtremeWarming/"
system("mkdir -p " + s_DirOut)
s_FileOut   = s_DirOut + "Regression_LocalTrendAndProjectedWarming_PseModel_" +\
            s_PseudoModel + "_" + s_Exp + s_4Models + ".nc"
system("rm -rf " + s_FileOut)
FileOut     = addfile(s_FileOut, "c")
GAttr                       = True
GAttr@history               = systemfunc("date") + "  " + systemfunc("pwd") +\
                                "/03.LocalTrendAndProjectedWarming_ForPMT.ncl"
GAttr@Des                   = "constrained loacl warming regresses onto local trend without GSAT warming across warming"
GAttr@Exp                   = s_Exp
GAttr@Models                = s_ModelSelected
    fileattdef(FileOut, GAttr)
;; Output 
        s_VarName2               := "RegInTimDim_" + s_Exp
        r_Reg_All@long_name       = s_VarName2
        FileOut->$s_VarName2$     = r_Reg_All
        ;
        s_VarName2                = "yintercept_" + s_Exp
        r_yintercept@long_name    = s_VarName2
        FileOut->$s_VarName2$     = r_yintercept
;
print(s_FileOut + "")
end











