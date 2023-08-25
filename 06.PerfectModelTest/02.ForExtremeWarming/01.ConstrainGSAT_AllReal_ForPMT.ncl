; =====================================================================
; Author: Ziming Chen 
; This code is for:
;   - constraining the global surface air temp (GSAT) in all the realization of each model 
;   - Output the constrained GSAT warming in different scenarios 
;   1. indicate the model, exp, periods and their name 
;   2. for model, get the name of realization based on the specific scenarios, 
; and then check whether this realization exists in historical and future projection 
;   3. if so, read the data, and then cal the present-day warming trend, and the 
; projected warming in the different projection terms 
;   4. store up and output the trend and warming of each realization, and then 
; cal the inter-member ensemble (IME)
;   5. cal the regression and correlation btw present-day trend and projected warming
; by using the IME 
;   6. read the obs trend, and then cal the constrained warming
; Attention: This code is for perfect mdoel test. So we should indicate the pseudo-obs model 
; =====================================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code/my_FunctionInPNNL.ncl"
;
begin
;
; =====================================================================
;   1. indicate the model, exp, periods and their name 
;
if (isdefined("s_PseudoModel") .eq. False) then
    s_PseudoModel       = "CanESM5" ; "MRI-ESM2-0"
end if
;
s_ModelSelect           = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)

s_4Models               = "_18Models"
s_4ModelsGSAT_Con       = "_23Models" ; Keep it as a largest ensemble
;
if (isdefined("s_Exp") .eq. False) then
    s_Exp               = "ssp245" ; "ssp126" ; "ssp370" ; "ssp585" ; 
end if
;
if (s_Exp .eq. "ssp126") then
    s_ModelSelected    := (/"ACCESS-CM2", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
end if
;
i_yrBaseState           = (/1995, 2014/)
i_yrConstrained         = (/1981, 2014/)
i_yrProj                = (/ (/2021, 2040/), (/2041, 2060/), (/2080, 2099/) /)
s_yrProj                = (/"near", "mid", "long"/)
;
s_VarName               = "tas"
;
;;;;;;; Attention: Here is only 22 Models now!
;; the present-day GSAT trend in the pseduo-obs model 
s_DirRead = "$WACCEM/my_data/CMIP6/tas/POST/GSAT_AndRegionalTrend_AllReal/"
s_FileRead                = s_DirRead + "GSAT_AndRegionalTrend_1981-2014" + \
                            s_4ModelsGSAT_Con + ".nc"
FilePD_Trend              = addfile(s_FileRead, "r")
;;;;;;; Attention: Here is only 22 Models now!
;
; =====================================================================
;   2. for model, get the name of realization based on the specific scenarios, 
; and then check whether this realization exists in historical and future projection 
;
;; for output file 
s_DirOut  = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/06.PerfectModelTest/02.ForExtremeWarming/"
system("mkdir -p " + s_DirOut)
s_FileOut               = s_DirOut + "ConstrainedGSAT_inThreeTerms_PseModel_" +\
                         s_PseudoModel + s_4Models + ".nc"

ls                      = systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileConstraint      = addfile(s_FileOut, "c")
    GAttr               = True
    GAttr@history       = systemfunc("date") + "  " + systemfunc("pwd") + \
                        "/01.ConstrainGSAT_AllReal_ForPMT.ncl"
    fileattdef(FileConstraint, GAttr)
else 
    FileConstraint      = addfile(s_FileOut, "w")
end if
;
;;;;;;; Attention: Here is only 22 Models now!
s_DirRead  = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/ProjectionOfAllReal/"
s_FileRead              = s_DirRead + "ProjectionInDifTerms_AllReal" + \
                            s_4ModelsGSAT_Con + ".nc"
FileProjChanges         = addfile(s_FileRead, "r")
;;;;;;; Attention: Here is only 22 Models now!
;
;; Read data 
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/2.5x2.5/AllReal/"
;
do iExp = 0, dimsizes(s_Exp)-1
    do iModel = 0, dimsizes(s_ModelSelect)-1
        l_list         := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "Amon")
        s_ExpModel_Path:= l_list[0]
        s_Real         := l_list[1]
        l_list_        := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "Amon")
        s_ExpModel_Path_:= l_list_[0]
        s_Real_         := l_list_[1]
        if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real) ) then
            s_Real      := s_Real_
            s_ExpModel_Path := s_ExpModel_Path_
        end if
        ;
        print("Models No.: " + iModel + ", " + s_ModelSelect(iModel) + \
              "; Exp: " + s_Exp(iExp))
        ;
        ; =================================================================
        ;   3. if so, read the data, and then cal the present-day warming trend, and the 
        ; projected warming in the different projection terms 
        print(s_ModelSelect(iModel) + "")
        do iReal = 0, dimsizes(s_Real)-1
            ;
            s_VarNameProj              := "tas_inDifTerms_" + \
                                    s_ModelSelect(iModel) + "_" + s_Exp(iExp) +\
                                    "_" + s_Real(iReal)
            s_VarNamePD_GSAT_Trend      = "GSAT_Trend_" + s_ModelSelect(iModel)\
                                    + "_historical_" + s_Real(iReal)
            if (s_ModelSelect(iModel) .ne. s_PseudoModel) then
                ;; Being constrained model 
                if (my_VarNameExist(FileProjChanges, s_VarNameProj) .eq. \
                    False) then
                    print(s_VarNameProj + "  Missing!")
                    continue
                end if
                r_ProjTerm             := FileProjChanges->$s_VarNameProj$
                ;
                r_PD_GSAT_Trend        := FilePD_Trend->$s_VarNamePD_GSAT_Trend$(1)
                ;
                ;; store up 
                if (isdefined("r_PD_GSAT_Trend_AllReal") .eq. False .or. \
                    iReal .eq. 0) then
                    r_PD_GSAT_Trend_AllReal    := new(dimsizes(s_Real), float)
                    r_PD_GSAT_Trend_AllReal!0   = "Real"
                    r_PD_GSAT_Trend_AllReal&Real= s_Real
                    ;
                    i_dims                     := array_append_record(\
                                    dimsizes(s_Real), dimsizes(r_ProjTerm), 0)
                    r_ProjTerm_AllReal         := new(i_dims, float)
                    copy_VarCoords(r_ProjTerm, r_ProjTerm_AllReal(0, :, :, :))
                    copy_VarCoords(r_PD_GSAT_Trend_AllReal, \
                                   r_ProjTerm_AllReal(:, 0, 0, 0))
                end if
                r_PD_GSAT_Trend_AllReal(iReal)  = r_PD_GSAT_Trend + 0
                r_ProjTerm_AllReal(iReal,:,:,:) = r_ProjTerm
                print(s_Real(iReal) + " Being Constrained")
            else
                ;; Pseudo model 
                if (my_VarNameExist(FileProjChanges, s_VarNameProj) .eq. \
                    False) then
                    print(s_VarNameProj + "  Missing!")
                    continue
                end if
                r_ProjTerm_PseObs      := FileProjChanges->$s_VarNameProj$                
                ;
                r_PD_GSAT_Trend_PseObs := FilePD_Trend->$s_VarNamePD_GSAT_Trend$(1)
                ;
                if (isdefined("r_PD_GSAT_Trend_AllReal_PseObs") .eq. False .or. \
                    iReal .eq. 0) then
                    r_PD_GSAT_Trend_AllReal_PseObs    := new(dimsizes(s_Real), float)
                    r_PD_GSAT_Trend_AllReal_PseObs!0   = "Real"
                    r_PD_GSAT_Trend_AllReal_PseObs&Real= s_Real
                    ;
                    i_dims                     := array_append_record(\
                            dimsizes(s_Real), dimsizes(r_ProjTerm_PseObs), 0)
                    r_ProjTerm_AllReal_PseObs  := new(i_dims, float)
                    copy_VarCoords(r_ProjTerm_PseObs, \
                                   r_ProjTerm_AllReal_PseObs(0, :, :, :))
                    copy_VarCoords(r_PD_GSAT_Trend_AllReal_PseObs, \
                                   r_ProjTerm_AllReal_PseObs(:, 0, 0, 0))
                end if
                r_PD_GSAT_Trend_AllReal_PseObs(iReal)  = r_PD_GSAT_Trend_PseObs + 0
                r_ProjTerm_AllReal_PseObs(iReal, :,:,:)= r_ProjTerm_PseObs
                ;
                print(s_Real(iReal) + " Pseudo")
            end if
            ;
        end do
        ; ===========================================================
        ;   4. store up and output the trend and warming of each realization, and then 
        ; cal the inter-member ensemble (IME)
        ;
        if (s_ModelSelect(iModel) .ne. s_PseudoModel) then
            if (isdefined("r_PD_GSAT_Trend_AllReal") .eq. False) then
                print("r_PD_GSAT_Trend_AllReal Missing in " + \
                      s_ModelSelect(iModel) + " in " + s_Exp(iExp))
                continue
            end if
            r_PD_GSAT_Trend_AllReal    := dim_avg_n_Wrap(\
                                                    r_PD_GSAT_Trend_AllReal, 0)
            r_ProjTerm_GSAT_AllReal    := my_AreaMean2([/r_ProjTerm_AllReal/], \
                            r_ProjTerm_AllReal&lat, r_ProjTerm_AllReal&lon, 0)
            r_ProjTerm_GSAT_AllReal    := dim_avg_n_Wrap(r_ProjTerm_GSAT_AllReal, 0)
            ;; store up the IME of all models     
            if (isdefined("r_PD_GSAT_Trend_IME_AllModel") .eq. False .or. iModel .eq. 0) then
                r_PD_GSAT_Trend_IME_AllModel := new(dimsizes(s_ModelSelect), float)
                r_PD_GSAT_Trend_IME_AllModel!0= "Models"
                r_PD_GSAT_Trend_IME_AllModel&Models = s_ModelSelect
                ;
                i_dims                       := array_append_record(\
                                            dimsizes(s_ModelSelect), \
                                            dimsizes(r_ProjTerm_GSAT_AllReal), 0)
                r_ProjTerm_GSAT_IME_AllModel := new(i_dims, float)
                copy_VarCoords(r_PD_GSAT_Trend_IME_AllModel, \
                               r_ProjTerm_GSAT_IME_AllModel(:, 0))
                copy_VarCoords(r_ProjTerm_GSAT_AllReal, \
                               r_ProjTerm_GSAT_IME_AllModel(0, :))
            end if
            r_PD_GSAT_Trend_IME_AllModel(iModel)    = r_PD_GSAT_Trend_AllReal
            r_ProjTerm_GSAT_IME_AllModel(iModel, :) = r_ProjTerm_GSAT_AllReal
            delete([/r_PD_GSAT_Trend_AllReal/])
        else
            if (isdefined("r_PD_GSAT_Trend_AllReal_PseObs") .eq. False) then
                print("Pse-obs model missing! " + s_ModelSelect(iModel) + \
                      " in " + s_Exp(iExp))
                exit()
            end if
            r_PD_GSAT_Trend_AllReal_PseObs   := avg(r_PD_GSAT_Trend_AllReal_PseObs)
            r_ProjTerm_AllReal_PseObs        := my_AreaMean2(\
                                            [/r_ProjTerm_AllReal_PseObs/], \
                                            r_ProjTerm_AllReal_PseObs&lat, \
                                            r_ProjTerm_AllReal_PseObs&lon, 0)
            r_ProjTerm_AllReal_PseObs        := dim_avg_n_Wrap(\
                                                r_ProjTerm_AllReal_PseObs, 0)
            r_TrendObs                       := r_PD_GSAT_Trend_PseObs
            print("Pseudo Trend: " + r_TrendObs)
        end if
        ;
        print(s_ModelSelect(iModel) + ", " + s_Exp(iExp) + " Done!")
        print("")
    end do
    ;
    ; ===========================================================
    ;   5. cal the regression and correlation btw present-day trend and projected warming
    ; by using the IME
    ;
    r_Coef                  := new((/2, dimsizes(s_yrProj)/), float)
    copy_VarCoords(r_ProjTerm_GSAT_IME_AllModel(0, :), r_Coef(0, :))
    r_Coef!0                 = "CorCoef_SigLev"
    r_Coef&CorCoef_SigLev    = (/0, 1/)
    r_RegYintercept         := new((/2, dimsizes(s_yrProj)/), float)
    copy_VarCoords(r_ProjTerm_GSAT_IME_AllModel(0, :), r_RegYintercept(0, :))
    r_RegYintercept!0        = "Reg_yintercept"
    r_RegYintercept&Reg_yintercept = (/0, 1/)
    do iTerm = 0, dimsizes(s_yrProj)-1
        r_Return            := my_CorCoef_SigLevel(\
                                r_PD_GSAT_Trend_IME_AllModel,\
                                r_ProjTerm_GSAT_IME_AllModel(:, iTerm), False)
        print(s_yrProj(iTerm) + "")
        print(s_ModelSelect + ": " + \
              sprintf("%4.2f", r_PD_GSAT_Trend_IME_AllModel) + ", " + \
              sprintf("%4.2f", r_ProjTerm_GSAT_IME_AllModel(:, iTerm)) )
        print("")
        ; (/"CorCoef", "Regression", "SigLevel", "Degree", "yintercept"/)
        r_Coef(0, iTerm)     = r_Return(0) + 0
        r_Coef(1, iTerm)     = r_Return(2) + 0
        r_RegYintercept(0, iTerm) = r_Return(1) + 0
        r_RegYintercept(1, iTerm) = r_Return(4) + 0
    end do
    print(s_yrProj + ": " + sprintf("%4.2f", r_Coef(0, :)) + \
          " (p < " + sprintf("%4.2f", r_Coef(1, :)) + ")")
    ;
    ; printVarSummary(r_RegYintercept) ; [Reg_yintercept | 2] x [Terms | 3]
    ; printVarSummary(r_Coef) ; [CorCoef_SigLev | 2] x [Terms | 3]
    ;
    ; ===========================================================
    ;   6. read the obs trend, and then cal the constrained warming
    r_ConstraintMME        := new(dimsizes(s_yrProj), float)
    copy_VarCoords(r_Coef(0, :), r_ConstraintMME)
    ;
    r_ConstraintModels     := new((/dimsizes(s_yrProj), \
                                   dimsizes(s_ModelSelect)/), float)
    copy_VarCoords(r_Coef(0, :), r_ConstraintModels(:, 0))
    r_ConstraintModels!1    = "Models"
    r_ConstraintModels&Models = s_ModelSelect
    print("")
    do iTerm = 0, dimsizes(s_yrProj)-1
        ;; Constrained MME 
        r_ConstraintMME(iTerm) = r_RegYintercept(0, iTerm) * r_TrendObs + \
                                r_RegYintercept(1, iTerm)
        print("Constrained GSAT: " + s_yrProj(iTerm) + ", " + \
              r_ConstraintMME(iTerm) + "K, Raw: " + \
              avg(r_ProjTerm_GSAT_IME_AllModel(:, iTerm)) + "K")
        ;
        ;; for each model 
        do iModel = 0, dimsizes(s_ModelSelect)-1
            if (s_ModelSelect(iModel) .eq. s_PseudoModel) then
                print(s_ModelSelect(iModel) + " is pseudo-obs model")
                print("Warming: " + r_ProjTerm_AllReal_PseObs(iTerm) + " K")
                print("Skip!")
                print("")
                continue
            end if
            r_delta_x      := r_PD_GSAT_Trend_IME_AllModel(iModel) - r_TrendObs
            r_delta_y      := r_RegYintercept(0, iTerm) * r_delta_x
            ;
            r_ConstraintModels(iTerm, iModel) = \
                        r_ProjTerm_GSAT_IME_AllModel(iModel, iTerm) - r_delta_y
            ;
        end do
    end do
    ; printVarSummary(r_ConstraintModels) ; [Terms | 3] x [Models | 5]
    ;
    ;; output 
    s_VarName2              = "ConstrainedGSAT_MME_" + s_Exp(iExp)
    FileConstraint->$s_VarName2$ = r_ConstraintMME
    ;
    s_VarName2              = "RegYintercept_" + s_Exp(iExp)
    FileConstraint->$s_VarName2$ = r_RegYintercept
    s_VarName2              = "CoefSigLev_" + s_Exp(iExp)
    FileConstraint->$s_VarName2$ = r_Coef
    ;
    do iModel = 0, dimsizes(s_ModelSelect)-1
        if (s_ModelSelect(iModel) .eq. s_PseudoModel) then
            print(s_ModelSelect(iModel) + " is pseudo-obs model")
            print("Skip!")
            print("")
            continue
        end if
            ;
        s_VarName2          = "ConstrainedGSAT_" + s_ModelSelect(iModel) + "_"\
                            + s_Exp(iExp)
        r_Temp             := r_ConstraintModels(:, iModel)
        r_Temp@Values       = r_Temp
        ; r_Temp@Terms       := s_yrProj
        FileConstraint->$s_VarName2$ = r_Temp
        ;
    end do
    ;
    delete([/r_PD_GSAT_Trend_IME_AllModel, r_PD_GSAT_Trend_AllReal_PseObs/])
    print("")
    print("")
end do
print(s_FileOut + "")
;
end














