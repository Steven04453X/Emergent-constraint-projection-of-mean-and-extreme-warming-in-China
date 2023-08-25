; =============================================================
; Author: Ziming Chen 
; This code is for: 
;   - cal the pearson correlation coef and rank correlation coef 
; btw the constrained warming and the specific variables 
;   - output as a nc file 
;   1. indicate the models, exp, var name 
;   2. for each model, read all existed realization and then make the 
; inter-member mean 
;   3. cal the two correlation coef and then output 
; =============================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code/my_FunctionInPNNL.ncl"
;
begin
;
; =============================================================
;   1. indicate the models, exp, var name 
s_ModelSelect           = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                            "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                            "CNRM-ESM2-1", "CanESM5", "CanESM5-CanOE", \
                            "EC-Earth3", "EC-Earth3-Veg", "FGOALS-f3-L", \
                            "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
                            "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8", \
                            "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR", \
                            "MRI-ESM2-0", "UKESM1-0-LL"/)
; s_ModelSelect        = (/"GFDL-CM4", "CNRM-CM6-1", "EC-Earth3"/)
s_4Models            = "_23Models" ; "_4Models" ; 
s_Exp                = (/"ssp126", "ssp245", "ssp370", "ssp585"/)
;
s_VarName            = (/"pr", "rlus", "rsds", "huss", "psl", "hfls", "ts", \
                         "tasmax"/)
; =============================================================
;   2. for each model, read all existed realization and then make the 
; inter-member mean 
do iExp = 0, dimsizes(s_Exp)-1
    s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/01.ForMeanWarming/04.ConEachRealForInterVarEst/"
    s_FileRead       = s_DirRead+"FurtherCorrectConstrainedWarming_ForEachReal_"\
                       + s_Exp(iExp) + s_4Models + ".nc"
    FileConstrainedMean:= addfile(s_FileRead, "r")
    ;
    s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/02.ForExtremeWarming/04.ConEachRealForInterVarEst/"
    s_FileRead       = s_DirRead+"FurtherCorrectConstrainedWarming_ForEachReal_"\
                        + s_Exp(iExp) + "_18Models.nc"
    FileConstrainedExt := addfile(s_FileRead, "r")
;   ;
    do ivar = 0, dimsizes(s_VarName)-1
        s_DirRead   := "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/08.RelatedToOtherVar/"
        s_FileRead  := s_DirRead + s_VarName(ivar) + \
                    "_Clim_STD_Trend_historical_1981-2014" + s_4Models + ".nc"
        FileVar     := addfile(s_FileRead, "r")
        ;
        do iModel = 0, dimsizes(s_ModelSelect)-1
            l_list         := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "Amon")
            s_ExpModel_Path:= l_list[0]
            s_Real         := l_list[1]
            l_list_        := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "Amon")
            s_ExpModel_Path_:= l_list_[0]
            s_Real_         := l_list_[1]
            if (all(ismissing(s_Real)) .or. dimsizes(s_Real_).gt.dimsizes(s_Real)) then
                s_Real      := s_Real_
                s_ExpModel_Path := s_ExpModel_Path_
            end if
        ;
            print("Models No.: " + iModel + ", " + s_ModelSelect(iModel) + \
              "; Exp: " + s_Exp(iExp))
            iReal_      := 0
            do iReal = 0, dimsizes(s_Real)-1
                ;; var 
                s_VarName2  := s_VarName(ivar) + "_" + s_ModelSelect(iModel) + \
                                "_historical_" + s_Real(iReal)
                L_Return    := my_VarNameExist(FileVar, s_VarName2)
                if (L_Return .eq. False) then
                    print(s_VarName2 + " Missing!")
                    continue
                end if
                r_Var       := FileVar->$s_VarName2$
                ;
                s_VarName2  := "ConstrainedWarming_" + s_ModelSelect(iModel) + \
                               "_" + s_Exp(iExp) + "_AllReal"
                r_MeanWarming := FileConstrainedMean->$s_VarName2$(iReal_, :, \
                                                            {-90:90}, {0:360})
                ;
                r_MeanWarming := my_SelChina([/r_MeanWarming/])
                r_MeanWarming := my_AreaMean2([/r_MeanWarming/], \
                                                r_MeanWarming&lat, \
                                                r_MeanWarming&lon, 0)
                r_Var         := my_AreaMean2([/r_Var/], r_Var&lat, r_Var&lon,0)
                ;
                if (isdefined("r_VarAll") .eq. False .or. iReal .eq. 0) then
                    i_dims  := array_append_record(dimsizes(s_Real), \
                                                   dimsizes(r_Var), 0)
                    r_VarAll:= new(i_dims, float)
                    copy_VarCoords(r_Var, r_VarAll(0, :))
                    r_VarAll!0 = "Real"
                    r_VarAll&Real  = s_Real
                    ;
                    i_dims  := array_append_record(dimsizes(s_Real), \
                                                   dimsizes(r_MeanWarming), 0)
                    r_MeanWarmingAll := new(i_dims, float)
                    copy_VarCoords(r_MeanWarming, r_MeanWarmingAll(0, :))
                    copy_VarCoords(r_VarAll(:, 0), r_MeanWarmingAll(:, 0))
                    ;
                    r_ConExtAll      := r_VarAll(:, 2)
                end if
                r_VarAll(iReal, :) = r_Var
                r_MeanWarmingAll(iReal, :) = r_MeanWarming
                ;
                ;; Constrained Extreme warming 
                s_VarName2   = "ConstrainedWarming_" + s_ModelSelect(iModel) + \
                                "_" + s_Exp(iExp) + "_AllReal"
                L_Return    := my_VarNameExist(FileConstrainedExt, s_VarName2)
                if (L_Return .eq. False) then
                    print(s_VarName2 + " Missing!")
                    continue 
                end if
                r_ConExt    := FileConstrainedExt->$s_VarName2$(iReal_, 2, \
                                                            {-90:90}, {0:360})
                r_ConExt    := my_SelChina([/r_ConExt/])
                r_ConExt    := my_AreaMean2([/r_ConExt/], r_ConExt&lat, \
                                            r_ConExt&lon, 0)
                r_ConExtAll(iReal) = r_ConExt
                ;
                iReal_       = iReal_ + 1
            ;
            end do 
            if (isdefined("r_VarAll")) then
                r_VarAll   := dim_avg_n_Wrap(r_VarAll, 0)
                r_MeanWarmingAll := dim_avg_n_Wrap(r_MeanWarmingAll, 0)
                r_ConExtAll      := avg(r_ConExtAll)
                ;
                if (isdefined("r_VarAllModels").eq.False .or. iModel.eq.0) then
                    i_dims       := array_append_record(dimsizes(s_ModelSelect), \
                                                        dimsizes(r_VarAll), 0)
                    r_VarAllModels := new(i_dims, float)
                    copy_VarCoords(r_VarAll, r_VarAllModels(0, :))
                    r_VarAllModels!0= "Models"
                    r_VarAllModels&Models = s_ModelSelect
                    ;
                    i_dims       := array_append_record(dimsizes(s_ModelSelect), \
                                                dimsizes(r_MeanWarmingAll), 0)
                    r_MeanWarmingAllModels := new(i_dims, float)
                    copy_VarCoords(r_MeanWarmingAll, r_MeanWarmingAllModels(0, :))
                    copy_VarCoords(r_VarAllModels(:, 0), r_MeanWarmingAllModels(:, 0))
                    ;
                    r_ConExtAllModels      := r_MeanWarmingAllModels(:, 2)
                end if
                r_VarAllModels(iModel, :)   = r_VarAll
                r_MeanWarmingAllModels(iModel, :) = r_MeanWarmingAll
                r_ConExtAllModels(iModel)   = r_ConExtAll
                delete([/r_VarAll/])
                delete([/r_MeanWarmingAll/])
                delete([/r_ConExtAll/])
            end if
        end do
        ; 
        ; printVarSummary(r_VarAllModels) ; [Models | 23] x [Clim_STD_Trend | 3]
        ; printVarSummary(r_MeanWarmingAllModels) ; [Models | 23] x [Terms | 3]
        ; print(sprintf("%4.2f", r_VarAllModels(:, 0)) + ", " + \
        ;       sprintf("%4.2f", r_VarAllModels(:, 1)) + ", " + \
        ;       sprintf("%4.2f", r_VarAllModels(:, 2)) + ", " + \
        ;       sprintf("%4.2f", r_MeanWarmingAllModels(:, 2)))
        ; exit()
        ; =============================================================
        ;   3. cal the two correlation coef and then output 
        if (isdefined("r_CorCoefAllVar") .eq. False) then
            i_dims                  := (/dimsizes(s_VarName), 3, 2/)
            r_CorCoefAllVar         := new(i_dims, float)
            r_CorCoefAllVar!0        = "Var"
            r_CorCoefAllVar&Var      = ispan(0, dimsizes(s_VarName)-1, 1)
            r_CorCoefAllVar@Var      = s_VarName
            copy_VarCoords(r_VarAllModels(0, :), r_CorCoefAllVar(0, :, 0))
            r_CorCoefAllVar!2        = "Coef_SigLev"
            r_CorCoefAllVar&Coef_SigLev = (/0, 1/)
            ;
            r_CorCoefAllVar_ForExt  := r_CorCoefAllVar
            ;
            r_RhoAllVar             := r_CorCoefAllVar
        end if
        do ivar2 = 0, 2
            r_Return                := my_CorCoef_SigLevel(\
                                                r_VarAllModels(:, ivar2), \
                                                r_MeanWarmingAllModels(:, 2), \
                                                False)
            r_CorCoefAllVar(ivar, ivar2, 0) = r_Return(0)
            r_CorCoefAllVar(ivar, ivar2, 1) = r_Return(2)
            ;
            r_Return                := my_CorCoef_SigLevel(\
                                                r_VarAllModels(:, ivar2), \
                                                r_ConExtAllModels, False)
            r_CorCoefAllVar_ForExt(ivar, ivar2, 0) = r_Return(0)
            r_CorCoefAllVar_ForExt(ivar, ivar2, 1) = r_Return(2)
            ;; spearman rank order correlation
            r_RhoAllVar(ivar, ivar2, 0)     = spcorr(r_VarAllModels(:, ivar2), \
                                                r_MeanWarmingAllModels(:, 2))
            r_RhoAllVar(ivar, ivar2, 1)     = rtest(r_RhoAllVar(ivar, ivar2, 0),\
                                    dimsizes(r_MeanWarmingAllModels(:, 2)), 0)
        end do
        print("-------------------------------------------------------")
        print("CorCoef: " + s_VarName(ivar) + ", in " + s_Exp(iExp))
        print("Mean CorCoef: " + sprintf("%4.2f", r_CorCoefAllVar(ivar,:,0)) + \
              " (p < " + sprintf("%4.2f", r_CorCoefAllVar(ivar, :, 1)) + ")" )
        print("Ext CorCoef: " + sprintf("%4.2f", r_CorCoefAllVar_ForExt(ivar,:,0)) + \
              " (p < " + sprintf("%4.2f", r_CorCoefAllVar_ForExt(ivar,:,1)) + ")" )
        ; print("Rho: " + r_RhoAllVar(ivar, :, 0) + " (p < " + \
        ;       sprintf("%4.2f", r_RhoAllVar(ivar, :, 1)) + ")" )
        print("-------------------------------------------------------")
        ;
    end do 
end do 
;

; =============================================================
;
end







