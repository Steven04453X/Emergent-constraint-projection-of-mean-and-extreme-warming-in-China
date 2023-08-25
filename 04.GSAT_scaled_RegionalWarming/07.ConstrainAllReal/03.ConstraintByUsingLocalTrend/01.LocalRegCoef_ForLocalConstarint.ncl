; =================================================================
; Author: Ziming Chen 
; This code is for:
;   - constraining the future local warming by using the regional trend 
;   - cal the regression and yintercept 
;   1. indicate the models, exp, terms, var name
;   2. for each model, read the present-day regional warming trend 
;   3. read the projected regional warmings in all of the realizations
;   4. make the inter-member ensemble, and then cal the regression coef and y intercept
;   5. output the regression coef and y intercept
; =================================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code/my_FunctionInPNNL.ncl"
;
begin
;
; ================================================================= 
;   1. indicate the models, exp, terms, var name
; s_ModelSelect           = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
;                             "BCC-CSM2-MR", "CAMS-CSM1-0", "CESM2-WACCM", \
;                             "CNRM-CM6-1", "CNRM-ESM2-1", "CanESM5", \
;                             "CanESM5-CanOE", "EC-Earth3", "EC-Earth3-Veg", \
;                             "FGOALS-f3-L", "FGOALS-g3", "GFDL-CM4", \
;                             "GFDL-ESM4", "GISS-E2-1-G", "HadGEM3-GC31-LL", \
;                             "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", \
;                             "MCM-UA-1-0", "MIROC6", "MIROC-ES2L", \
;                             "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", \
;                             "NESM3", "NorESM2-LM", "UKESM1-0-LL"/)
s_ModelSelect           = (/"ACCESS-CM2", "CNRM-CM6-1", "EC-Earth3", \
                            "GFDL-CM4", "MPI-ESM1-2-LR"/)
s_4Models               = "_5Models"
;
s_Exp                   = (/"ssp126", "ssp245", "ssp370", "ssp585"/)
i_yrBaseState           = (/1995, 2014/)
i_yrConstrained         = (/1981, 2014/)
i_yrProj                = (/ (/2021, 2040/), (/2041, 2060/), (/2080, 2099/) /)
s_yrProj                = (/"near", "mid", "long"/)
;
s_VarName               = "tas"
;
s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/03.ConstraintByUsingLocalTrend/"
system("mkdir -p " + s_DirOut)
s_FileOut               = s_DirOut + "LocalRegCoef_ForLocalConstraint" + \
                        s_4Models + ".nc"
system("rm -rf " + s_FileOut)
FileOut                 = addfile(s_FileOut, "c")
GAttr                   = True
GAttr@history           = systemfunc("date") + "  " + systemfunc("pwd") + \
                        "/01.LocalRegCoef_ForLocalConstraint.ncl"
; printVarSummary(i_yrProj)
GAttr@Terms             = s_yrProj + ": " + i_yrProj(:, 0) + " to " + \
                        i_yrProj(:, 1)
fileattdef(FileOut, GAttr)
;
;; Obs Warming Trend 
s_DirRead2 = "/global/cfs/cdirs/m1867/zmchen/my_data/Obs/TEMP/POST/YearMon/"
s_FileRead              = s_DirRead2 + \
                "SurfaceAirTemp_Post_YearMon_1900-2014.nc"
FileObs                 = addfile(s_FileRead, "r")
;
s_DirGPCP               = "$WACCEM/my_data/Obs/Precip/GPCP-v2.3/"
s_FileGPCP              = s_DirGPCP + "precip.mon.mean.nc"
FileGPCP                = addfile(s_FileGPCP, "r")
r_latGPCP               = FileGPCP->lat({-90:90})
r_lonGPCP               = FileGPCP->lon({0:360})
;
; s_VarNameObs                    = (/ "Cowtan_WayV2", "GISTEMP" /)
s_VarNameObs            = (/ "Cowtan_WayV2", "BEST", \
                                     "GISTEMP", "NOAAGlobalTemp" /)
; r_TrendObs              = new(dimsizes(s_VarNameObs), float)
do iObs = 0, dimsizes(s_VarNameObs)-1
    s_VarName2          = "tas_" + s_VarNameObs(iObs)
    r_HadCRUT4         := FileObs->$s_VarName2$(\
                                {i_yrConstrained(0):i_yrConstrained(1)}, :,:,:)
    r_HadCRUT4         := dim_avg_n_Wrap(r_HadCRUT4, 1)
    r_HadCRUT4!1        = "lat"
    r_HadCRUT4!2        = "lon"
    ;
    r_Temp_Trend        = trend_manken(r_HadCRUT4, False, 0)
    copy_VarCoords(r_HadCRUT4(0, :, :), r_Temp_Trend(0, :, :))
    ;
    ;; interpolation 
    opt                  = True
    opt@SrcGridLon      := r_Temp_Trend&lon
    opt@SrcGridLat      := r_Temp_Trend&lat
    opt@DstGridLon       = r_lonGPCP
    opt@DstGridLat       = r_latGPCP
    opt@DstReginal       = True
    opt@InterpMethod     = "conserve"
    opt@ForceOverwrite   = True
    r_Temp_Trend        := ESMF_regrid(r_Temp_Trend, opt)
    system("rm -rf ./destination_grid_file.nc ./PET0.RegridWeightGen.Log ./source_grid_file.nc ./weights_file.nc")

    if (isdefined("r_TrendObs") .eq. False .or. iObs .eq. 0) then
        i_dims          := array_append_record(dimsizes(s_VarNameObs)+1, \
                                            dimsizes(r_Temp_Trend(0, :, :)), 0)
        r_TrendObs       = new(i_dims, float)
        copy_VarCoords(r_Temp_Trend(0, :, :), r_TrendObs(0, :, :)) 
    end if
    r_TrendObs(iObs, :, :)= r_Temp_Trend(1, :, :)
end do 
;
;; CN05.1 
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/Obs/CN05.1/"
s_FileRead                = s_DirRead + "CN05.1_Tm_1961_2015_month_025x025.nc"
FileCN05                  = addfile(s_FileRead, "r")
i_Ind                     = my_TimeInd2(FileCN05->time, i_yrConstrained)
r_tm                      = FileCN05->tm(i_Ind(0):i_Ind(1), :, :)
r_tm                     := my_month_to_annual([/r_tm/])
r_tm_Trend                = trend_manken(r_tm, False, 0)
copy_VarCoords(r_tm(0, :, :), r_tm_Trend(0, :, :))
r_tm_Trend               := r_tm_Trend(1, :, :)
;
;; put it in the data pool 
r_lat_CN                  = r_tm&lat 
r_lon_CN                  = r_tm&lon 
r_Temp                   := r_TrendObs(0, {min(r_lat_CN):max(r_lat_CN)}, \
                                       {min(r_lon_CN):max(r_lon_CN)})
r_Temp                    = r_tm_Trend({r_Temp&lat}, {r_Temp&lon})
r_TrendObs(4, {r_Temp&lat}, {r_Temp&lon}) = r_Temp
;
r_TrendObs                = r_TrendObs * 10.
; print(r_TrendObs + "")
; r_TrendSTD                      = stddev(r_TrendObs)
; print("Obs STD: " + r_TrendSTD)
; exit()
r_TrendObs                     := dim_avg_n_Wrap(r_TrendObs, 0)
r_TrendObs@units                = "K per decade"
;
; =================================================================
;   3. read the projected regional warmings in all of the realizations
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/GSAT_AndRegionalTrend_AllReal/"
s_FileRead              = s_DirRead + "GSAT_AndRegionalTrend_" + \
                i_yrConstrained(0) + "-" + i_yrConstrained(1) + s_4Models + ".nc"
FileRegionalTrend_PD    = addfile(s_FileRead, "r")
;
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/ProjectionOfAllReal/"
s_FileRead              = s_DirRead + "ProjectionInDifTerms_AllReal" + \
                        s_4Models + ".nc"
FileProjWarming         = addfile(s_FileRead, "r")
do iExp = 0, dimsizes(s_Exp)-1
    do iModel = 0, dimsizes(s_ModelSelect)-1
        l_list         := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "Amon")
        s_Real         := l_list[1]
        l_list_        := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "Amon")
        s_Real_        := l_list_[1]
        if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real)) then
            s_Real     := s_Real_
        end if
        ;
        do iReal = 0, dimsizes(s_Real)-1
            s_VarName2 := "tas_inDifTerms_" + s_ModelSelect(iModel) + "_" + \
                        s_Exp(iExp) + "_" + s_Real(iReal)
            ; print(my_VarNameExist(FileProjWarming, s_VarName2))
            if (my_VarNameExist(FileProjWarming, s_VarName2) .eq. False) then
                print(s_VarName2 + " Missing!")
                continue
            end if
            r_ProjWarm := FileProjWarming->$s_VarName2$(:, {-90:90}, {0:360})
            ;
            ; =================================================================
            ;   2. for each model, read the present-day regional warming trend 
            s_VarName2 := "TrendPattern_" + s_ModelSelect(iModel) + \
                        "_historical_" + s_Real(iReal)
            if (my_VarNameExist(FileRegionalTrend_PD, s_VarName2) .eq. False) then
                print(s_VarName2 + " Missing!")
                continue
            end if
            r_PD_RegionalTrend = FileRegionalTrend_PD->$s_VarName2$(1, :, :)
            ;
            if (isdefined("r_PD_RegionalTrend_AllReal") .eq. False .or. \
                iReal .eq. 0) then
                i_dims        := array_append_record(dimsizes(s_Real), \
                                                dimsizes(r_PD_RegionalTrend), 0)
                r_PD_RegionalTrend_AllReal := new(i_dims, float)
                copy_VarCoords(r_PD_RegionalTrend, \
                               r_PD_RegionalTrend_AllReal(0, :, :))
                r_PD_RegionalTrend_AllReal!0= "Real"
                r_PD_RegionalTrend_AllReal&Real = s_Real
                ;
                i_dims        := array_append_record(dimsizes(s_Real), \
                                                     dimsizes(r_ProjWarm), 0)
                r_ProjWarmAllReal := new(i_dims, float)
                copy_VarCoords(r_ProjWarm, r_ProjWarmAllReal(0, :, :, :))
                copy_VarCoords(r_PD_RegionalTrend_AllReal(:, 0, 0), \
                               r_ProjWarmAllReal(:, 0, 0, 0))
            end if
            r_PD_RegionalTrend_AllReal(iReal, :, :) = r_PD_RegionalTrend
            r_ProjWarmAllReal(iReal, :, :, :)       = r_ProjWarm
        end do
        ;
        if (isdefined("r_PD_RegionalTrend_AllReal") .eq. False .or. \
            all(ismissing(r_PD_RegionalTrend_AllReal))) then
            print(s_ModelSelect(iModel) + " " + s_Exp(iExp) + "  Missing!")
            continue
        end if
        ;
        ; ====================================================================
        ;   4. make the inter-member ensemble, and then cal the regression coef and y intercept
        if (isdefined("r_ProjWarm_IME") .eq. False .or. iModel .eq. 0) then
            i_dims            := array_append_record(dimsizes(s_ModelSelect), \
                                                     dimsizes(r_ProjWarm), 0)
            r_ProjWarm_IME    := new(i_dims, float)
            copy_VarCoords(r_ProjWarm, r_ProjWarm_IME(0, :, :, :))
            r_ProjWarm_IME!0   = "Models"
            r_ProjWarm_IME&Models = s_ModelSelect
            ;
            i_dims            := array_append_record(dimsizes(s_ModelSelect), \
                                                dimsizes(r_PD_RegionalTrend), 0)
            r_PD_RegionalTrend_IME := new(i_dims, float)
            copy_VarCoords(r_PD_RegionalTrend, r_PD_RegionalTrend_IME(0, :, :))
            copy_VarCoords(r_ProjWarm_IME(:, 0, 0, 0), \
                           r_PD_RegionalTrend_IME(:, 0, 0))
        end if
        r_PD_RegionalTrend_IME(iModel, :, :) = dim_avg_n_Wrap(\
                                                r_PD_RegionalTrend_AllReal, 0)
        r_ProjWarm_IME(iModel, :, :, :)      = dim_avg_n_Wrap(\
                                                        r_ProjWarmAllReal, 0)
        ;
        delete([/r_PD_RegionalTrend_AllReal, r_ProjWarmAllReal/])
    end do
    ; ===============================================================
    ;   5. output the regression coef and y intercept
    r_RegAllTerms             := r_ProjWarm_IME(0, :, :, :)
    r_SigLevAllTerms          := r_ProjWarm_IME(0, :, :, :)
    r_yinterceptAllTerms      := r_ProjWarm_IME(0, :, :, :)
    do iTerm = 0, dimsizes(s_yrProj)-1
        r_Reg                 := regCoef_n(r_PD_RegionalTrend_IME, \
                                           r_ProjWarm_IME(:, iTerm, :, :), 0, 0)
        copy_VarCoords(r_PD_RegionalTrend_IME(0, :, :), r_Reg)
        r_tval                := r_Reg@tval
        r_tval                := reshape(r_tval, dimsizes(r_Reg))
        copy_VarCoords(r_Reg, r_tval)
        ; df                    := toint(r_tval)
        df                    := r_Reg@nptxy - 2
        df                    := reshape(df, dimsizes(r_Reg))
        copy_VarCoords(r_tval, df)
        r_Temp                := r_tval
        r_Temp                 = 0.5
        r_Sig                 := betainc(df/(df+r_tval^2), df/2.0, r_Temp)
        copy_VarCoords(r_Reg, r_Sig)
        ;
        r_yintercept          := r_Reg@yintercept 
        r_yintercept          := reshape(r_yintercept, dimsizes(r_Reg))
        copy_VarCoords(r_Reg, r_yintercept)
        ;
        r_RegAllTerms(iTerm, :, :) = r_Reg + 0
        r_SigLevAllTerms(iTerm,:,:)= r_Sig + 0
        r_yinterceptAllTerms(iTerm, :, :) = r_yintercept + 0
    end do
    s_VarName2                := "Reg_RegionalPD_Trend_AndWarming_" +s_Exp(iExp)
    r_RegAllTerms@long_name    = s_VarName2
    FileOut->$s_VarName2$      = r_RegAllTerms
    ;
    s_VarName2                := "SigLev_" + s_Exp(iExp)
    r_SigLevAllTerms@long_name = s_VarName2
    FileOut->$s_VarName2$      = r_SigLevAllTerms
    ;
    s_VarName2                 = "yintercept_" + s_Exp(iExp)
    r_yinterceptAllTerms@long_name = s_VarName2
    FileOut->$s_VarName2$      = r_yinterceptAllTerms
    ;
    print(s_Exp(iExp) + " Done!")
    ;
    ;; Scale 
    ; printVarSummary(r_TrendObs) ; [lat | 72] x [lon | 144]
    r_Temp_TrendObs           := conform(r_RegAllTerms, r_TrendObs,(/1,2/))
    copy_VarCoords(r_RegAllTerms, r_Temp_TrendObs)
    r_ConstrainedWarming_ByRegionalTrend := r_RegAllTerms
    r_ConstrainedWarming_ByRegionalTrend  = r_Temp_TrendObs * r_RegAllTerms + \
                                            r_yinterceptAllTerms
    s_VarName2                 = "ConstrainedRegionalWarming_MME_" + s_Exp(iExp)
    r_ConstrainedWarming_ByRegionalTrend@long_name = s_VarName2
    FileOut->$s_VarName2$      = r_ConstrainedWarming_ByRegionalTrend
    ;
    do iModel = 0, dimsizes(s_ModelSelect)-1
        r_delta_X             := r_PD_RegionalTrend_IME(iModel, :, :)
        r_delta_X              = r_PD_RegionalTrend_IME(iModel, :, :) - r_TrendObs
        ;
        r_Adjusted            := r_RegAllTerms
        r_delta_X             := conform(r_RegAllTerms, r_delta_X, (/1, 2/))
        copy_VarCoords(r_Adjusted, r_delta_X)
        r_Adjusted             = r_delta_X * r_Adjusted
        ;
        r_ConstrainedWarmingForModel = r_ProjWarm_IME(iModel,:,:,:) - r_Adjusted
        copy_VarCoords(r_Adjusted, r_ConstrainedWarmingForModel)
        ;
        s_VarName2             = "ConstrainedRegionalWarming_" + \
                                s_ModelSelect(iModel) + "_" + s_Exp(iExp)
        r_ConstrainedWarmingForModel@long_name = s_VarName2
        FileOut->$s_VarName2$  = r_ConstrainedWarmingForModel
    end do
    print(s_Exp(iExp) + "  Done!")
    ;
end do
;
end




















