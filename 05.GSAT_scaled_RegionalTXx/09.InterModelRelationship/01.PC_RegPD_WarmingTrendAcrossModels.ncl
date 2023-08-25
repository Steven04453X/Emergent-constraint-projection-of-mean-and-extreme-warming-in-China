; ======================================================================
; Author: Ziming Chen 
; This code is for:
;   - regressing the present-day local warming trend onto PC across models 
;   1. indicate the model, exp, var name of PD trend 
;   2. for each exp, read the PC at first, and then read the trend
;   3. make the inter-member ensemble for the PD trend if necessary 
;   4. regression 
;   5. output 
; ======================================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code//my_FunctionInPNNL.ncl"
;
undef("my_Reg_SigLev")
function my_Reg_SigLev(GSAT[*]:numeric, SSP[*][*][*]:numeric)
local r_GSAT2, r_SSP2
begin
;
r_GSAT2                       = GSAT 
r_SSP2                        = SSP 
;
        r_Reg                := regCoef_n(r_GSAT2, r_SSP2, 0, 0)
        copy_VarCoords(r_SSP2(0, :, :), r_Reg)
    ;
        r_tval               := r_Reg@tval 
        r_tval               := reshape(r_tval, dimsizes(r_Reg))
        copy_VarCoords(r_SSP2(0, :, :), r_tval)
        r_df                 := r_Reg@nptxy
        r_df                 := reshape(r_df, dimsizes(r_Reg))
        copy_VarCoords(r_tval, r_df)
    ;
        r_df                  = r_df - 2
        r_Temp               := r_tval
        r_Temp                = 0.5
        r_Sig                := betainc(r_df/(r_df + r_tval^2), r_df/2., r_Temp)
        copy_VarCoords(r_Temp, r_Sig)
    ;
        r_Temp               := r_Reg@yintercept
        r_Temp               := reshape(r_Temp, dimsizes(r_Reg))
        copy_VarCoords(r_Sig, r_Temp)
        ;
i_dims                        = array_append_record(2, dimsizes(r_Reg), 0)
r_Reg_All                     = new(i_dims, float)
copy_VarCoords(r_Reg, r_Reg_All(0, :, :))
r_Reg_All!0                   = "Reg_SigLev"
r_Reg_All&Reg_SigLev          = (/0, 1/)
r_Reg_All(0, :, :)            = r_Reg + 0
r_Reg_All(1, :, :)            = r_Sig + 0
;
r_Reg_All@long_name           = "regression btw GSAT and local warming"
r_Temp@long_name              = "y intercept of regression"
return([/r_Reg_All, r_Temp/])
;
end
;
begin
; ======================================================================
;   1. indicate the model, exp, var name of PD trend 
s_ModelSelected = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MRI-ESM2-0", "UKESM1-0-LL"/)
; s_ModelSelected = (/"ACCESS-CM2", "EC-Earth3", "CanESM5" /)
s_4Models       = "_17Models" ; "_3Models" ; 
s_Exp           = (/"ssp245"/)
s_VarName       = "TrendPattern"
;
s_RunAve        = "_20yr" ; or ""
; ======================================================================
;   2. for each exp, read the PC at first, and then read the trend
do iExp  = 0, dimsizes(s_Exp)-1
    s_DirRead   = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tasmax/POST/TXxRegGSAT/InterModelEOF/"
    s_FileRead  = s_DirRead + "Inter-modelEOF_Regression_" + s_Exp(iExp) + \
                    s_RunAve + s_4Models + ".nc"
    FilePC     := addfile(s_FileRead, "r")
    s_VarName2  = "RegInTimDim_PC_" + s_Exp(iExp)
    r_PC       := FilePC->$s_VarName2$(0, :, :)
    ;
    s_DirRead   = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/GSAT_AndRegionalTrend_AllReal/"
    s_FileRead  = s_DirRead + "GSAT_AndRegionalTrend_1981-2014_22Models.nc"
    FileTrend   = addfile(s_FileRead, "r")
    ;
    do iModel = 0, dimsizes(s_ModelSelected)-1
        ;
        l_list                  := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp(iExp), "Amon")
        s_ExpModel_Path         := l_list[0]
        s_Real                  := l_list[1]
        l_list_                 := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp(iExp), "Amon")
        s_ExpModel_Path_        := l_list_[0]
        s_Real_                 := l_list_[1]
        if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real) ) then
            s_Real              := s_Real_
            s_ExpModel_Path     := s_ExpModel_Path_
        end if
        ;
        do iReal = 0, dimsizes(s_Real)-1
            s_VarName2          := "TrendPattern_" + s_ModelSelected(iModel) + \
                                "_historical_" + s_Real(iReal)
            if (my_VarNameExist(FileTrend, s_VarName2) .eq. False) then
                print(s_VarName2 + "  Missing!")
                continue
            end if
            r_Trend             := FileTrend->$s_VarName2$(1, {-90:90}, {0:360}) 
            ;
            if (iReal .eq. 0 .or. isdefined("r_TrendAllReal") .eq. False) then
                i_dims          := array_append_record(dimsizes(s_Real), \
                                                       dimsizes(r_Trend), 0)
                r_TrendAllReal  := new(i_dims, float)
                copy_VarCoords(r_Trend, r_TrendAllReal(0, :, :))
            end if
            r_TrendAllReal(iReal, :, :) = r_Trend
        end do
        ; ============================================================
        ;   3. make the inter-member ensemble for the PD trend if necessary 
        r_TrendAllReal          := dim_avg_n_Wrap(r_TrendAllReal, 0)
        ;
        if (iModel .eq. 0 .or. isdefined("r_TrendAllModel") .eq. False) then
            i_dims              := array_append_record(dimsizes(s_ModelSelected), \
                                                       dimsizes(r_TrendAllReal), 0)
            r_TrendAllModel     := new(i_dims, float)
            copy_VarCoords(r_TrendAllReal, r_TrendAllModel(0, :, :))
            r_TrendAllModel!0    = "Models"
            r_TrendAllModel&Models = s_ModelSelected
        end if
        r_TrendAllModel(iModel, :, :) = r_TrendAllReal
        delete([/r_TrendAllReal/])
        ;
    end do
    ;
    ; ==============================================================
    ;   4. regression 
    i_dims                     := array_append_record(3, \
                                        dimsizes(r_TrendAllModel(0, :, :)), 0)
    i_dims                     := array_append_record(2, i_dims, 0)
    r_Reg_All                  := new(i_dims, float)
    copy_VarCoords(r_TrendAllModel(0, :, :), r_Reg_All(0, 0, :, :))
    ; printVarSummary(r_PC) ; [ProjectionTerm | 3] x [Models | 17]
    ; printVarSummary(r_Reg_All) ; [2] x [3] x [lat | 72] x [lon | 144]
    ;
    copy_VarCoords(r_PC(:, :), r_Reg_All(0, :, 0, 0))
    r_Reg_All!0                 = "RegAndSigLevel"
    r_Reg_All&RegAndSigLevel    = (/0, 1/)
    ;
    r_yintercept               := r_Reg_All(0, :, :, :)
    do iPeriod = 0, 2
        r_Temp                 := r_PC(iPeriod, :)
        r_Temp                  = dim_standardize_n_Wrap(r_Temp, 0, 0)
        ; print(r_PC(iPeriod, :) + ", " + r_Temp)
        ; exit()
        l_list                 := my_Reg_SigLev(r_Temp, \
                                                r_TrendAllModel(:, :, :))
        r_Reg_All(:, iPeriod, :, :) = l_list[0] + 0
        r_yintercept(iPeriod, :, :) = l_list[1] + 0
    end do
    ;
    ; =============================================================
    ;   5. output 
    s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/05.GSAT_scaled_RegionalTXx/09.InterModelRelationship/"
    system("mkdir -p " + s_DirOut)
    s_FileOut                   = s_DirOut + "TXxPC_RegPD_WarmingTrendAcrossModels_"\
                                + s_Exp(iExp) + s_RunAve + s_4Models + ".nc"
    system("rm -rf " + s_FileOut)
    FileOut                    := addfile(s_FileOut, "c")
    GAttr                       = True
    GAttr@history               = systemfunc("date") + "  " + systemfunc("pwd")\
                                + "/01.PC_RegPD_WarmingTrendAcrossModels.ncl"
    GAttr@Models                = s_ModelSelected
    fileattdef(FileOut, GAttr)
    s_VarName2                  = "PC_Reg_PD_WarmingTrend_" + s_Exp(iExp)
    r_Reg_All@long_name         = "TXx " + s_VarName2
    FileOut->$s_VarName2$       = r_Reg_All
    s_VarName2                  = "yintercept_" + s_Exp(iExp)
    r_Reg_All@long_name         = "TXx " + s_VarName2
    FileOut->$s_VarName2$       = r_yintercept
    ;
    print(s_Exp(iExp) + "  Done!")
    print("")
end do
;




end 






