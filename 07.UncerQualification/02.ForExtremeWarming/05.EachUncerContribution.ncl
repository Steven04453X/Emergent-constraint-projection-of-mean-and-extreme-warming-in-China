; =======================================================
; Author: Ziming Chen 
; This code is for:
;   - cal the area-mean value of each uncertainty component
;   - estimate the remainder variance by using the total variance minusing
; the var related to EC and the var related to internal var 
;   1. indicate the models, exp
;   2. open the data files that keep the uncertainty, and then read
; each uncertainty var, and then cal the area mean over China 
;   3. read the constrained projection for each model, and then 
; cal the inter-model var 
;   4. cal the remainder var 
; =======================================================
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
; =======================================================
;   1. indicate the models, exp
s_ModelSelected         = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                        "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                        "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
s_4Models               = "_18Models"
s_4ModelsGSAT_Con       = "_23Models" ; Keep it as a largest ensemble
s_ModelLENS                     = (/"CanESM5", "MIROC6"/)
s_Exp                           = (/"ssp126", "ssp245", "ssp370", "ssp585"/)
;
s_DirOut = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/02.ForExtremeWarming/"
system("mkdir -p " + s_DirOut)
s_FileOut                       = s_DirOut + "EachUncerContribution.nc"
ls                              = systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileOut                     = addfile(s_FileOut, "c")
    GAttr                       = True
    GAttr@history               = systemfunc("date") + "  " + systemfunc("pwd")\
                                + "/05.EachUncerContribution.ncl"
    GAttr@Uncertainty           = "Var(ConstrainedAcrossModels) = " + \
                                "Var(EC) + Var(InterVar) + Var(Remainder)"
    fileattdef(FileOut, GAttr)
else
    FileOut                     = addfile(s_FileOut, "w")
end if
; =======================================================
;   2. open the data files that keep the uncertainty, and then read
; each uncertainty var, and then cal the area mean over China 
do iExp = 0, dimsizes(s_Exp)-1
    ;; Var Related to EC 
    s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/02.ForExtremeWarming/"
    s_FileRead                  = s_DirRead + \
            "Dis_of_ConstrainedLocalWarmingByLocalTrend_Yc2_" + s_Exp(iExp) + \
            s_4Models + ".nc"
    FileEC                     := addfile(s_FileRead, "r")
    s_VarName2                 := "STD_FurConByLocalTrend_Yc2_" + s_Exp(iExp)
    r_EC                       := FileEC->$s_VarName2$
    r_EC                        = r_EC^2
    ; for area-mena uncertainty 
    r_EC                       := my_AreaMean2([/r_EC/], r_EC&lat, r_EC&lon, 0)
    ;
    ;; Var Related to Internal Var 
    s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/02.ForExtremeWarming/04.ConEachRealForInterVarEst/"
    s_FileRead                  = s_DirRead + "EstInterVar_fromEachLENS_1000times.nc"
    FilePD_InterVar            := addfile(s_FileRead, "r")
    s_FileRead                  = s_DirRead + "InterVarInProj_" + s_Exp(iExp) +\
                            s_4Models + ".nc"
    FileProj_InterVar          := addfile(s_FileRead, "r")
    ;
    ;; Emergent constrained results 
    s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/05.GSAT_scaled_RegionalTXx/02.ConstrainAllReal/02.ConstraintByUsingLocalTrend/"
    s_FileRead                  = s_DirRead + \
                            "FurtherCorrectConstrainedTXxWarming_" + s_Exp(iExp) +\
                            s_4Models + ".nc"
    FileConstrarinedResults    := addfile(s_FileRead, "r")
    ;
    do iLENS = 0, dimsizes(s_ModelLENS)-1
        s_VarName2             := "InterVar_STD_ByLENS_" + s_ModelLENS(iLENS) + \
                            "_" + s_Exp(iExp)
        r_InterVar             := FilePD_InterVar->$s_VarName2$
        r_InterVar              = r_InterVar^2
        ;
        ;; proj IV 
        s_VarName2             := "InterVar_STD_Proj_ByLENS_" + \
                                s_ModelLENS(iLENS) + "_" + s_Exp(iExp)
        r_Proj_InterVar        := FilePD_InterVar->$s_VarName2$
        r_Proj_InterVar         = r_Proj_InterVar^2
        ;
        ;; Area Mean at first 
        r_InterVar             := my_SelChina([/r_InterVar/])
        r_InterVar             := my_AreaMean2([/r_InterVar/], r_InterVar&lat, \
                                               r_InterVar&lon, 0)
        r_Proj_InterVar        := my_SelChina([/r_Proj_InterVar/])
        r_Proj_InterVar        := my_AreaMean2([/r_Proj_InterVar/], r_Proj_InterVar&lat, \
                                               r_Proj_InterVar&lon, 0)
        if (iLENS .eq. 0 .or. isdefined("r_InterVarAll") .eq. False) then
            i_dims             := array_append_record(dimsizes(s_ModelLENS), \
                                                      dimsizes(r_InterVar), 0)
            r_InterVarAll      := new(i_dims, float)
            copy_VarCoords(r_InterVar, r_InterVarAll(0, :))
            r_InterVarAll!0     = "LENS"
            r_InterVarAll&LENS  = s_ModelLENS
            ;
            r_InterVarProjAll   = r_InterVarAll ;(:, :, 0, 0)
        end if
        r_InterVarAll(iLENS, :) = r_InterVar
        r_InterVarProjAll(iLENS, :) = r_Proj_InterVar
    end do
    r_PD_InterVarAllModels     := r_InterVarAll
    r_Proj_InterVarAllModels   := r_InterVarProjAll
    ; do iModel = 0, dimsizes(s_ModelSelected)-1
    ;     ; s_VarName2              = "ConstrainedWarming_STD_" + \
    ;     ;                         s_ModelSelected(iModel) + "_" + s_Exp(iExp) + \
    ;     ;                         "_AllReal"
    ;     ; r_PD_InterVar          := FilePD_InterVar->$s_VarName2$
    ;     ; r_PD_InterVar           = r_PD_InterVar^2
    ;     ; r_PD_InterVar          := my_SelChina([/r_PD_InterVar/])
    ;     ; r_PD_InterVar          := my_AreaMean2([/r_PD_InterVar/], \
    ;     ;                                        r_PD_InterVar&lat, \
    ;     ;                                        r_PD_InterVar&lon, 0)
    ;     ;
    ;     ;; Projected Inter Var 
    ;     s_VarName2             := "InterVarInProj_" + s_ModelSelected(iModel) + \
    ;                             "_" + s_Exp(iExp)
    ;     if (my_VarNameExist(FileProj_InterVar, s_VarName2) .eq. False) then
    ;         print(s_VarName2 + " Missing!")
    ;         continue
    ;     end if
    ;     r_Proj_InterVar        := FileProj_InterVar->$s_VarName2$
    ;     r_Proj_InterVar         = r_Proj_InterVar^2
    ;     r_Proj_InterVar        := my_SelChina([/r_Proj_InterVar/])
    ;     r_Proj_InterVar        := my_AreaMean2([/r_Proj_InterVar/], \
    ;                                            r_Proj_InterVar&lat, \
    ;                                            r_Proj_InterVar&lon, 0)
    ;     ;
    ;     ;
    ;     if (iModel .eq. 0 .or. \
    ;         isdefined("r_Proj_InterVarAllModels") .eq. False) then
    ;         i_dims             := array_append_record(dimsizes(s_ModelSelected), \
    ;                                             dimsizes(r_Proj_InterVar), 0)
    ;         r_Proj_InterVarAllModels     := new(i_dims, float)
    ;         ;
    ;     end if
    ;     ; ; for local-cal uncertainty
    ;     ; r_PD_InterVarAllModels(iModel, :, :, :) = r_PD_InterVar
    ;     ; r_Proj_InterVarAllModels(iModel, :, :, :)= r_Proj_InterVar
    ;     ;
    ;     ; for area-mena uncertainty 
    ;     ; r_PD_InterVarAllModels(iModel, :) = r_PD_InterVar
    ;     r_Proj_InterVarAllModels(iModel, :)= r_Proj_InterVar
    ; end do
    ; ;
    ; ;; constrained values
    ; do iModel = 0, dimsizes(s_ModelSelected)-1
    ;     ; =======================================================
    ;     ;   3. read the constrained projection for each model, and then 
    ;     ; cal the inter-model var 
    ;     s_VarName2             := "ConstrainedWarming_" + \
    ;                             s_ModelSelected(iModel) + "_" + s_Exp(iExp)
    ;     r_ConProj              := FileConstrarinedResults->$s_VarName2$
    ;     ;
    ;     if (iModel .eq. 0 .or. isdefined("r_ConProjAllModels") .eq. False) then
    ;         i_dims             := array_append_record(dimsizes(s_ModelSelected), \
    ;                                                   dimsizes(r_ConProj), 0)
    ;         r_ConProjAllModels := new(i_dims, float)
    ;         copy_VarCoords(r_ConProj, r_ConProjAllModels(0, :, :, :))
    ;         ; copy_VarCoords(r_PD_InterVarAllModels(:, 0, 0, 0), \
    ;         ;                r_ConProjAllModels(:, 0, 0, 0))
    ;     end if
    ;     r_ConProjAllModels(iModel, :, :, :) = r_ConProj
    ; end do
    ;
    ;; spread across all models and realizations
    s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/02.ForExtremeWarming/04.ConEachRealForInterVarEst/"
    s_FileRead                 := s_DirRead + \
                            "FurtherCorrectConstrainedWarming_ForEachReal_" +\
                            s_Exp(iExp) + s_4Models + ".nc"
    FileConReal                := addfile(s_FileRead, "r")
    do iModel = 0, dimsizes(s_ModelSelected)-1
        ; =======================================================
        ;   3. read the constrained projection for each model, and then 
        ; cal the inter-model var 
        s_VarName2             := "ConstrainedWarming_" + \
                                s_ModelSelected(iModel) + "_" + s_Exp(iExp) + \
                                "_AllReal"
        if (my_VarNameExist(FileConReal, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue
        end if
        r_ConProj              := FileConReal->$s_VarName2$
        r_ConProj!0             = "Real"
        ;
        ;; all member together 
        if (iModel .eq. 0 .or. isdefined("r_ConProjAllModels") .eq. False) then
            r_ConProjAllModels := r_ConProj
        else 
            r_ConProjAllModels := array_append_record(r_ConProjAllModels, \
                                                      r_ConProj, 0)
        end if
        ;
        ; ;; inter-member ensemble at first 
        ; if (iModel .eq. 0 .or. isdefined("r_ConProjAllModels") .eq. False) then
        ;     i_dims             := array_append_record(dimsizes(s_ModelSelected), \
        ;                                     dimsizes(r_ConProj(0, :, :, :)), 0)
        ;     r_ConProjAllModels := new(i_dims, float)
        ;     copy_VarCoords(r_ConProj(0, :, :, :), r_ConProjAllModels(0, :, :, :))
        ;     ; copy_VarCoords(r_PD_InterVarAllModels(:, 0, 0, 0), \
        ;     ;                r_ConProjAllModels(:, 0, 0, 0))
        ; end if
        ; r_ConProjAllModels(iModel, :, :, :) = dim_avg_n_Wrap(r_ConProj, 0)
    end do
    ; ; r_ConProjAllModels&Real    := ispan(1, dimsizes(r_ConProjAllModels(:,0,0,0)), 1)
    ;
    r_PD_InterVarAllModels     := dim_avg_n_Wrap(r_PD_InterVarAllModels, 0)
    r_Proj_InterVarAllModels   := dim_avg_n_Wrap(r_Proj_InterVarAllModels, 0)
    ;; inter-model var in constrained projection 
    ; printMinMax(r_ConProjAllModels, False)
    ; print(r_ConProjAllModels(:, 2, 10, 10) + "")
    ; exit()
    ;
    r_ConProjAllModels         := my_SelChina([/r_ConProjAllModels/])
    ; s_FileTest = "./Test.nc"
    ; system("rm -rf " + s_FileTest)
    ; FileTest   = addfile(s_FileTest, "c")
    ; FileTest->Test = r_ConProjAllModels(0, :, :, :)
    ; exit()
    ;
    r_Temp                     := my_AreaMean2([/r_ConProjAllModels/], \
                                                r_ConProjAllModels&lat, \
                                                r_ConProjAllModels&lon, 0)
    ; printVarSummary(r_Temp)
    ; print(s_ModelSelected + ": " + r_Temp(:, 2))
    ; print("")
    ; print("STD: " + stddev(r_Temp(:, 2)))
    ; r_Temp                     := r_Temp(:, 2)
    ; r_TempMME                  := avg(r_Temp)
    ; i_Ind                      := ind(ismissing(r_Temp) .eq. False)
    ; r_Temp                     := my_QuaterNum2(r_Temp(i_Ind), 5, 95)
    ; print("MME: " + r_TempMME)
    ; print(r_Temp + "")
    ; exit()
    ;
    ; printMinMax(r_ConProjAllModels, False)
    r_ConProjAllModels         := dim_stddev_n_Wrap(r_ConProjAllModels, 0)
    r_ConProjAllModels          = r_ConProjAllModels^2
    r_ConProjAllModels         := my_SelChina([/r_ConProjAllModels/])
    ;
    r_ConProjAllModels         := my_AreaMean2([/r_ConProjAllModels/], \
                                               r_ConProjAllModels&lat, \
                                               r_ConProjAllModels&lon, 0)
    ;
    ;; remainder
    ; =========================================================
    ;   4. cal the remainder var 
    r_RemainderVar             := r_ConProjAllModels
    ; r_RemainderVar              = r_ConProjAllModels - r_EC - \
    ;                         r_PD_InterVarAllModels ;- r_Proj_InterVarAllModels
    r_RemainderVar              = r_ConProjAllModels - r_EC - \
                            r_PD_InterVarAllModels - r_Proj_InterVarAllModels
    ;
    ; r_RemainderVar             := my_AreaMean2([/r_RemainderVar/], \
    ;                                            r_RemainderVar&lat, \
    ;                                            r_RemainderVar&lon, 0)
    ; print(sprintf("%4.2f", r_ConProjAllModels) + ", " + \
    ;       sprintf("%4.2f", r_RemainderVar))
    ; printMinMax(r_ConProjAllModels, False)
    ;
    ; exit()
    print(s_Exp(iExp) + "")
    print("Tot Var, EC Var, PD_IV, Proj_IV, Remainder:")
    print(sprintf("%4.2f", r_ConProjAllModels) + ", " + sprintf("%4.2f", r_EC)\
          + ", " + sprintf("%4.2f", r_PD_InterVarAllModels) + ", " + \
          sprintf("%4.2f", r_Proj_InterVarAllModels) + ", " + \
          sprintf("%4.2f", r_RemainderVar))
    ; print("Tot Var, EC Var, PD_IV, Remainder:")
    ; print(sprintf("%4.2f", r_ConProjAllModels) + ", " + sprintf("%4.2f", r_EC)\
    ;       + ", " + sprintf("%4.2f", r_PD_InterVarAllModels) + ", " + \
    ;       sprintf("%4.2f", r_RemainderVar))
    print("")
    ;
    s_VarName2                  = "TotalUncer_" + s_Exp(iExp)
    r_ConProjAllModels@long_name= "Var of constrained projection across models"
    r_ConProjAllModels@Values   = r_ConProjAllModels
    FileOut->$s_VarName2$       = r_ConProjAllModels
    ;
    s_VarName2                  = "EC_Uncer_" + s_Exp(iExp)
    r_EC@long_name              = "Var of emergent constraint"
    r_EC@Values                 = r_EC
    FileOut->$s_VarName2$       = r_EC
    ;
    s_VarName2                  = "IV_Uncer_" + s_Exp(iExp)
    r_PD_InterVarAllModels@long_name     = "Var of internal variability"
    r_PD_InterVarAllModels@Values        = r_PD_InterVarAllModels
    FileOut->$s_VarName2$       = r_PD_InterVarAllModels
    ;
    s_VarName2                  = "IV_UncerInProj_" + s_Exp(iExp)
    r_Proj_InterVarAllModels@long_name   = "Var of internal variability in the projection"
    r_Proj_InterVarAllModels@Values      = r_Proj_InterVarAllModels
    FileOut->$s_VarName2$       = r_Proj_InterVarAllModels
    ;
    s_VarName2                  = "Re_Uncer_" + s_Exp(iExp)
    r_RemainderVar@long_name    = "Var of remainder"
    r_RemainderVar@Values       = r_RemainderVar
    FileOut->$s_VarName2$       = r_RemainderVar
    ;
    delete([/r_PD_InterVarAllModels, r_Proj_InterVarAllModels, r_ConProjAllModels/])
end do

;
end






