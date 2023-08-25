; =======================================================
; Author: Ziming Chen 
; This code is for:
;   - making the distribution of pattern scaling coefficient 
;   - constrain the local warming by using the bootstrap members 
;   1. indicate the models, exp, projection terms
;   2. read the constrained GSAT, and the raw GSAT warming
;   3. read the local projected warming for each model
;   5. Select the data over China
;   6. for each regression coefficient, make 1000 synthetic members for the 
; regression coefficient btw local warming and GSAT warming by using bootstrap method
;   7. constrain the local warming, and output the PDF of constrained projection
; =======================================================
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
;   1. indicate the models, exp, projection terms
s_ModelSelect           = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                            "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                            "CNRM-ESM2-1", "CanESM5", "CanESM5-CanOE",\
                            "EC-Earth3", "EC-Earth3-Veg", "FGOALS-f3-L", \
                            "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
                            "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8",\
                            "IPSL-CM6A-LR", "MPI-ESM1-2-LR", "MRI-ESM2-0", \
                            "UKESM1-0-LL"/)
s_4Models               = "_22Models"
if (isdefined("s_Exp") .eq. False) then
    s_Exp               = "ssp245"
end if
i_yrBaseState           = (/1995, 2014/)
i_yrProj                = (/ (/2021, 2040/), (/2041, 2060/), (/2080, 2099/) /)
s_yrProj                = (/"near", "mid", "long"/)
s_VarName               = "tas"
;
; =======================================================
;   2. read the constrained GSAT, and the raw GSAT warming
s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/01.ForMeanWarming/"
s_FileRead              = s_DirRead + "Dis_of_pc_and_constrainedGSAT_" +\
                                s_Exp + s_4Models + ".nc"
FileConGSAT_PDF         = addfile(s_FileRead, "r")
;
s_DirRead  = "$WACCEM/my_data/CMIP6/tas/POST/GMSAT/MultiMembers/"
s_FileRead                = s_DirRead + "GMSAT_Series_MultiMember_"+s_Exp + ".nc"
FileProjGSAT              = addfile(s_FileRead, "r")
s_FileRead                = s_DirRead + "GMSAT_Series_MultiMember_historical.nc"
FileHistGSAT              = addfile(s_FileRead, "r")
;
s_DirRead  = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/ProjectionOfAllReal/"
s_FileRead                = s_DirRead + "ProjectionInDifTerms_AllReal" + \
                            s_4Models + ".nc"
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
        if (my_VarNameExist(FileHistGSAT, s_VarName2) .eq. False) then
            print(s_VarName2 + "  Missing!")
            continue
        end if;
        r_GSAT_Hist            := FileHistGSAT->$s_VarName2$(\
                                        {i_yrBaseState(0):i_yrBaseState(1)}, :)
        r_GSAT_Hist            := dim_avg_n_Wrap(r_GSAT_Hist, 1)
        r_GSAT_Hist            := dim_avg_n_Wrap(r_GSAT_Hist, 0)
        ;
        r_GSAT                  = r_GSAT - r_GSAT_Hist
        ;
        ; =======================================================
        ;   3. read the local projected warming for each model
        ;; tas projection 
        s_VarNameProj          := "tas_inDifTerms_" + s_ModelSelect(iModel) + \
                                    "_" + s_Exp + "_" + s_Real(iReal)
        if (my_VarNameExist(FileTasProjChanges, s_VarNameProj) .eq. False) then
            print(s_VarNameProj + "  Missing!")
            continue
        end if
        r_ProjTerm             := FileTasProjChanges->$s_VarNameProj$
        ;
        if (isdefined("r_ProjTerm") .eq. False .or. iReal .eq. 0) then
            i_dims             := array_append_record(dimsizes(s_Real), \
                                                      dimsizes(r_GSAT), 0)
            r_GSAT_AllReal     := new(i_dims, float)
            copy_VarCoords(r_GSAT, r_GSAT_AllReal(0, :))
            r_GSAT_AllReal!0    = "Real"
            r_GSAT_AllReal&Real = s_Real
            ;
            i_dims             := array_append_record(\
                                    dimsizes(s_Real), dimsizes(r_ProjTerm), 0)
            r_ProjTerm_AllReal := new(i_dims, float)
            copy_VarCoords(r_ProjTerm, r_ProjTerm_AllReal(0, :, :, :))
            copy_VarCoords(r_GSAT_AllReal(:, 0), r_ProjTerm_AllReal(:,0,0,0))
        end if
        r_GSAT_AllReal(iReal, :)        = r_GSAT
        r_ProjTerm_AllReal(iReal,:,:,:) = r_ProjTerm
    end do
    if (isdefined("r_GSAT_AllReal") .eq. False) then
        print("r_GSAT_AllReal Missing!")
        continue
    end if
    ; r_PD_GSAT_Trend_AllReal    := dim_avg_n_Wrap(r_PD_GSAT_Trend_AllReal, 0)
    r_GSAT_AllReal             := dim_avg_n_Wrap(r_GSAT_AllReal, 0)
    r_ProjTerm_AllReal         := dim_avg_n_Wrap(r_ProjTerm_AllReal, 0)
    ;
    ;; Constrained Warming 
    s_VarName2                  = "rcBoot_ConstrainedGSAT_" + s_Exp + "_" + \
                                    s_ModelSelect(iModel)
    r_ConstrainedGSAT          := FileConGSAT_PDF->$s_VarName2$
    ;
    if (isdefined("r_GSAT_All") .eq. False) then
        i_dims                 := array_append_record(dimsizes(s_ModelSelect), \
                                                      3, 0)
        r_GSAT_All             := new(i_dims, float)
        r_GSAT_All!0            = "Model"
        r_GSAT_All&Model        = s_ModelSelect
        r_GSAT_All!1            = "Terms"
        r_GSAT_All&Terms        = ispan(0, 2, 1)
        ;
        i_dims                 := array_append_record(dimsizes(s_ModelSelect), \
                                        dimsizes(r_ProjTerm_AllReal), 0)
        r_ProjTerm_All         := new(i_dims, float)
        copy_VarCoords(r_GSAT_All(:, 0), r_ProjTerm_All(:, 0, 0, 0))
        copy_VarCoords(r_ProjTerm_AllReal, r_ProjTerm_All(0, :, :, :))
        ;
        i_dims                 := array_append_record(dimsizes(s_ModelSelect), \
                                            dimsizes(r_ConstrainedGSAT), 0)
        r_ConstrainedGSAT_All   = new(i_dims, float)
        copy_VarCoords(r_ConstrainedGSAT, r_ConstrainedGSAT_All(0, :, :))
        copy_VarCoords(r_GSAT_All(:, 0), r_ConstrainedGSAT_All(:, 0, 0))
    end if
    r_GSAT_All(iModel, :)       = my_3TermsMean(r_GSAT_AllReal, i_yrProj)
    r_ProjTerm_All(iModel,:,:,:)= r_ProjTerm_AllReal
    r_ConstrainedGSAT_All(iModel, :, :) = r_ConstrainedGSAT
end do 
;
; =========================================================================
;   5. Select the data over China
r_ProjTerm_All                  := my_SelChina([/r_ProjTerm_All/])
; 
; =========================================================================
;   6. for each regression coefficient, make 1000 synthetic members for the 
; regression coefficient btw local warming and GSAT warming by using bootstrap method
; printVarSummary(r_ConstrainedGSAT_All) ; (model, boot, Terms)
;
s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/01.ForMeanWarming/"
system("mkdir -p " + s_DirOut)
setfileoption("nc", "Format", "LargeFile")
s_FileOut                        = s_DirOut + \
                        "Dis_of_constrainedLocalWarming_Yc1_" + s_Exp + \
                        s_4Models + ".nc"
system("rm -rf " + s_FileOut)
FileOut                          = addfile(s_FileOut, "c")
GAttr                            = True 
GAttr@history                    = systemfunc("date") + "  " + systemfunc("pwd")\
                        + "/02.Dis_of_pc1_and_constrainedLocalWarming_Yc1.ncl"
fileattdef(FileOut, GAttr)
;
nBoot                            = 100
nDim                             = 0
opt                              = True
opt@rseed3                       = "clock"
do iBoot = 0, dimsizes(r_ConstrainedGSAT_All(0, :, 0))-1
    ;
    print("Boot of GSAT: " + iBoot)
    i_dims                      := array_append_record(nBoot, \
                                                dimsizes(r_ProjTerm_All), 0)
    r_ConLocWarmingPattern      := new(i_dims, float)
    copy_VarCoords(r_ProjTerm_All, r_ConLocWarmingPattern(0, :, :, :, :))
    r_ConLocWarmingPattern!0     = "boot"
    r_ConLocWarmingPattern&boot  = ispan(0, nBoot-1, 1)
    ;
    i_dims                      := array_append_record(nBoot, \
                                                dimsizes(r_ProjTerm_All), 0)
    r_ConLocWarming             := new(i_dims, float)
    copy_VarCoords(r_ProjTerm_All, r_ConLocWarming(0, :, :, :, :))
    r_ConLocWarming!0            = "boot"
    r_ConLocWarming&boot         = ispan(0, nBoot-1, 1)
    ;
    do iTerm = 0, dimsizes(s_yrProj)-1
        print(s_yrProj(iTerm) + "")
        r_X                     := r_GSAT_All(:, iTerm)
        r_delta_X               := r_ConstrainedGSAT_All(:, iBoot, iTerm)
        r_delta_X                = r_GSAT_All(:, iTerm) - \
                                    r_ConstrainedGSAT_All(:, iBoot, iTerm)
        ;
        BootStrap               := bootstrap_regcoef(r_X, \
                                r_ProjTerm_All(:,iTerm, :, :), nBoot, nDim, opt)
        rcBoot                  := BootStrap[0] ; Bootstrapped regression coefficients in ascending order
        ; printVarSummary(rcBoot) ; [boot | 1000] x [lat | 72] x [lon | 144]
        ;
        do iModel = 0, dimsizes(s_ModelSelect)-1
            r_ConLocWarming     := rcBoot
            r_Temp              := conform(rcBoot, \
                                r_ProjTerm_All(iModel, iTerm, :, :), (/1, 2/))
            copy_VarCoords(rcBoot, r_Temp)
            r_ConLocWarming      = r_Temp - rcBoot * r_delta_X(iModel)
            r_ConLocWarmingPattern(:, iModel, iTerm, :, :) = r_ConLocWarming
        end do
        ; printMinMax(r_ConLocWarmingPattern, False)
        ; exit()
        ;;
    end do
    ;
    ; =========================================================================
    ;   7. constrain the local warming, and output the PDF of constrained projection
    do iModel = 0, dimsizes(s_ModelSelect)-1
        s_VarName2                := "ConstrainedWarming_Yc1_" + \
                                s_ModelSelect(iModel) + "_In_pc_boot" + iBoot
        FileOut->$s_VarName2$      = r_ConLocWarmingPattern(:, iModel, :, :, :)
    end do
end do
;
print(s_FileOut + "")
;
end



















