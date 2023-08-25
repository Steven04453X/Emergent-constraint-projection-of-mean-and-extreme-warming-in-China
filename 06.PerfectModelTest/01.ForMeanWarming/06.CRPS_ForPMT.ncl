; ===================================================
; Author: Ziming Chen 
; This code is for:
;   - cal the continous ranked probability score (CRPS)
;   - CRPS(P, x_a) = int([P(x) - P_a(x)]^2 dx), -inf to inf 
;   - P & P_a are cumulative distributions:
;   -   P(x)       = int(p(y) dy), -inf to x 
;   -   P_a(x)     = H(x - x_a)
;   -   H(x) = 0, for x < 0; 1, for x >= 0
;   1. indicate the pseudo model, selected models, Exp 
;   2. for each pseudo model, read the pseudo-obs warming at first, 
; and then read the raw projection 
;   3. read the constrained projection for all model except for the 
; pseudo-obs model 
;   4. cal the CRPS (self-defined function?)
;   5. output the results
; ===================================================
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
undef("my_CRPS")
function my_CRPS(x_a[*]:numeric, x[*]:numeric, PDF_Min, PDF_Max, DeltaBin)
;   - CRPS(P, x_a) = int([P(x) - P_a(x)]^2 dx), -inf to inf 
;   - P & P_a are cumulative distributions:
;   -   P(x)       = int(p(y) dy), -inf to x 
;   -   P_a(x)     = H(x - x_a)
;   -   H(x) = 0, for x < 0; 1, for x >= 0
local r_x_a, r_x 
begin
;
r_x_a                           = x_a
r_x                             = x 
;
;; PDF 
r_PDF_x                         = my_Historgram([/r_x/], PDF_Min, PDF_Max, \
                                                DeltaBin)
; print(r_PDF_x&bins + ": " + r_PDF_x)
;
r_P_a                           = r_PDF_x&bins
r_P_a                           = where(r_P_a .ge. r_x_a, 1., 0.)
;
r_P                             = r_PDF_x
r_P                             = 0
do ibins = 0, dimsizes(r_PDF_x)-1
    r_P(ibins)                  = sum(r_PDF_x(:ibins)) / 100.
end do
; print(r_PDF_x&bins + ": " + r_P + ", " + r_P_a)
; exit()
r_P                             = r_P - r_P_a
r_P                             = r_P^2
;
r_CRPS                          = dim_sum_n_Wrap(r_P, 0)
return(r_CRPS)
;
end
;
begin
; ===================================================
;   1. indicate the pseudo model, selected models, Exp 
if (isdefined("s_PseudoModel") .eq. False) then
    s_PseudoModel               = (/ "ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                         "BCC-CSM2-MR", "CAMS-CSM1-0", "CanESM5", \
                         "CanESM5-CanOE", "CNRM-CM6-1", "CNRM-ESM2-1", \
                         "EC-Earth3", "EC-Earth3-Veg", "FGOALS-f3-L", \
                         "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
                         "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8", \
                         "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR", \
                         "MRI-ESM2-0", "UKESM1-0-LL"/)
end if
;
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
if (isdefined("s_Exp") .eq. False) then
    s_Exp                       = (/"ssp126", "ssp245"/) ; (/"ssp370", "ssp585"/) ; 
end if
s_Terms                         = (/"Near", "Mid", "Long"/)
;
; ===================================================
;   2. for each pseudo model, read the pseudo-obs warming at first, 
; and then read the raw projection 
s_DirRead  = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/ProjectionOfAllReal/"
s_FileRead              = s_DirRead + "ProjectionInDifTerms_AllReal" + \
                            s_4Models + ".nc"
FileProjChanges         = addfile(s_FileRead, "r")
;
do iPse = 0, dimsizes(s_PseudoModel)-1
    print("Pseudo Model: " + s_PseudoModel(iPse))
    do iExp = 0, dimsizes(s_Exp)-1
        ;
        s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/06.PerfectModelTest/01.ForMeanWarming/"
        s_FileRead              = s_DirRead + \
                                "FurtherCorrectConstrainedWarming_PseModel_" + \
                                s_PseudoModel(iPse) + "_" + s_Exp(iExp) + \
                                s_4Models + ".nc"
        ls                     := systemfunc("ls " + s_FileRead)
        if (ismissing(ls(0))) then
            print(s_PseudoModel(iPse) + ", " + s_Exp(iExp) + " Missing!")
            continue
        end if
        FileCon                := addfile(s_FileRead, "r")
        ;
        s_VarName2              = "PseudoWarming_" + s_Exp(iExp)
        if (my_VarNameExist(FileCon, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue
        end if
        r_Pseudo               := FileCon->$s_VarName2$
        ;
        ;; Constrained 
        s_VarName2              = "ConstrainedWarming_" + s_Exp(iExp)
        r_ConWarming           := FileCon->$s_VarName2$
        ;
        ;; Raw projection 
        do iModel = 0, dimsizes(s_ModelSelected)-1
            l_list             := my_CMIP6Model_Exp_Path_ForAllReal(\
                                s_ModelSelected(iModel), s_Exp(iExp), "Amon")
            s_Real             := l_list[1]
            l_list_            := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                s_ModelSelected(iModel), s_Exp(iExp), "Amon")
            s_ExpModel_Path_   := l_list_[0]
            s_Real_            := l_list_[1]
            if (all(ismissing(s_Real)) .or. \
                dimsizes(s_Real_) .gt. dimsizes(s_Real) ) then
                s_Real         := s_Real_
            end if
            ;
            do iReal = 0, dimsizes(s_Real)-1
            ;
                s_VarNameProj  := "tas_inDifTerms_" + s_ModelSelected(iModel) +\
                                 "_" + s_Exp(iExp) + "_" + s_Real(iReal)
                if (s_ModelSelected(iModel) .ne. s_PseudoModel(iPse)) then
                    ;; Being constrained model 
                    if (my_VarNameExist(FileProjChanges, s_VarNameProj) .eq. \
                        False) then
                        ; print(s_VarNameProj + "  Missing!")
                        continue
                    end if
                    r_ProjTerm := FileProjChanges->$s_VarNameProj$
                ;
                    ;; store up 
                    if (isdefined("r_ProjTerm_AllReal") .eq. False .or. \
                        iReal .eq. 0) then
                        i_dims                     := array_append_record(\
                                    dimsizes(s_Real), dimsizes(r_ProjTerm), 0)
                        r_ProjTerm_AllReal         := new(i_dims, float)
                        copy_VarCoords(r_ProjTerm, r_ProjTerm_AllReal(0, :,:,:))
                        r_ProjTerm_AllReal!0        = "Real"
                        r_ProjTerm_AllReal&Real     = s_Real
                    end if
                    r_ProjTerm_AllReal(iReal,:,:,:) = r_ProjTerm
                end if 
                ;
            end do 
            ;
            if (s_ModelSelected(iModel) .ne. s_PseudoModel(iPse)) then 
                if (isdefined("r_ProjTerm_AllReal") .eq. False) then
                    print("r_ProjTerm_AllReal Missing! " + \
                          s_ModelSelected(iModel) + " " + s_Exp(iExp))
                    continue
                end if
                r_ProjTerm_AllReal := dim_avg_n_Wrap(r_ProjTerm_AllReal, 0)
            ;
                ; ===================================================
                ;   3. read the constrained projection for all model except for the 
                ; pseudo-obs model 
                s_VarName2 := "ConstrainedWarming_" + s_ModelSelected(iModel)\
                                + "_" + s_Exp(iExp)
                r_TempCon  := FileCon->$s_VarName2$
            ;
                if (isdefined("r_ConWarmingAllModel") .eq. False .or. \
                    iReal .eq. 0) then
                    i_dims := array_append_record(dimsizes(s_ModelSelected), \
                                                      dimsizes(r_TempCon), 0)
                    r_ConWarmingAllModel := new(i_dims, float)
                    copy_VarCoords(r_TempCon, r_ConWarmingAllModel(0, :, :, :))
                    r_ConWarmingAllModel!0= "Model"
                    r_ConWarmingAllModel&Model = s_ModelSelected
                ;
                    r_ProjTerm_AllModel  := r_ConWarmingAllModel
                end if
                r_ConWarmingAllModel(iModel, :, :, :) = r_TempCon
                r_ProjTerm_AllModel(iModel, :, :, :)  = r_ProjTerm_AllReal
            ;
                delete([/r_ProjTerm_AllReal/])
            end if 
        end do 
        ;
        if (isdefined("r_ProjTerm_AllModel") .eq. False) then
            print("r_ProjTerm_AllModel, " + s_ModelSelected(iModel) + ", " + \
                  s_Exp(iExp))
            continue
        end if
        ;; China Area 
        r_ProjTerm_AllModel          := my_SelChina([/r_ProjTerm_AllModel/])
        r_ConWarmingAllModel         := my_SelChina([/r_ConWarmingAllModel/])
        r_ConWarming                 := my_SelChina([/r_ConWarming/])
        r_Pseudo                     := my_SelChina([/r_Pseudo/])
        ;
        ;; Check the AreaMean 
        ; r_ConWarming_AreaMean        := my_AreaMean2([/r_ConWarming/], \
        ;                             r_ConWarming&lat, r_ConWarming&lon, 0)
        ; r_ConWarmingAllModel_AreaMean:= my_AreaMean2([/r_ConWarmingAllModel/], \
        ;                 r_ConWarmingAllModel&lat, r_ConWarmingAllModel&lon, 0)
        ; print(s_PseudoModel(iPse) + "")
        ; print(r_ConWarming_AreaMean + ", " + dim_avg_n(r_ConWarmingAllModel_AreaMean, 0))
        ; print("")
        ;
        ; =========================================================
        ;   4. cal the CRPS (self-defined function?)
        r_ProjTerm_CRPS              := r_ProjTerm_AllModel(0, :, :, :)
        r_ConWarmingAllModel_CRPS    := r_ProjTerm_AllModel(0, :, :, :)
        i_dims                       := dimsizes(r_ProjTerm_AllModel(0, :, :,:))
        do iTerm = 0, i_dims(0)-1
            do ilat = 0, i_dims(1)-1
                do ilon = 0, i_dims(2)-1
                    if (all(ismissing(r_ProjTerm_AllModel(:, iTerm, \
                                                          ilat, ilon))) ) then
                        continue
                    end if
                    r_ProjTerm_CRPS(iTerm, ilat, ilon) = my_CRPS(\
                                    r_Pseudo(iTerm, ilat, ilon),\
                                    r_ProjTerm_AllModel(:, iTerm, ilat, ilon), \
                                    0, 6, 0.1)
                    r_ConWarmingAllModel_CRPS(iTerm, ilat, ilon) = my_CRPS(\
                                    r_Pseudo(iTerm, ilat, ilon),\
                                    r_ConWarmingAllModel(:, iTerm, ilat, ilon),\
                                    0, 6, 0.1)
                    ; print(r_ProjTerm_CRPS(iTerm, ilat, ilon) + ", " + \
                    ;       r_ConWarmingAllModel_CRPS(iTerm, ilat, ilon))
                    ; exit()
                end do
            end do
        end do
        r_ProjTerm_CRPS             := my_AreaMean2([/r_ProjTerm_CRPS/], \
                                                    r_ProjTerm_CRPS&lat, \
                                                    r_ProjTerm_CRPS&lon, 0)
        r_ConWarmingAllModel_CRPS   := my_AreaMean2([/r_ConWarmingAllModel_CRPS/], \
                                            r_ConWarmingAllModel_CRPS&lat, \
                                            r_ConWarmingAllModel_CRPS&lon, 0)
        print("pseudo-obs: " + s_PseudoModel(iPse) + "")
        print("CRPS: Con, Raw")
        print(r_ConWarmingAllModel_CRPS + ", " + r_ProjTerm_CRPS)
        print("")
        ;
        delete([/r_ConWarmingAllModel/])
        ;
        ; ===================================================
        ;   5. output the results
        s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/06.PerfectModelTest/01.ForMeanWarming/"
        system("mkdir -p " + s_DirOut)
        s_FileOut               = s_DirOut + "CRPS_ForPMT_AllPseModel_" + \
                                s_Exp(iExp) + s_4Models + ".nc"
        ls                      = systemfunc("ls " + s_FileOut)
        if (ismissing(ls(0))) then
            FileOut            := addfile(s_FileOut, "c")
            GAttr               = True
            GAttr@history       = systemfunc("date") + "  " + systemfunc("pwd")\
                                 + "/06.CRPS_ForPMT.ncl"
            fileattdef(FileOut, GAttr)
    ;
        else 
            FileOut            := addfile(s_FileOut, "w")
        end if
        ;
        s_VarName2             := "ConCRPS_ForPseudo_" + s_PseudoModel(iPse) + \
                                "_" + s_Exp(iExp)
        r_ConWarmingAllModel_CRPS@long_name = s_VarName2
        r_ConWarmingAllModel_CRPS@Values    = r_ConWarmingAllModel_CRPS
        FileOut->$s_VarName2$   = r_ConWarmingAllModel_CRPS
        ;
        s_VarName2             := "RawCRPS_ForPseudo_" + s_PseudoModel(iPse) + \
                                "_" + s_Exp(iExp)
        r_ProjTerm_CRPS@long_name = s_VarName2
        r_ProjTerm_CRPS@Values    = r_ProjTerm_CRPS
        FileOut->$s_VarName2$     = r_ProjTerm_CRPS
    end do
    ;
end do
;


;
end














