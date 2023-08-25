; ==================================================
; 本程序基于GSAT与局地升温的回归系数，对中国地区的升温进行约束
;   1. 读取回归系数，然后读取约束后的GSAT，两者相乘和加上截距
;   2. 读取约束前的结果，然后约束前后作差
; Attention: This code only show the projected warming in the constrained projection,
; and the raw projection, and the 5~95th as well
;
; Attention: This code cal the Low-likehood-high-impact results
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
begin
;
; ==================================================
;   1. 读取回归系数，然后读取约束后的GSAT，两者相乘和加上截距
s_4Models                 = "_18Models" ; "_5Models" ; "_30Models" ; "" ; 
; s_ModelSelected           = (/"ACCESS-CM2", "CNRM-CM6-1", "EC-Earth3", \
;                             "GFDL-CM4", "MPI-ESM1-2-LR"/)
s_ModelSelected           = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
s_Month                   = "AnMean"
s_Exp                     = "ssp585" ; "ssp370" ; "ssp126" ; "ssp245" ; 
;
if (s_Exp .eq. "ssp126") then
    s_ModelSelected      := (/"ACCESS-CM2", "AWI-CM-1-1-MR", \
                        "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                        "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
end if
i_yrProj                  = (/ (/2021, 2040/), (/2041, 2060/), \
                                   (/2080, 2099/) /)
s_yrProj                  = (/"near", "mid", "long"/)
;
; s_DirOut = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/03.ConstraintByUsingLocalTrend/"
; s_FileOut                 = s_DirOut + "ConstrainedWarmingInChina_" + s_Exp + \
;                             s_4Models + ".nc"
; system("rm -rf " + s_FileOut)
;     FileOut               = addfile(s_FileOut, "c")
;     GAttr                 = True
;     GAttr@history         = systemfunc("date") + "  " + systemfunc("pwd") + \
;                             "/07.1.1.ShowFurtherCorrectConstraint_withDif.ncl"
;     fileattdef(FileOut, GAttr)
;
; ;;;;; for Pattern 
; s_VarName                 = "tas" ; "pr" ; 
; s_DirRead = "$WACCEM/my_data/CMIP6/tas/POST/ProjectionOfAllReal/"
; s_FileRead                = s_DirRead + "ProjectionInDifTerms_AllReal" + \
;                             s_4Models + ".nc"
; FileRead_Raw              = addfile(s_FileRead, "r")
; ;;;;; for Pattern 
;
; s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/03.ConstraintByUsingLocalTrend/"
; s_FileRead                = s_DirRead + "FurtherCorrectConstrainedWarming_" + \
;                             s_Exp + s_4Models + ".nc"
; FileConstrained           = addfile(s_FileRead, "r")

; s_VarName2              = "ConstrainedWarming_" + s_Exp
; r_Warming               = FileConstrained->$s_VarName2$
; ;
; ;;; Area-Mean over China 
; r_China                := my_SelChina([/r_Warming/])
; ;
; r_China                := my_AreaMean2([/r_China/], r_China&lat, r_China&lon, 0)
; print("3 Term Warming:")
; print(r_China + " K")
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/02.ForExtremeWarming/04.ConEachRealForInterVarEst/"
s_FileRead              = s_DirRead + \
                                "FurtherCorrectConstrainedWarming_ForEachReal_"\
                                + s_Exp + s_4Models + ".nc"
FileConstrained         = addfile(s_FileRead, "r")
;
do iModel = 0, dimsizes(s_ModelSelected)-1
    s_VarName2          = "ConstrainedWarming_" + s_ModelSelected(iModel) + "_"\
                             + s_Exp + "_AllReal"
    L_Return            = my_VarNameExist(FileConstrained, s_VarName2)
    if (L_Return .eq. False) then
        print(s_VarName2 + " Missing!")
        continue
    end if
    r_Reg_Warming      := FileConstrained->$s_VarName2$
    r_Warming0         := r_Reg_Warming
    ;
    r_Reg_Warming      := my_SelChina([/r_Reg_Warming/])
    ;
    r_Reg_Warming      := my_AreaMean2([/r_Reg_Warming/], r_Reg_Warming&lat, \
                                       r_Reg_Warming&lon, 0)
    r_Reg_Warming      := dim_avg_n_Wrap(r_Reg_Warming, 0)
    if (isdefined("r_Reg_Warming_All") .eq. False) then
        r_Reg_Warming_All        = new((/dimsizes(s_ModelSelected), 3/), float)
        copy_VarCoords(r_Reg_Warming, r_Reg_Warming_All(0, :))
        r_Reg_Warming_All!0      = "Models"
        r_Reg_Warming_All&Models = s_ModelSelected
        ;
        i_dims                  := array_append_record(dimsizes(s_ModelSelected), \
                                        dimsizes(r_Warming0(0, :, :, :)), 0)
        r_Warming                = new(i_dims, float)
        copy_VarCoords(r_Warming0(0, :, :, :), r_Warming(0, :, :, :))
        r_Warming!0              = "Models"
        r_Warming&Models         = s_ModelSelected
    end if
    r_Reg_Warming_All(iModel, :) = r_Reg_Warming
    r_Warming(iModel, :, :, :)   = dim_avg_n_Wrap(r_Warming0, 0)
end do
print(r_Reg_Warming_All&Models + ": " + r_Reg_Warming_All(:, 2))
r_China                          = dim_avg_n_Wrap(r_Reg_Warming_All, 0)
r_Warming                       := dim_avg_n_Wrap(r_Warming, 0)
; r_Temp                 := my_QuaterNum_withMean3(r_Reg_Warming_All(:, 2), 5, 95)
; print("")
; print("MME: " + avg(r_Reg_Warming_All(:, 2)))
; print("5th to 95th:")
; print(r_Temp + "")
; exit()
r_Reg_Warming_Spread             = dim_stddev_n_Wrap(r_Reg_Warming_All, 0)
print("STD: " + r_Reg_Warming_Spread(2))
;;
; s_VarNameAll                     = getfilevarnames(FileRead_Raw)
s_VarName                        = "TXx" ; "pr" ; 
s_DirRead = "$WACCEM/my_data/CMIP6/tasmax/POST/TXx/AllReal/"
    do iModel = 0, dimsizes(s_ModelSelected)-1
        ;
        l_list                  := my_CMIP6Model_Exp_Path_ForAllReal(\
                                        s_ModelSelected(iModel), s_Exp, "day")
        s_Real                  := l_list[1]
        l_list_                 := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                        s_ModelSelected(iModel), s_Exp, "day")
        s_Real_                 := l_list_[1]
        if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real)) then
            s_Real              := s_Real_
        end if
        ;
        if (all(ismissing(s_Real)) ) then
            print(s_ModelSelected(iModel) + ", " + s_Exp + " Missing!")
            continue
        end if
        s_FileRead               = s_DirRead + s_VarName + "_" + s_Month + "_"+\
                                s_Exp + "_" + s_ModelSelected(iModel) + \
                                "_CMIP6_AllReal.nc"
        FileRaw                 := addfile(s_FileRead, "r")
        s_FileRead               = s_DirRead + s_VarName + "_" + s_Month + \
                                "_historical_" + s_ModelSelected(iModel) + \
                                "_CMIP6_AllReal.nc"
        FileRawHist             := addfile(s_FileRead, "r")
        do iReal = 0, dimsizes(s_Real)-1
            s_VarName2          := s_VarName + "_" + s_ModelSelected(iModel) + \
                                "_" + s_Exp + "_" + s_Real(iReal)
            if (my_VarNameExist(FileRaw, s_VarName2) .eq. False) then
                print(s_VarName2 + " Missing!")
                continue 
            end if
            r_Temp              := FileRaw->$s_VarName2$(:, {-90:90}, {0:360})
            ;
            s_VarName2          := s_VarName + "_" + s_ModelSelected(iModel) + \
                                    "_historical_" + s_Real(iReal)
            if (my_VarNameExist(FileRawHist, s_VarName2)  .eq. False) then
                print(s_VarName2 + " Missing!")
                continue 
            end if
            r_TempHist          := FileRawHist->$s_VarName2$({1995:2014}, \
                                                             {-90:90}, {0:360})
            r_TempHist          := dim_avg_n_Wrap(r_TempHist, 0)
            r_TempHist          := conform(r_Temp, r_TempHist, (/1, 2/))
            copy_VarCoords(r_Temp, r_TempHist)
            r_Temp               = r_Temp - r_TempHist
            ;
            if (isdefined("r_TempAllReal") .eq. False .or. iReal .eq. 0) then
                i_dims          := array_append_record(dimsizes(s_Real), \
                                                    dimsizes(r_Temp(:2,:,:)), 0)
                r_TempAllReal   := new(i_dims, float)
                copy_VarCoords(r_Temp(0, :, :), r_TempAllReal(0, 0, :, :))
                r_TempAllReal!0  = "Real"
                r_TempAllReal&Real  = s_Real
                r_TempAllReal!1  = "Terms"
                r_TempAllReal&Terms = (/"near", "mid", "long"/)
            end if
            ;
            do iTerm = 0, dimsizes(i_yrProj(:, 0))-1
                i_TermYr        := i_yrProj(iTerm, :)
                r_TempAllReal(iReal, iTerm, :, :) = dim_avg_n_Wrap(\
                                    r_Temp({i_TermYr(0):i_TermYr(1)}, :, :), 0)
            end do
        end do 
    ;
        if (isdefined("r_TempAllReal") ) then 
            if (isdefined("r_VarRaw") .eq. False) then
                i_dims          := array_append_record(dimsizes(s_ModelSelected), \
                                        dimsizes(r_TempAllReal(0, :2, :, :)), 0)
                r_VarRaw        := new(i_dims, float)
                copy_VarCoords(r_TempAllReal(0, :, :, :), r_VarRaw(0, :, :, :))
        ;
                r_VarRaw!0       = "Models"
                r_VarRaw&Models  = s_ModelSelected
        ;
            end if   
            r_VarRaw(iModel, :, :, :)   = dim_avg_n_Wrap(r_TempAllReal, 0)
            delete([/r_TempAllReal/])
        end if 
    ;
    end do
;
    r_VarRawMME        := dim_avg_n_Wrap(r_VarRaw, 0)
    ;
r_lonForDif            := (/min(r_Warming&lon), max(r_Warming&lon)/)
r_latForDif            := (/min(r_Warming&lat), max(r_Warming&lat)/)
r_VarRawMMEForDif      := r_VarRawMME(:, {r_latForDif(0)-1:r_latForDif(1)+1}, \
                                   {r_lonForDif(0)-1:r_lonForDif(1)+1})
;
;; Area Mean of projection warming 
r_VarRawChina           = my_SelChina([/r_VarRaw/])
r_VarRawChina          := my_AreaMean2([/r_VarRawChina/], r_VarRawChina&lat, \
                                       r_VarRawChina&lon, 0)
r_VarRawMMEChina        = my_SelChina([/r_VarRawMME/])
r_VarRawMMEChina       := my_AreaMean2([/r_VarRawMMEChina/], \
                                r_VarRawMMEChina&lat, r_VarRawMMEChina&lon, 0)
; printVarSummary(r_VarRawChina) ; [Models | 22] x [Terms | 3]
;; constrained results 
;;; r_China, r_Reg_Warming_All
;
; Dif
r_WarmingDif            = r_Warming
r_WarmingDif            = r_WarmingDif - r_VarRawMME
;; Area-mean of difference over China 
r_ChinaDif             := my_SelChina([/r_WarmingDif/])
r_ChinaDif             := my_AreaMean2([/r_ChinaDif/], r_ChinaDif&lat, r_ChinaDif&lon, 0)
print("3 Term Warming Dif:")
print(r_ChinaDif + " K")
print("")
;
print(s_Exp + "")
do iTerm = 0, 2
    ; By quater number 
    i_Ind              := ind(ismissing(r_VarRawChina(:, iTerm)) .eq. False)
    r_Temp             := r_VarRawChina(i_Ind, iTerm)
    r_Temp             := my_QuaterNum_withMean3(r_Temp, 95, 100)
    ;
    ; ;; By STD 
    ; r_TempSTD          := stddev(r_VarRawChina(:, iTerm))
    ; r_TempSTD          := (/r_VarRawMMEChina(iTerm) - r_TempSTD, \
    ;                         r_VarRawMMEChina(iTerm) + r_TempSTD/)
    ;
    ; ;; Raw SNR 
    ; r_SNR_Raw          := abs(r_VarRawMMEChina(iTerm)) / (r_Temp(1) - r_Temp(0))
    ;
    ;; Constrained 
    ;; By quater number 
    r_Temp_Con         := r_Reg_Warming_All(i_Ind, iTerm)
    r_Temp_Con         := my_QuaterNum_withMean3(r_Temp_Con, 95, 100)
    ; ;; By STD 
    ; r_Temp_ConSTD      := stddev(r_Reg_Warming_All(:, iTerm))
    ; r_Temp_ConSTD      := (/r_China(iTerm) - r_Temp_ConSTD, \
    ;                         r_China(iTerm) + r_Temp_ConSTD/)
    ;
    if (iTerm .eq. 0) then
        r_Range        := new((/3, 3/), float)
        r_Range!0       = "Terms"
        r_Range&Terms   = ispan(0, 2, 1)
        r_Range!1       = "MME_95th_100th"
        r_Range&$r_Range!1$ = ispan(0, 2, 1)
        ;
        r_RangeRaw      = r_Range
        ;
        r_RangeWithSTD    = r_Range
        r_RangeWithSTD!1  = "MME_Plus_Minus_STD"
        r_RangeRawWithSTD = r_RangeWithSTD
    end if
    ; r_Range(iTerm, 0)   = r_China(iTerm)
    r_Range(iTerm, :)   = r_Temp_Con + 0
    ;
    ; r_RangeRaw(iTerm, 0)= r_VarRawMMEChina(iTerm)
    r_RangeRaw(iTerm, :)= r_Temp + 0
    ;
    ; r_RangeWithSTD(iTerm, 0) = r_China(iTerm) + 0
    ; r_RangeWithSTD(iTerm, 1:)= r_Temp_ConSTD + 0
    ; ;
    ; r_RangeRawWithSTD(iTerm, 0) = r_VarRawMMEChina(iTerm) + 0
    ; r_RangeRawWithSTD(iTerm, 1:)= r_TempSTD + 0
    ;
    ; ;; SNR 
    ; r_SNR              := abs(r_China(iTerm)) / (r_Temp_Con(1) - r_Temp_Con(0))
    ; ;
    ; ;; Range dif 
    ; r_RangeDif         := r_Temp_Con(1) - r_Temp_Con(0)
    ; r_RangeDif          = r_RangeDif - (r_Temp(1) - r_Temp(0))
    ; r_RangeDif          = r_RangeDif / (r_Temp(1) - r_Temp(0)) * 100.
    ; ;
    ; ; Dif 
    ; r_Dif_Range        := r_Reg_Warming_All(i_Ind, iTerm) - \
    ;                         r_VarRawChina(i_Ind, iTerm)
    ; r_Dif_Range        := my_QuaterNum_withMean3(r_Dif_Range, 5, 95)
    ; ;
    ; ;; By STD 
    ; r_Dif_Range        := r_Reg_Warming_All(:, iTerm) - r_VarRawChina(:, iTerm)
    ; r_Dif_Range        := stddev(r_Dif_Range)
    ; r_Dif_Range        := (/r_ChinaDif(iTerm) - r_Dif_Range, \
    ;                         r_ChinaDif(iTerm) + r_Dif_Range/)
    print(s_yrProj(iTerm) + "")
    print("Raw Proj & SNR:")
    print(sprintf("%4.2f", r_Temp(0)) + " (" + sprintf("%4.2f", r_Temp(1)) + \
          "~" + sprintf("%4.2f", r_Temp(2)) + "), ")
    print("Con Proj & SNR:")
    print(sprintf("%4.2f", r_Temp_Con(0)) + " (" + \
          sprintf("%4.2f", r_Temp_Con(1)) + "~" + \
          sprintf("%4.2f", r_Temp_Con(2)) + ")")
    ; print("Dif:")
    ; print(sprintf("%4.2f", r_ChinaDif(iTerm)) + " (" + \
    ;       sprintf("%4.2f", r_Dif_Range(0)) + "~" + \
    ;       sprintf("%4.2f", r_Dif_Range(1)) + ")")
    print("")
    ;
end do
; s_VarName2              = "ConstrainedWarmingChina_" + s_Exp
; FileOut->$s_VarName2$   = r_Range
; s_VarName2              = "RawWarmingChina_" + s_Exp
; FileOut->$s_VarName2$   = r_RangeRaw
; ;
; s_VarName2              = "ConstrainedWarmingChina_WithSTD_" + s_Exp
; FileOut->$s_VarName2$   = r_RangeWithSTD
; s_VarName2              = "RawWarmingChina_WithSTD_" + s_Exp
; FileOut->$s_VarName2$   = r_RangeRawWithSTD
;
;;
end













