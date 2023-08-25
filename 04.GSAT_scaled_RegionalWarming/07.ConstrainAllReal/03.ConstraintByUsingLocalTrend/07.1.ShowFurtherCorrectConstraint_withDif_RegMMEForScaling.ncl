; ==================================================
; 本程序基于GSAT与局地升温的回归系数，对中国地区的升温进行约束
;   1. 读取回归系数，然后读取约束后的GSAT，两者相乘和加上截距
;   2. 读取约束前的结果，然后约束前后作差
; Attention: This code only show the projected warming in the constrained projection,
; and the raw projection, and the 5~95th as well
; ==================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code/my_FunctionInPNNL.ncl"
;
undef("my_Add_SouthChinaSea")
function my_Add_SouthChinaSea(map, Draw:numeric, DrawSigLev:numeric, \
                              SigLev:numeric, wks)
; 本程序单独绘制指定图形右下角的南海区域
;; map表示已经绘制完成的plot变量
;; Draw表示被绘制的矩阵
;; DrawSigLev表示显著性水平的矩阵，其维度必须与Draw一致，如果没有，全部设为0即可
;; SigLev表示打点的显著性水平，缺测值设为-999.
;; wks表示画布
local r_Draw
begin
;
r_Draw                              = Draw
r_DrawSigLev                        = DrawSigLev
r_DrawSigLev@_FillValue             = -999.
r_SigLev                            = SigLev
r_SigLev@_FillValue                 = -999. 

r_Draw!0                            = "lat"
r_Draw!1                            = "lon"
r_DrawSigLev!0                      = "lat"
r_DrawSigLev!1                      = "lon"

; --------------- add South China Sea --------------
    nhres                           = True
    nhres@gsnFrame                  = False
    nhres@gsnDraw                   = False

    nhres@vpHeightF                 = 0.13
    nhres@vpWidthF                  = 0.13

    nhres@mpMinLatF                 = 2.
    nhres@mpMaxLatF                 = 23.
    nhres@mpMinLonF                 = 105.
    nhres@mpMaxLonF                 = 123.

    getvalues map
    "mpDataSetName"             : nhres@mpDataSetName
    "mpDataBaseVersion"         : nhres@mpDataBaseVersion
    "mpFillOn"                  : nhres@mpFillOn
    "mpFillDrawOrder"           : nhres@mpFillDrawOrder
    "mpAreaMaskingOn"           : nhres@mpAreaMaskingOn
    "mpMaskAreaSpecifiers"      : nhres@mpMaskAreaSpecifiers ;;;;;; very important
    "mpOutlineSpecifiers"       : nhres@mpOutlineSpecifiers
    "mpOutlineBoundarySets"     : nhres@mpOutlineBoundarySets
    "mpLandFillColor"           : nhres@mpLandFillColor
    "mpOceanFillColor"          : nhres@mpOceanFillColor
    "mpInlandWaterFillColor"    : nhres@mpInlandWaterFillColor
    "mpNationalLineColor"       : nhres@mpNationalLineColor
    "mpProvincialLineColor"     : nhres@mpProvincialLineColor
    end getvalues

    getvalues map@contour
    "cnFillOn"                  : nhres@cnFillOn
    "cnLevelSelectionMode"      : nhres@cnLevelSelectionMode
    "cnLevels"                  : nhres@cnLevels
    "cnFillColors"              : nhres@cnFillColors
    "cnFillDrawOrder"           : nhres@cnFillDrawOrder
    "cnLinesOn"                 : nhres@cnLinesOn
    "cnLineLabelsOn"            : nhres@cnLineLabelsOn
    end getvalues

    nhres@lbLabelBarOn              = False
    nhres@tmXBOn                    = False
    nhres@tmYLOn                    = False
    nhres@tmYROn                    = False
    nhres@tmXTOn                    = False
    nhres@gsnLeftString             = ""
    nhres@gsnRightString            = ""

    if ( (min(r_Draw&lon) .gt. 1 .and. max(r_Draw&lon) .lt. 359) .or.\
         (min(r_Draw&lon)*max(r_Draw&lon) .gt. 0) ) then
        nhres@gsnAddCyclic          = False
    end if

    map_nanhai = gsn_csm_contour_map(wks, r_Draw, nhres)
    ;
    if (all(ismissing(r_DrawSigLev)) .eq. False .and. \
        ismissing(r_SigLev) .eq. False) then 
        sres    = True
        sres@gsnDraw                = False
        sres@gsnFrame               = False
        sres@cnLinesOn              = False
        sres@cnFillOn               = False
        sres@cnLineLabelsOn         = False
        sres@cnInfoLabelOn          = False
        sres@gsnLeftString          = ""
        sres@gsnCenterString        = ""
        sres@gsnRightString         = ""
        sres@cnFillDrawOrder        = "PreDraw"
        opt1    = True
        opt1@gsnShadeFillType       = "pattern"
        ; opt1@gsnShadeHigh              = 3
        opt1@gsnShadeLow            = 17
        sres@gsnAddCyclic           = True
        opt1@gsnAddCyclic           = True
        if ( min(r_DrawSigLev&lon) .gt. 1 .and. max(r_DrawSigLev&lon) .lt. 359 \
            .or. (min(r_Draw&lon)*max(r_Draw&lon) .gt. 0) ) then
            sres@gsnAddCyclic       = False
            opt1@gsnAddCyclic       = False
        end if

        map_stipple                 = gsn_csm_contour(wks, r_DrawSigLev,\
                                                     sres)
        plotSig                     = gsn_contour_shade(map_stipple, r_SigLev, \
                                                       1.6, opt1)
        overlay(map_nanhai, plotSig)
    else 
        plotSig                     = False
    end if 
    ;
    ;;; 另一种overlay
    adres                           = True
    adres@amParallelPosF            = 0.495
    adres@amOrthogonalPosF          = 0.49
    adres@amJust                    = "BottomRight"

    plotnh = gsn_add_annotation(map, map_nanhai, adres)
    ; plotSig = gsn_add_annotation(map, plotSig, adres)
    return([/plotnh, plotSig/])
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
; ==================================================
;   1. 读取回归系数，然后读取约束后的GSAT，两者相乘和加上截距
s_4Models                 = "_23Models" ; "_5Models" ; "_30Models" ; "" ; 
; s_ModelSelected           = (/"ACCESS-CM2", "CNRM-CM6-1", "EC-Earth3", \
;                             "GFDL-CM4", "MPI-ESM1-2-LR"/)
s_ModelSelected           = (/ "ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                         "BCC-CSM2-MR", "CAMS-CSM1-0", "CanESM5", \
                         "CanESM5-CanOE", "CNRM-CM6-1", "CNRM-ESM2-1", \
                         "EC-Earth3", "EC-Earth3-Veg", "FGOALS-f3-L", \
                         "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
                         "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8", \
                         "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR", "MRI-ESM2-0", \
                         "UKESM1-0-LL"/)
; s_ModelSel= (/ "ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR", \
;                 "CAMS-CSM1-0", "CanESM5", "CanESM5-CanOE", "CESM2-WACCM", \
;                 "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3", "EC-Earth3-Veg", \
;                 "FGOALS-f3-L", "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
;                 "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8", "INM-CM5-0", \
;                 "IPSL-CM6A-LR", "MCM-UA-1-0", "MIROC6", "MIROC-ES2L", \
;                 "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NESM3", \
;                 "NorESM2-LM", "UKESM1-0-LL"/)

s_Month                   = "AnMean"
s_Exp                     = "ssp245" ; "ssp126" ; "ssp370" ; "ssp585" ; 

i_yrProj                  = (/ (/2021, 2040/), (/2041, 2060/), \
                                   (/2080, 2099/) /)
s_yrProj                  = (/"near", "mid", "long"/)
;
s_DirOut = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/03.ConstraintByUsingLocalTrend/"
s_FileOut                 = s_DirOut + "ConstrainedWarmingInChina" + s_4Models \
                            + ".nc"
ls                        = systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileOut               = addfile(s_FileOut, "c")
    GAttr                 = True
    GAttr@history         = systemfunc("date") + "  " + systemfunc("pwd") + \
        "/07.1.07.1.ShowFurtherCorrectConstraint_withDif_RegMMEForScaling.ncl"
    fileattdef(FileOut, GAttr)
else 
    FileOut               = addfile(s_FileOut, "w")
end if
;
;;;;; for Pattern 
s_VarName                 = "tas" ; "pr" ; 
s_DirRead = "$WACCEM/my_data/CMIP6/tas/POST/ProjectionOfAllReal/"
s_FileRead                = s_DirRead + "ProjectionInDifTerms_AllReal" + \
                            s_4Models + ".nc"
FileRead_Raw              = addfile(s_FileRead, "r")
;;;;; for Pattern 

s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/03.ConstraintByUsingLocalTrend/"
s_FileRead                = s_DirRead + "FurtherCorrectConstrainedWarming_" + \
                            s_Exp + s_4Models + ".nc"
FileConstrained           = addfile(s_FileRead, "r")

s_VarName2              = "ConstrainedWarming_" + s_Exp
r_Warming               = FileConstrained->$s_VarName2$
;
;;; Area-Mean over China 
r_China                := my_SelChina([/r_Warming/])
;
r_China                := my_AreaMean2([/r_China/], r_China&lat, r_China&lon, 0)
print("3 Term Warming:")
print(r_China + " K")
;
do iModel = 0, dimsizes(s_ModelSelected)-1
    s_VarName2          = "ConstrainedWarming_" + s_ModelSelected(iModel) + "_"\
                             + s_Exp
    L_Return            = my_VarNameExist(FileConstrained, s_VarName2)
    r_Reg_Warming      := FileConstrained->$s_VarName2$
    ;
    r_Reg_Warming      := my_SelChina([/r_Reg_Warming/])
    ;
    r_Reg_Warming      := my_AreaMean2([/r_Reg_Warming/], r_Reg_Warming&lat, \
                                       r_Reg_Warming&lon, 0)
    if (isdefined("r_Reg_Warming_All") .eq. False) then
        r_Reg_Warming_All        = new((/dimsizes(s_ModelSelected), 3/), float)
        copy_VarCoords(r_Reg_Warming, r_Reg_Warming_All(0, :))
        r_Reg_Warming_All!0      = "Models"
        r_Reg_Warming_All&Models = s_ModelSelected
    end if
    r_Reg_Warming_All(iModel, :) = r_Reg_Warming
end do
print(r_Reg_Warming_All&Models + ": " + r_Reg_Warming_All(:, 2))
r_Temp                 := my_QuaterNum2(r_Reg_Warming_All(:, 2), 5, 95)
print("")
print("MME: " + avg(r_Reg_Warming_All(:, 2)))
print("5th to 95th:")
print(r_Temp + "")
exit()
r_Reg_Warming_Spread             = dim_stddev_n_Wrap(r_Reg_Warming_All, 0)
print("STD: " + r_Reg_Warming_Spread(2))
; exit()
;
; r_China_MME_Spread               = new((/3, 2/), float)
; copy_VarMeta(r_China, r_China_MME_Spread(:, 0))
; r_China_MME_Spread!1             = "MME_STD"
; r_China_MME_Spread&MME_STD       = (/0, 1/)
; do iPer = 0, 2
;     r_China_MME_Spread(iPer, 0)  = r_China(iPer)
;     r_China_MME_Spread(iPer, 1)  = r_Reg_Warming_Spread(iPer)
; end do
; ; printVarSummary(r_China_MME_Spread) ; [ProjectionTerm | 3] x [MME_STD | 2]
;;
s_VarNameAll                     = getfilevarnames(FileRead_Raw)
    do iModel = 0, dimsizes(s_ModelSelected)-1
        ;
        l_list                  := my_CMIP6Model_Exp_Path_ForAllReal(\
                                        s_ModelSelected(iModel), s_Exp, "Amon")
        s_Real                  := l_list[1]
        l_list_                 := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                        s_ModelSelected(iModel), s_Exp, "Amon")
        s_Real_                 := l_list_[1]
        if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real)) then
            s_Real              := s_Real_
        end if
        ;
        do iReal = 0, dimsizes(s_Real)-1
            s_VarName2          := s_VarName + "_inDifTerms_" + \
                                s_ModelSelected(iModel) + "_" + s_Exp + "_" + \
                                s_Real(iReal)
            i_Ind               := ind(s_VarName2 .eq. s_VarNameAll)
            if (ismissing(i_Ind(0))) then
                print(s_VarName2 + " Missing!")
                continue 
            end if
    ;
            r_Temp              := FileRead_Raw->$s_VarName2$(:, {-90:90}, \
                                                              {0:360})
            if (isdefined("r_TempAllReal") .eq. False .or. iReal .eq. 0) then
                i_dims          := array_append_record(dimsizes(s_Real), \
                                                       dimsizes(r_Temp), 0)
                r_TempAllReal   := new(i_dims, float)
                copy_VarCoords(r_Temp, r_TempAllReal(0, :, :, :))
                r_TempAllReal!0  = "Real"
                r_TempAllReal&Real= s_Real
            end if
            r_TempAllReal(iReal, :, :, :) = r_Temp
        end do 
    ;
        if (isdefined("r_TempAllReal") ) then 
            if (isdefined("r_VarRaw") .eq. False) then
                i_dims          := array_append_record(\
                            dimsizes(s_ModelSelected), dimsizes(r_Temp), 0)
                r_VarRaw        := new(i_dims, float)
                copy_VarCoords(r_Temp, r_VarRaw(0, :, :, :))
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
    ;; By quater number 
    ; i_Ind              := ind(ismissing(r_VarRawChina(:, iTerm)) .eq. False)
    ; r_Temp             := r_VarRawChina(i_Ind, iTerm)
    ; r_Temp             := my_QuaterNum2(r_Temp, 5, 95)
    ;
    ;; By STD 
    r_Temp             := stddev(r_VarRawChina(:, iTerm))
    r_Temp             := (/r_VarRawMMEChina(iTerm) - r_Temp, \
                            r_VarRawMMEChina(iTerm) + r_Temp/)
    ;
    ;; Raw SNR 
    r_SNR_Raw          := abs(r_VarRawMMEChina(iTerm)) / (r_Temp(1) - r_Temp(0))
    ;
    ;; Constrained 
    ; ;; By quater number 
    ; r_Temp_Con         := r_Reg_Warming_All(i_Ind, iTerm)
    ; r_Temp_Con         := my_QuaterNum2(r_Temp_Con, 5, 95)
    ;; By STD 
    r_Temp_Con         := stddev(r_Reg_Warming_All(:, iTerm))
    r_Temp_Con         := (/r_China(iTerm) - r_Temp_Con, r_China(iTerm) + r_Temp_Con/)
    ;
    if (iTerm .eq. 0) then
        r_Range        := new((/3, 3/), float)
        r_Range!0       = "Terms"
        r_Range&Terms   = ispan(0, 2, 1)
        r_Range!1       = "MME_5th_95th"
        r_Range&$r_Range!1$ = ispan(0, 2, 1)
        ;
        r_RangeRaw      = r_Range
    end if
    r_Range(iTerm, 0)   = r_China(iTerm)
    r_Range(iTerm, 1:)  = r_Temp_Con
    ;
    r_RangeRaw(iTerm, 0)= r_VarRawMMEChina(iTerm)
    r_RangeRaw(iTerm,1:)= r_Temp
    ;
    ;; SNR 
    r_SNR              := abs(r_China(iTerm)) / (r_Temp_Con(1) - r_Temp_Con(0))
    ;
    ;; Dif 
    ; r_Dif_Range        := r_Reg_Warming_All(i_Ind, iTerm) - \
    ;                         r_VarRawChina(i_Ind, iTerm)
    ; r_Dif_Range        := my_QuaterNum2(r_Dif_Range, 5, 95)
    ;
    ;; By STD 
    r_Dif_Range        := r_Reg_Warming_All(:, iTerm) - r_VarRawChina(:, iTerm)
    r_Dif_Range        := stddev(r_Dif_Range)
    r_Dif_Range        := (/r_ChinaDif(iTerm) - r_Dif_Range, \
                            r_ChinaDif(iTerm) + r_Dif_Range/)
    print(s_yrProj(iTerm) + "")
    print("Raw Proj & SNR:")
    print(sprintf("%4.2f", r_VarRawMMEChina(iTerm)) + " (" + \
          sprintf("%4.2f", r_Temp(0)) + "~" + sprintf("%4.2f", r_Temp(1)) + "), "\
          + sprintf("%4.2f", r_SNR_Raw))
    print("Con Proj & SNR:")
    print(sprintf("%4.2f", r_China(iTerm)) + " (" + \
          sprintf("%4.2f", r_Temp_Con(0)) + "~" + \
          sprintf("%4.2f", r_Temp_Con(1)) + "), " + sprintf("%4.2f", r_SNR))
    print("Dif:")
    print(sprintf("%4.2f", r_ChinaDif(iTerm)) + " (" + \
          sprintf("%4.2f", r_Dif_Range(0)) + "~" + \
          sprintf("%4.2f", r_Dif_Range(1)) + ")")
    print("")
    ;
end do
s_VarName2              = "ConstrainedWarmingChina_" + s_Exp
FileOut->$s_VarName2$   = r_Range
s_VarName2              = "RawWarmingChina_" + s_Exp
FileOut->$s_VarName2$   = r_RangeRaw
;
;;
end















