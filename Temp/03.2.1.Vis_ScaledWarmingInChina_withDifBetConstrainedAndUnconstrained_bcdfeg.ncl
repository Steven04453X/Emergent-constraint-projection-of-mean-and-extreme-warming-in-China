; ==================================================
; 本程序基于GSAT与局地升温的回归系数，对中国地区的升温进行约束
;   1. 读取回归系数，然后读取约束后的GSAT，两者相乘和加上截距
;   2. 读取约束前的结果，然后约束前后作差
;   3. 进行可视化
; 注意本程序在旁边附上约束前后的差异
; ==================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/WORK/chenzm/Work/Others/Special_Code/my_Function.ncl"
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
    s_DirRead         = "/WORK/pengdd/Mydata/Obs/ChinaSTN/CN05.1/"
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
    r_lonTaiwan       = (/119, 125/)
    r_latTaiwan       = (/20, 26/)
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
s_4Models                       = "_30Models" ; "" ; 
; s_ModelSelected                 = (/"BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1",\
;                                    "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
;                                    "EC-Earth3-Veg", "FGOALS-f3-L", "FGOALS-g3",\
;                                    "GFDL-CM4", "GFDL-ESM4", "INM-CM5-0", \
;                                    "IPSL-CM6A-LR", "MCM-UA-1-0", "MIROC6", \
;                                    "MIROC-ES2L", "MRI-ESM2-0", "NESM3", \
;                                    "UKESM1-0-LL"/)
s_ModelSelected = (/ "ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR", \
                "CAMS-CSM1-0", "CanESM5", "CanESM5-CanOE", "CESM2-WACCM", \
                "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3", "EC-Earth3-Veg", \
                "FGOALS-f3-L", "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
                "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8", "INM-CM5-0", \
                "IPSL-CM6A-LR", "MCM-UA-1-0", "MIROC6", "MIROC-ES2L", \
                "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NESM3", \
                "NorESM2-LM", "UKESM1-0-LL"/)

s_Month                         = "AnMean"
s_Exp                           = "ssp245" ; "ssp370" ; "ssp126" ; "ssp585" ; 
s_ExpTitle                      = "SSP2-4.5" ; "SSP3-7.0" ; "SSP1-2.6" ; "SSP5-8.5" ; 
s_DirOut  = "/WORK/chenzm/Work/OtherWork/2020/ProjectionInChina/Fig/04.GSAT_scaled_RegionalWarming/"
system("mkdir -p " + s_DirOut)
s_FigOut                        = s_DirOut + \
    "03.2.1.Vis_ScaledWarmingInChina_withDifBetConstrainedAndUnconstrained_bcdefg" + \
    s_Exp + s_4Models
wks                             = gsn_open_wks("pdf", s_FigOut)

i_yrProj                  = (/ (/2021, 2040/), (/2041, 2060/), \
                                   (/2080, 2099/) /)

;;;;; for Pattern 
s_VarName                 = "tas" ; "pr" ; 
s_DirRead = "/WORK/chenzm/Work/OtherWork/2020/ProjectionInChina/Data/01.ChinaMean/"
s_FileRead                = s_DirRead + s_VarName + "_ProjectionInChinaPattern_" + \
                        s_Month + ".nc"
FileRead_Raw              = addfile(s_FileRead, "r")
;;;;; for Pattern 

s_DirRead   = "/data1_99/CMIP6/chenzm/my_data/CMIP6/tas/POST/GSAT_scaled_RegionalWarming/"
s_FileRead              = s_DirRead + \
        "GSAT_reg_RegionalWarming_AnMean_" + s_Exp + s_4Models + ".nc"
FileRead                = addfile(s_FileRead, "r")

r_Reg                   = FileRead->Reg(0, :, :, :)
r_yintercept            = FileRead->yintercept

;; Constrained GSAT 
s_DirRead   = "/WORK/chenzm/Work/OtherWork/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/"
s_FileRead              = s_DirRead + "ObsX_ObsY" + s_4Models + "_ObsCorr.nc"
; s_FileRead              = s_DirRead + "ObsX_ObsY" + s_4Models + ".nc"
FileConstrainedGSAT     = addfile(s_FileRead, "r")
s_VarName2              = "ObsXY_" + s_Exp
r_GSAT                  = FileConstrainedGSAT->$s_VarName2$(:, 1)
print("Constrained GSAT: " + r_GSAT + "")
; exit()
; r_GSAT                  = (/0.72, 1.41, 3.21/)  ; units: K

r_Warming               = r_Reg
do iW = 0, 2
    r_Warming(iW, :, :) = r_Reg(iW, :, :) * r_GSAT(iW) + r_yintercept(iW, :, :)
end do
;
;;; Area-Mean over China 
r_China                := my_SelChina([/r_Warming/])
r_China                := my_AreaMean2([/r_China/], r_China&lat, r_China&lon, 0)
print("3 Term Warming:")
print(r_China + " K")
;
;; 单独计算中国地区的区域平均，以及spread
do iModel = 0, dimsizes(s_ModelSelected)-1
    s_VarName2          = "ConstrainedGSAT_" + s_ModelSelected(iModel) + "_" + \
                            s_Exp
    L_Return            = my_VarNameExist(FileConstrainedGSAT, s_VarName2)
    r_Reg_Warming      := FileConstrainedGSAT->$s_VarName2$
    ; print(s_ModelSelected(iModel) + ": " + r_Reg_Warming + "")
    ;
    r_Reg_Warming      := conform(r_Reg, r_Reg_Warming, (/0/))
    copy_VarCoords(r_Reg, r_Reg_Warming)
    r_Reg_Warming       = r_Reg * r_Reg_Warming + r_yintercept(:, :, :)
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
; print(r_Reg_Warming_All&Models + ": " + r_Reg_Warming_All(:, 2))
r_Reg_Warming_Spread             = dim_stddev_n_Wrap(r_Reg_Warming_All, 0)
; print("STD: " + r_Reg_Warming_Spread(2))
; exit()

r_China_MME_Spread               = new((/3, 2/), float)
copy_VarMeta(r_China, r_China_MME_Spread(:, 0))
r_China_MME_Spread!1             = "MME_STD"
r_China_MME_Spread&MME_STD       = (/0, 1/)
do iPer = 0, 2
    r_China_MME_Spread(iPer, 0)  = r_China(iPer)
    r_China_MME_Spread(iPer, 1)  = r_Reg_Warming_Spread(iPer)
end do
; printVarSummary(r_China_MME_Spread) ; [ProjectionTerm | 3] x [MME_STD | 2]
;
;; Output 
s_DirOut = "/WORK/chenzm/Work/OtherWork/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/"
s_FileOut                        = s_DirOut + "ConstrainedWarmingInChina" + \
                                s_4Models + ".nc"
ls                               = systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileOut                      = addfile(s_FileOut, "c")
    GAttr                        = True
    GAttr@history                = systemfunc("date") + ";  " + \
                        systemfunc("pwd") + "/03.2.Vis_ScaledWarmingInChina.ncl"
    GAttr@Modles                 = s_ModelSelected
    fileattdef(FileOut, GAttr)
else 
    FileOut                      = addfile(s_FileOut, "w")
    FileOut@history              = systemfunc("date") + ";  " + \
                        systemfunc("pwd") + "/03.2.Vis_ScaledWarmingInChina.ncl"
end if
s_VarName2                       = "ConstrainedWarmingChina_" + s_Exp
FileOut->$s_VarName2$            = r_China_MME_Spread
; exit()
;
;; difference between constrained and unconstrained projection 
;
s_VarNameAll                     = getfilevarnames(FileRead_Raw)
    do iModel = 0, dimsizes(s_ModelSelected)-1
        s_VarName2              := s_VarName + "_" + s_ModelSelected(iModel) + \
                                "_" + s_Exp
        i_Ind                   := ind(s_VarName2 .eq. s_VarNameAll)
        if (ismissing(i_Ind(0))) then
            print(s_VarName2 + " Missing!")
            continue 
        end if
    ;
        r_Temp                  := FileRead_Raw->$s_VarName2$(:, {-90:90}, \
                                                              {0:360})
    ;
        ;; 计算降水的变化率
        if (s_VarName .eq. "pr") then
            s_VarName2          := s_VarName + "_" + s_ModelSelected(iModel) + \
                                    "_historical"
            r_Clim              := FileRead_Raw->$s_VarName2$({-90:90}, {0:360})
            r_Clim              := conform(r_Temp, r_Clim, (/1, 2/))
            copy_VarCoords(r_Temp, r_Clim)
            r_Clim               = where(r_Clim .eq. 0, r_Clim@_FillValue, \
                                         r_Clim)
            r_Temp               = r_Temp / r_Clim * 100.
        end if
    ;
        if (isdefined("r_VarRaw") .eq. False) then
            i_dims              := array_append_record(\
                                dimsizes(s_ModelSelected), dimsizes(r_Temp), 0)
            r_VarRaw            := new(i_dims, float)
            copy_VarCoords(r_Temp, r_VarRaw(0, :, :, :))
        ;
            r_VarRaw!0           = "Models"
            r_VarRaw&Models      = s_ModelSelected
        ;
        end if
    ;   
        r_VarRaw(iModel, :, :, :)   = r_Temp
    ;
    end do
;
    r_VarRawMME        := dim_avg_n_Wrap(r_VarRaw, 0)
;
r_WarmingDif            = r_Warming
r_WarmingDif            = r_WarmingDif - r_VarRawMME
;
;; Area-mean of difference over China 
r_China                := my_SelChina([/r_WarmingDif/])
r_China                := my_AreaMean2([/r_China/], r_China&lat, r_China&lon, 0)
print("3 Term Warming:")
print(r_China + " K")
; exit()
;
; ==================================================
;   2. 进行可视化
res2                               = True
res2@cnInfoLabelOn                 = False
res2@cnFillOn                      = True
res2@cnLinesOn                     = False
res2@cnLineLabelsOn                = False
res2@gsnDraw                       = False
res2@gsnFrame                      = False
;
res2@gsnLeftStringFontHeightF      = 0.028
res2@gsnCenterStringFontHeightF    = 0.013
res2@gsnRightStringFontHeightF     = 0.024
res2@gsnLeftStringOrthogonalPosF   = 0.01
res2@gsnRightStringOrthogonalPosF  = 0.01
res2@gsnCenterStringOrthogonalPosF = 0.03

res2@pmTickMarkDisplayMode         = "Always"
res2@tmXBLabelFontHeightF          = 0.023
res2@tmYLLabelFontHeightF          = 0.023

res2@lbLabelBarOn                  = True
res2@lbLabelFontHeightF            = 0.02
res2@pmLabelBarWidthF              = 0.75
res2@pmLabelBarOrthogonalPosF      = 0.1

; --------------------------------------------------------
    ; set for the map 
    res2@mpMinLatF       = 15.
    res2@mpMaxLatF       = 55.
    res2@mpMinLonF       = 70.
    res2@mpMaxLonF       = 140.

    res2@mpFillOn        = True
    res2@mpDataSetName   = "/WORK/chenzm/my_data/topo/NCL-Chinamap-master/database/Earth..4"
    res2@mpDataBaseVersion           = "MediumRes" ; or "Ncarg4_1"
    res2@mpAreaMaskingOn             = True
    res2@mpMaskAreaSpecifiers        = (/ "China" /)
    res2@mpOutlineSpecifiers         = (/ "China", "China:Provinces" /)

    res2@mpLandFillColor             = "white"
    res2@mpInlandWaterFillColor      = "white"
    res2@mpOceanFillColor            = "white"
    res2@mpOutlineBoundarySets       = "NoBoundaries"
    res2@mpFillBoundarySets          = "NoBoundaries"
    res2@mpNationalLineColor         = "black"
    res2@mpProvincialLineColor       = "gray" ;"black"
    res2@mpGeophysicalLineColor      = "black"
    res2@mpNationalLineThicknessF    = 2
    res2@mpProvincialLineThicknessF  = 1
    res2@cnFillDrawOrder             = "PreDraw"

sres    = True
sres@gsnDraw                       = False
sres@gsnFrame                      = False
sres@cnLinesOn                     = False
sres@cnFillOn                      = False
sres@cnLineLabelsOn                = False
sres@cnInfoLabelOn                 = False
sres@gsnLeftString                 = ""
sres@gsnCenterString               = ""
sres@gsnRightString                = ""
sres@cnFillDrawOrder               = "PreDraw"
opt1        = True
opt1@gsnShadeFillType              = "pattern"
; opt1@gsnShadeHigh                  = 3
opt1@gsnShadeLow                   = 17
sres@gsnAddCyclic                  = True
opt1@gsnAddCyclic                  = True

; 2.25, 
r_cnLinesLevel_c1                 := (/1., 1.2, 1.6, 2, 2.5, 3, 3.5, 4, 4.5, 5/)
r_cnLinesLevel_c2                 := (/-1.2, -1., -0.8, -0.5, -0.4, -0.3, -0.2, -0.1/)
if (s_Exp .eq. "ssp126") then
    r_cnLinesLevel_c1             := (/0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.6/)
    r_cnLinesLevel_c2             := (/-0.9, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0/)
else if (s_Exp .eq. "ssp245") then
    ; r_cnLinesLevel_c1             := (/0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2/) + 0.5
    ; r_cnLinesLevel_c2             := (/-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0/)
    r_cnLinesLevel_c1             := (/0.7, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3./)
    r_cnLinesLevel_c2             := (/-0.9, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1/)
else if (s_Exp .eq. "ssp370") then
    r_cnLinesLevel_c1             := (/0.6, 0.8, 1.2, 1.4, 1.6, 1.8, 2.4, 2.8, 3.2, 3.6/) 
    r_cnLinesLevel_c2             := (/-1.5, -1, -0.7, -0.5, -0.4, -0.3, -0.2, -0.1/)
end if
end if
end if

s_LeftString                       = (/"(b)", "(c)", "(d)", "(e)", "(f)", "(g)"/)
plot1                              = new(6, graphic)
plotSig                            = plot1
plotnh                             = plot1
; plotnhSig                          = plotnh
r_Temp                            := r_Warming(0, :, :)
r_Temp                             = -999.
iSub2                              = 0
do iSub = 0, 2
    res2@gsnLeftString             = s_LeftString(iSub2) + " " + \
                                    i_yrProj(iSub, 0) + "-" + i_yrProj(iSub, 1)
    res2@gsnCenterString           = "" 
    res2@gsnRightString            = s_ExpTitle + ": " + s_Month
    ; ; -------------- Set the colorbar
    ; cmap   := (/ (/103, 0, 31/), (/178, 24, 43/), (/214, 96, 77/), \
    ;          (/244, 165, 130/), (/253, 219, 219/), (/240, 240, 240/),\ ; (/247, 247, 247/), \
    ;          (/209, 229, 240/), (/146, 197, 222/), (/67, 147, 195/), \
    ;          (/33, 102, 172/), (/5, 48, 97/) /) / 255. ; 

    ; cmap                          := cmap(::-1, :)
    ; res2@cnFillColors             := cmap
    ; res2@cnLevelSelectionMode      = "ExplicitLevels"
    
    i_dimColor                    := dimsizes(r_cnLinesLevel_c1)
    cmap                          := read_colormap_file("cmocean_matter")
    i_dimColor                     = dimsizes(cmap(:, 0)) / (i_dimColor+1)
    cmap                          := cmap(::i_dimColor, :)
    res2@cnFillColors             := cmap
    res2@cnLevelSelectionMode      = "ExplicitLevels"
    res2@cnLevels                 := r_cnLinesLevel_c1
    ;
    ; res2@cnLevelSelectionMode      = "ExplicitLevels"
    ; gsn_define_colormap(wks, "MPL_YlOrRd")
    ; i_cmap:= ispan(9, 127, 1)
    ; i_Temp:= (/0, 4/)
    ; i_cmap:= array_append_record(i_Temp, i_cmap, 0)
    ; r_cnLevels      := (/0.5, 1, 1.25, 1.5, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5/)
    ; ; (/0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5/)
    ; ; i_Temp          := toint(dimsizes(i_cmap)/dimsizes(r_cnLevels))-1
    ; res2@cnFillColors := i_cmap(::8)
    ; res2@cnLevels := r_cnLevels
    ;
    plot1(iSub2)                   = gsn_csm_contour_map(wks, \
                                                    r_Warming(iSub, :, :), res2)
    l_list                        := my_Add_SouthChinaSea(plot1(iSub2), \
                                                    r_Warming(iSub, :, :), \
                                                    r_Temp, \
                                                    -999., wks)
    plotnh(iSub2)                  = l_list[0]
    iSub2                          = iSub2 + 1
    ; plotnhSig(iSub)                = l_list[1]
    ;
    ;;;;;;;;;;;;;;;;;;;;;;;
    res2@gsnLeftString             = s_LeftString(iSub2) + \
                                    " Constrained - Unconstrained"
    res2@gsnRightString            = ""
    cmap                          := read_colormap_file("posneg_1")
    cmap                          := cmap(3:12, :)
    res2@cnLevels                 := r_cnLinesLevel_c2         
    res2@cnFillColors             := cmap
    ;
    plot1(iSub2)                   = gsn_csm_contour_map(wks, \
                                            r_WarmingDif(iSub, :, :), res2)
    l_list                        := my_Add_SouthChinaSea(plot1(iSub2), \
                                                    r_WarmingDif(iSub, :, :), \
                                                    r_Temp, \
                                                    -999., wks)
    plotnh(iSub2)                  = l_list[0]
    iSub2                          = iSub2 + 1
end do 

resP                               = True
resP@gsnMaximize                   = False
gsn_panel(wks, plot1, (/3, 2/), resP)
;
end