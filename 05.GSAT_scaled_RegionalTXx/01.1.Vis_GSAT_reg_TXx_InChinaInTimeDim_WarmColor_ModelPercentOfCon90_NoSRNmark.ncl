; ================================================================
; 本程序绘制不同预估时期内，中国地区的TXx与GSAT之间的回归系数
;   1. 指定数据，指定时期，然后读取回归系数和显著性水平
;   2. 绘图，分成三个图绘制
; 补充99%的信度检验，将colorbar改成暖色调
; 注意该脚本是在时间维上进行回归，对每个模式分别进行回归，然后计算MME
; ================================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code//my_FunctionInPNNL.ncl"
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
    s_DirRead         = "/global/cfs/cdirs/m1867/zmchen/my_data/Obs/CN05.1/"
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
        r_VarChina_({r_latTaiwan(0):r_latTiwan(1)}, \
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
undef("my_Add_SouthChinaSea")
function my_Add_SouthChinaSea(map, Draw:numeric, DrawSigLev:numeric, \
                              DrawSigLev_99:numeric, SigLev:numeric, wks)
; 本程序单独绘制指定图形右下角的南海区域
;; map表示已经绘制完成的plot变量
;; Draw表示被绘制的矩阵
;; DrawSigLev表示显著性水平的矩阵，其维度必须与Draw一致，如果没有，全部设为0即可
;; SigLev表示打点的显著性水平
;; wks表示画布
local r_Draw
begin
;
r_Draw                              = Draw
r_DrawSigLev                        = DrawSigLev
r_DrawSigLev_99                     = DrawSigLev_99
r_SigLev                            = SigLev
; r_SNR_Pattern                       = SNR_Pattern
; r_SNR_Thre                          = SNR_Thre

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
    sres    = True
    sres@gsnDraw                    = False
    sres@gsnFrame                   = False
    sres@cnLinesOn                  = False
    sres@cnFillOn                   = False
    sres@cnLineLabelsOn             = False
    sres@cnInfoLabelOn              = False
    sres@gsnLeftString              = ""
    sres@gsnCenterString            = ""
    sres@gsnRightString             = ""
    sres@cnFillDrawOrder            = "PreDraw"
    opt1        = True
    opt1@gsnShadeFillType           = "pattern"
    opt1@gsnShadeHigh                  = 17
    ; opt1@gsnShadeLow                = 17
    sres@gsnAddCyclic               = True
    opt1@gsnAddCyclic               = True

    plotSig                         = new(2, graphic)
    sres@cnFillColor               := "black"
    r_Temp                         := where(r_DrawSigLev .ge. r_SigLev, 1.1, 0)
    copy_VarCoords(r_DrawSigLev, r_Temp)
    r_Temp(0, 0)                    = 0
    map_stipple                     = gsn_csm_contour(wks, r_Temp, sres)
    plotSig(0)                      = gsn_contour_shade(map_stipple, 0, 1., opt1)
    overlay(map_nanhai, plotSig(0))
    ;
    sres@cnFillColor               := "white"
    r_Temp                         := where(r_DrawSigLev_99 .ge. r_SigLev, 1.1, 0)
    copy_VarCoords(r_DrawSigLev_99, r_Temp)
    r_Temp(0, 0)                    = 0
    map_stipple                     = gsn_csm_contour(wks, r_Temp, sres)
    plotSig(1)                      = gsn_contour_shade(map_stipple, 0, 1., opt1)
    overlay(map_nanhai, plotSig(1))
    ; if (dimsizes(r_SigLev) .eq. 1) then
    ;     r_Temp                     := where(r_DrawSigLev .ge. r_SigLev, 1.1, 0)
    ;     map_stipple                 = gsn_csm_contour(wks, r_Temp, sres)
    ;     plotSig                     = gsn_contour_shade(map_stipple, 0, 1., opt1)
    ;     overlay(map_nanhai, plotSig)
    ; else 
    ;     plotSig                     = new(dimsizes(r_SigLev), graphic)
    ;     r_ColorStippling            = (/ "black", "white" /)
    ;     do iStippling = 0, dimsizes(r_SigLev)-1
    ;         sres@cnFillColor       := r_ColorStippling(iStippling)
    ;         r_Temp                 := where(r_DrawSigLev .ge. r_SigLev(iStippling), \
    ;                                         1.1, 0)
    ;         copy_VarCoords(r_DrawSigLev, r_Temp)
    ;         r_Temp(0, 0)            = 0
    ;         map_stipple            := gsn_csm_contour(wks, r_Temp, sres)
    ;         plotSig(iStippling)     = gsn_contour_shade(map_stipple, 0, 1.,opt1)
    ;         overlay(map_nanhai, plotSig(iStippling))
    ;     end do
    ; end if
    
    ; delete([/opt1@gsnShadeHigh/])
    ; opt1@gsnShadeLow                = 3
    ; ; For SNR 
    ; r_SNR_Pattern(0, 0)             = 0
    ; if (dimsizes(r_SNR_Thre) .eq. 1) then
    ;     map_stipple                 = gsn_csm_contour(wks, r_SNR_Pattern, sres)
    ;     plotSNR                     = gsn_contour_shade(map_stipple, \
    ;                                                     r_SNR_Thre, 9999., opt1)
    ;     overlay(map_nanhai, plotSNR)
    ; else 
    ;     plotSNR                     = new(dimsizes(r_SNR_Thre), graphic)
    ;     r_ColorStippling            = (/ "gray", "black"/)
    ;     do iStippling = 0, dimsizes(r_SNR_Thre)-1
    ;         sres@cnFillColor       := r_ColorStippling(iStippling)
    ;         map_stipple            := gsn_csm_contour(wks, r_SNR_Pattern, sres)
    ;         plotSNR(iStippling)     = gsn_contour_shade(map_stipple, \
    ;                                         r_SNR_Thre(iStippling), 9999., opt1)
    ;         overlay(map_nanhai, plotSNR(iStippling))
    ;     end do
    ; end if    
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

begin
;
; ================================================================
;   1. 指定数据，指定时期，然后读取回归系数和显著性水平
i_yrProj                        = (/ (/2021, 2040/), (/2041, 2060/), \
                                   (/2080, 2099/) /)
s_Month                         = "AnMean"
s_Exp                           = "ssp585" ; "ssp245" ; "ssp126" ; "ssp370" ; 
s_ExpTitle                      = "SSP5-8.5" ; "SSP2-4.5" ; "SSP1-2.6" ; "SSP3-7.0" ; 
s_4Models                       = "_18Models"
s_ModelSelected    = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
if (s_Exp .eq. "ssp126") then
    s_ModelSelected      := (/"ACCESS-CM2", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
end if

i_yrRun                         = 20

r_SigLev_ModelPercentage        = 70 ; in percentage 

s_DirOut  = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Fig/05.GSAT_scaled_RegionalTXx/"
system("mkdir -p " + s_DirOut)
s_FigOut                        = s_DirOut + \
            "01.1.Vis_GSAT_reg_TXx_InChinaInTimeDim_WarmColor_ModelPercentOfCon90_NoSRNmark_" + s_Exp + \
        s_4Models
if (isdefined("i_yrRun")) then
    s_FigOut                    = s_DirOut + \
        "01.1.Vis_GSAT_reg_TXx_InChinaInTimeDim_WarmColor_ModelPercentOfCon90_NoSRNmark" + s_Exp + \
        "_" + i_yrRun + "yr" + s_4Models
end if
wks                             = gsn_open_wks("pdf", s_FigOut)

s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tasmax/POST/TXxRegGSAT/InTimeDim/"
s_FileRead                      = s_DirRead + "GSAT_reg_TXxWarmingInTimeDim_" +\
                                s_Month + "_" + s_Exp + "_" + i_yrRun + "yr" + \
                                s_4Models + ".nc"
FileRead                        = addfile(s_FileRead, "r")
;
s_VarName2                      = "Reg_MME_AcrossModel_" + s_Exp
r_Reg                           = FileRead->$s_VarName2$
;
r_RegChina                      = my_SelChina([/r_Reg(:, :, :)/])
r_RegChina                     := my_AreaMean2([/r_RegChina/], r_RegChina&lat, \
                                               r_RegChina&lon, 0)
print(r_RegChina&ProjectionTerm + ": " + r_RegChina)
; 
; ;; Confidence Level 
s_VarName2                      = "ModelPercentOfCon90_" + s_Exp
r_SigLevMME                     = FileRead->$s_VarName2$
s_VarName2                      = "ModelPercentOfCon99_" + s_Exp
r_SigLevMME_99                  = FileRead->$s_VarName2$
; s_VarName2                      = "RegInTimDim_MME_" + s_Exp
; r_SigLevMME                     = FileRead->$s_VarName2$(1, :, :, :)
;
; do iModel = 0, dimsizes(s_ModelSelected)-1
;     s_VarName2                  = "RegInTimDim_" + s_ModelSelected(iModel) + \
;                                 "_" + s_Exp
;     if (my_VarNameExist(FileRead, s_VarName2) .eq. False) then
;         print(s_VarName2 + " Missing!")
;         continue
;     end if
;     r_Temp                     := FileRead->$s_VarName2$(0, :, {-90:90}, {0:360})
;     ;
;     if (isdefined("r_RegAllModels") .eq. False .or. iModel .eq. 0) then
;         i_dims                 := array_append_record(dimsizes(s_ModelSelected), \
;                                                       dimsizes(r_Temp), 0)
;         r_RegAllModels          = new(i_dims, float)
;         copy_VarMeta(r_Temp, r_RegAllModels(0, :, :, :))
;         r_RegAllModels!0        = "Models"
;         r_RegAllModels&Models   = s_ModelSelected
;     end if
;     r_RegAllModels(iModel,:,:,:)= r_Temp
; end do
;
; ;; likely range of Reg 
; printVarSummary(r_RegAllModels) ; [Models] x [Periods] x [lat] x [lon]
; i_dims                         := dimsizes(r_RegAllModels)
; r_LikelyRange                   = r_RegAllModels(:1, :, :, :)
; r_LikelyRange                   = r_LikelyRange@_FillValue
; do iPeriod = 0, i_dims(1)-1
;     do ilat = 0, i_dims(2)-1
;         do ilon = 0, i_dims(3)-1
;             r_Temp             := r_RegAllModels(:, iPeriod, ilat, ilon)
;             i_Ind              := ind(ismissing(r_Temp) .eq. False)
;             r_Temp             := r_Temp(i_Ind)
;             r_Temp             := my_QuaterNum2(r_Temp, 5, 95)
;             r_LikelyRange(:, iPeriod, ilat, ilon) = r_Temp + 0
;         end do
;     end do
; end do
; ; printMinMax(r_LikelyRange, False)
; ;
; r_SNR                           = r_LikelyRange(0, :, :, :)
; r_SNR                           = r_LikelyRange(1, :, :, :) - \
;                                     r_LikelyRange(0, :, :, :)
; r_SNR                           = r_Reg / r_SNR
; i_SNR_Threhold                  = (/1.5, 1./)
; ================================================================
;   2. 绘图，分成三个图绘制
res2                               = True
res2@cnInfoLabelOn                 = False
res2@cnFillOn                      = True
res2@cnLinesOn                     = False
res2@cnLineLabelsOn                = False
res2@gsnDraw                       = False
res2@gsnFrame                      = False
;
; -------------- Set the colorbar
; cmap    = (/ (/103, 0, 31/), (/178, 24, 43/), (/214, 96, 77/), \
;              (/244, 165, 130/), (/253, 219, 219/), \;(/247, 247, 247/), \
;              (/209, 229, 240/), (/146, 197, 222/), (/67, 147, 195/), \
;              (/33, 102, 172/) /) / 255.

; cmap                               = cmap(::-1, :)

; r_cnLinesLevel_c                   = (/0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2/)

r_cnLinesLevel_c                  := (/0.5, 1, 1.25, 1.5, 1.75, 2, 3, 5/)
; if (s_Exp .eq. "ssp126") then
;     r_cnLinesLevel_c              := (/0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3/)
; end if
i_dimColor                         = dimsizes(r_cnLinesLevel_c)

cmap                               = read_colormap_file("cmocean_matter")
i_dimColor                         = dimsizes(cmap(:, 0)) / (i_dimColor+1)
cmap                              := cmap(::i_dimColor, :)
res2@cnFillColors                  = cmap
res2@cnLevelSelectionMode          = "ExplicitLevels"
res2@cnLevels                      = r_cnLinesLevel_c
;
res2@gsnLeftStringFontHeightF      = 0.02
res2@gsnCenterStringFontHeightF    = 0.013
res2@gsnRightStringFontHeightF     = 0.02
res2@gsnLeftStringOrthogonalPosF   = 0.01
res2@gsnRightStringOrthogonalPosF  = 0.01
res2@gsnCenterStringOrthogonalPosF = 0.03

res2@pmTickMarkDisplayMode         = "Always"
res2@tmXBLabelFontHeightF          = 0.013
res2@tmYLLabelFontHeightF          = 0.013

res2@lbLabelBarOn                  = True
res2@lbLabelFontHeightF            = 0.02

; --------------------------------------------------------
    ; set for the map 
    res2@mpMinLatF       = 15.
    res2@mpMaxLatF       = 55.
    res2@mpMinLonF       = 70.
    res2@mpMaxLonF       = 140.

    res2@mpFillOn        = True
    res2@mpDataSetName   = "$WACCEM/my_data/topo/NCL-Chinamap-master/database/Earth..4"
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
opt1@gsnShadeHigh                  = 17
; opt1@gsnShadeLow                   = 17
sres@gsnAddCyclic                  = True
opt1@gsnAddCyclic                  = True

;
s_LeftString                       = (/ "(d)", "(e)", "(f)" /)
plot1                              = new(3, graphic)
plotSig                            = plot1
plotSig99                          = plot1
plotnh                             = plot1
plotnhSig                          = plotnh
plotnhSig_                         = plotnh
plotSigSNR                         = plot1
plotSigSNR1                        = plot1
plotSigSNR_nh                      = plot1
plotSigSNR1_nh                     = plot1
; plotnhSig                          = new((/3, 2/), graphic)
do iSub = 0, 2
    res2@gsnLeftString             = s_LeftString(iSub) + \
                                    " TXx Reg GSAT"
    res2@gsnCenterString           = "" 
    res2@gsnRightString            = s_ExpTitle + ": " + \
                                    i_yrProj(iSub, 0) + "-" + i_yrProj(iSub, 1)
    ;
    plot1(iSub)                    = gsn_csm_contour_map(wks, \
                                                    r_Reg(iSub, :, :), res2)
    
    opt1@gsnShadeHigh              = 17
    delete([/opt1@gsnShadeLow/])
    sres@cnFillColor               = "black"
    r_Temp                        := where(r_SigLevMME(iSub, :, :) \
                                    .ge. r_SigLev_ModelPercentage/100., 1.1, 0)
    copy_VarCoords(r_Reg(iSub, :, :), r_Temp)
    if (s_Exp .eq. "ssp126" .and. iSub .eq. 2) then
        r_Temp                     = where(r_SigLevMME(iSub,:,:) .ge. 50/100., \
                                           1.1, 0)
    end if
    r_Temp(0, 0)                   = 0
    plotSig0                      := gsn_csm_contour(wks, r_Temp, sres)
    plotSig(iSub)                  = gsn_contour_shade(plotSig0, 0, 1., opt1)
    overlay(plot1(iSub), plotSig(iSub))
    ; ;
    sres@cnFillColor               = "white"
    r_Temp                        := where(r_SigLevMME_99(iSub, :, :) .ge. \
                                    r_SigLev_ModelPercentage/100., 1.1, 0)
    if (s_Exp .eq. "ssp126" .and. iSub .eq. 2) then
        r_Temp                     = where(r_SigLevMME_99(iSub,:,:) .ge. 50/100.,\
                                           1.1, 0)
    end if
    r_Temp(0, 0)                   = 0
    copy_VarCoords(r_Reg(iSub, :, :), r_Temp)
    plotSig0                      := gsn_csm_contour(wks, r_Temp, sres)
    plotSig99(iSub)                = gsn_contour_shade(plotSig0, 0, 1., opt1)
    overlay(plot1(iSub), plotSig99(iSub))
    ;
    delete([/opt1@gsnShadeHigh/])
    ; ;; SNR .le. 1.5
    ; r_Temp                        := r_SNR(iSub, :, :)
    ; opt1@gsnShadeLow               = 3
    ; sres@cnFillColor               = "gray"
    ; plotSig0                      := gsn_csm_contour(wks, r_Temp, sres)
    ; plotSigSNR(iSub)               = gsn_contour_shade(plotSig0, \
    ;                                                 i_SNR_Threhold(0), 2, opt1)
    ; overlay(plot1(iSub), plotSigSNR(iSub))
    ; sres@cnFillColor               = "black"
    ; plotSig0                      := gsn_csm_contour(wks, r_Temp, sres)
    ; plotSigSNR1(iSub)              = gsn_contour_shade(plotSig0, \
    ;                                                 i_SNR_Threhold(1), 2, opt1)
    ; overlay(plot1(iSub), plotSigSNR1(iSub))

    l_list                        := my_Add_SouthChinaSea(plot1(iSub), \
                        r_Reg(iSub, :, :), r_SigLevMME(iSub, :, :), \
                        r_SigLevMME_99(iSub, :, :), r_SigLev_ModelPercentage/100., wks)
    if (s_Exp .eq. "ssp126" .and. iSub .eq. 2) then
        l_list                    := my_Add_SouthChinaSea(plot1(iSub), \
                        r_Reg(iSub, :, :), r_SigLevMME(iSub, :, :), \
                        r_SigLevMME_99(iSub, :, :), 50/100., wks)
    end if
    plotnh(iSub)                   = l_list[0]
    g_Temp                        := l_list[1]
    plotnhSig(iSub)                = g_Temp(0)
    plotnhSig_(iSub)               = g_Temp(1)
    ; g_Temp                        := l_list[2]
    ; plotSigSNR_nh(iSub)            = g_Temp(0)
    ; plotSigSNR1_nh(iSub)           = g_Temp(1)
end do 

resP                               = True
resP@gsnMaximize                   = False
gsn_panel(wks, plot1, (/3, 1/), resP)

print(s_FigOut + "")
;
end