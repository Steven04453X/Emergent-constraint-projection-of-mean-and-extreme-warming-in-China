; =================================================================
; Author: Ziming Chen 
; This code is for:
;   - vis the regression btw present-day residual trend and constrained warming in the projection 
; =================================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
undef("my_Add_SouthChinaSea")
function my_Add_SouthChinaSea(map, Draw:numeric, DrawSigLev:numeric, wks)
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
; r_DrawSigLev_99                     = DrawSigLev_99
; r_SigLev                            = SigLev
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
    sres@cnFillDotSizeF             = 0.002
    opt1        = True
    opt1@gsnShadeFillType           = "pattern"
    opt1@gsnShadeHigh               = 17
    ; opt1@gsnShadeLow                = 17
    sres@gsnAddCyclic               = True
    opt1@gsnAddCyclic               = True

    plotSig                         = new(2, graphic)
    sres@cnFillColor               := "white"
    r_Temp                         := where(r_DrawSigLev .le. 0.1, 1.1, 0)
    copy_VarCoords(r_DrawSigLev, r_Temp)
    r_Temp(0, 0)                    = 0
    map_stipple                     = gsn_csm_contour(wks, r_Temp, sres)
    plotSig(0)                      = gsn_contour_shade(map_stipple, 0, 1.,opt1)
    overlay(map_nanhai, plotSig(0))
    ;
    ; sres@cnFillColor               := "white"
    ; r_Temp                         := where(r_DrawSigLev .le. 0.01, 1.1, 0)
    ; copy_VarCoords(r_DrawSigLev, r_Temp)
    ; r_Temp(0, 0)                    = 0
    ; map_stipple                     = gsn_csm_contour(wks, r_Temp, sres)
    ; plotSig(1)                      = gsn_contour_shade(map_stipple, 0, 1.,opt1)
    ; overlay(map_nanhai, plotSig(1))
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
begin
;
s_VarName       = "RegInModelDim"
s_Exp           = "ssp370" ; "ssp585" ; "ssp245" ; "ssp126" ; 
s_ExpTitle      = "SSP3-7.0" ; "SSP5-8.5" ; "SSP2-4.5" ; "ssp1-2.6" ; 
s_4Models       = "_23Models"
i_yrProj        = (/ (/2021, 2040/), (/2041, 2060/), (/2080, 2099/) /)
;
s_DirOut = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Fig/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/03.ConstraintByUsingLocalTrend/"
system("mkdir -p " + s_DirOut)
s_FigOut        = s_DirOut + "05.VisRegressionBtwLocalRTrend_ConWarming_" + \
                s_Exp + s_4Models
wks             = gsn_open_wks("X11", s_FigOut)
;
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/03.ConstraintByUsingLocalTrend/"
s_FileRead      = s_DirRead + "Regression_LocalTrendAndProjectedWarming_" + \
                s_Exp + s_4Models + ".nc"
FileRead        = addfile(s_FileRead, "r")
;
s_VarName2      = s_VarName + "_" + s_Exp
r_Reg           = FileRead->$s_VarName2$
;
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
; r_cnLinesLevel_c                   = (/0.2, 0.4, 0.6, 0.8, 1., 1.5, 2, 3/)
; r_cnLinesLevel_c                  := (/0.5, 1, 1.25, 1.5, 1.75, 2, 3, 5/)
r_cnLinesLevel_c                   = (/-3, -2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2, 3/)
if (s_Exp .eq. "ssp585" .or. s_Exp .eq. "ssp370") then
    r_cnLinesLevel_c              := (/-5, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 5/)
end if
i_dimColor                         = dimsizes(r_cnLinesLevel_c)

; cmap                               = read_colormap_file("cmocean_matter")
; i_dimColor                         = dimsizes(cmap(:, 0)) / (i_dimColor+1)
; cmap                              := cmap(::i_dimColor, :)
cmap  = (/(/103, 0, 31/), (/178, 24, 43/), (/214, 96, 77/), (/244, 165, 130/), \
          (/253, 219, 199/), (/247, 247, 247/), (/209, 229, 240/), \
          (/146, 197, 222/), (/67, 147, 195/), (/33, 102, 172/), (/5, 48, 97/)/) / 255.
cmap                               = cmap(::-1, :)

res2@cnFillColors                  = cmap
res2@cnLevelSelectionMode          = "ExplicitLevels"
res2@cnLevels                      = r_cnLinesLevel_c
;
res2@gsnLeftStringFontHeightF      = 0.02
res2@gsnCenterStringFontHeightF    = 0.013
res2@gsnRightStringFontHeightF     = 0.018
res2@gsnLeftStringOrthogonalPosF   = 0.09
res2@gsnRightStringOrthogonalPosF  = 0.01
res2@gsnCenterStringOrthogonalPosF = 0.03

res2@pmTickMarkDisplayMode         = "Always"
res2@tmXBLabelFontHeightF          = 0.013
res2@tmYLLabelFontHeightF          = 0.013
; res2@gsnAddCyclic                  = True

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
sres@cnFillDotSizeF                = 0.002
opt1        = True
opt1@gsnShadeFillType              = "pattern"
opt1@gsnShadeHigh                  = 17
; opt1@gsnShadeLow                   = 17
; sres@gsnAddCyclic                  = True
; opt1@gsnAddCyclic                  = True

;
s_LeftString                       = (/ "(a)", "(b)", "(c)" /)
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
                    " Trend~F10B~res~NF21~ Reg Constrained Mean Warming across Models"
    res2@gsnCenterString           = "" 
    res2@gsnRightString            = s_ExpTitle + ": " + \
                                    i_yrProj(iSub, 0) + "-" + i_yrProj(iSub, 1)
    ;
    plot1(iSub)                    = gsn_csm_contour_map(wks, \
                                                    r_Reg(0, iSub, :, :), res2)
    
    opt1@gsnShadeHigh              = 17
    delete([/opt1@gsnShadeLow/])
    sres@cnFillColor               = "white"
    r_Temp                        := where(r_Reg(1, iSub, :, :).le.0.1, 1.1, 0)
    copy_VarCoords(r_Reg(1, iSub, :, :), r_Temp)
    r_Temp(0, 0)                   = 0
    plotSig0                      := gsn_csm_contour(wks, r_Temp, sres)
    plotSig(iSub)                  = gsn_contour_shade(plotSig0, 0, 1., opt1)
    overlay(plot1(iSub), plotSig(iSub))
    ; ;
    ; sres@cnFillColor               = "white"
    ; r_Temp                        := where(r_Reg(1, iSub, :, :).le.0.01, 1.1, 0)
    ; r_Temp(0, 0)                   = 0
    ; copy_VarCoords(r_Reg(1, iSub, :, :), r_Temp)
    ; plotSig0                      := gsn_csm_contour(wks, r_Temp, sres)
    ; plotSig99(iSub)                = gsn_contour_shade(plotSig0, 0, 1., opt1)
    ; overlay(plot1(iSub), plotSig99(iSub))
    ;
    delete([/opt1@gsnShadeHigh/])

    l_list                        := my_Add_SouthChinaSea(plot1(iSub), \
                                r_Reg(0, iSub, :, :), r_Reg(1, iSub, :, :), wks)
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