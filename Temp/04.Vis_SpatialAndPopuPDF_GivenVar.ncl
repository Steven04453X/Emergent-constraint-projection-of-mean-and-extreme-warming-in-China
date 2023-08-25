; ===========================================================
; 本程序绘制指定空间PDF和人口PDF的结果
;   1. 读取数据，注意先读取约束前和约束后的MME的空间PDF和人口PDF
;   2. 参考之前的程序，进行绘图，绘制成近、中、远三个时期的结果，上中下结构，然后右侧是人口PDF的结果
; ===========================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
undef("my_DoublePDF_HistogramAndShading")
function my_DoublePDF_HistogramAndShading(wks, PDF1, PDF2, resO)
; 本程序绘制输入的两条PDF曲线，注意PDF1叠加在PDF2上面
; 两个PDF上，第一维包含三个元素，分别是MME，不确定性的下界和上界
; 特殊的res:
;   res2@my_LegendLabels：图例的名称
;   res2@my_trXMinMaxF_trYMinMaxF：x和y轴的范围，显示x轴，然后是y轴，从小到大的一维数组
;   res2@my_DrawLegend: 控制是否绘制图中的Legend
;   res2@my_LegendLabelsPos: 控制legend的位置，先竖直位置然后水平
;   r_PDF2@my_percentage和r_PDF1@my_percentage: 显著的PDF面积或人口占比
local r_PDF1, r_PDF2, res2 
begin
;
r_PDF1                              = PDF1
r_PDF2                              = PDF2
res2                                = resO
if (res2 .eq. False) then
    res2                            = True
end if

res2                                = True
res2@gsnDraw                        = False
res2@gsnFrame                       = False
res2@xyLineThicknesses              = (/3.0, 3./)
res2@xyDashPatterns                 = 0

;;; 梯级图
res2@gsnXYBarChart                  = True
res2@gsnXYBarChartOutlineOnly       = True
res2@gsnXYBarChartBarWidth          = 0.2
res2@gsnYRefLine                    = 0.
res2@gsnYRefLineDashPattern         = 1
;;; 梯级图
    
if (ismissing(res2@tiXAxisString) .eq. True) then 
    res2@tiXAxisString              = "Precip. Changes (%)"
end if 
if (ismissing(res2@tiYAxisString) .eq. True) then 
    res2@tiYAxisString              = "Land Fraction (%)"
end if 
res2@tiYAxisFontHeightF             = 0.011
res2@tiXAxisFontHeightF             = 0.011

if (ismissing(res2@gsnLeftString) .eq. True) then 
    res2@gsnLeftString              = "(a) " 
end if 
if (ismissing(res2@gsnCenterString) .eq. True) then 
    res2@gsnCenterString            = ""
end if 
if (ismissing(res2@gsnRightString) .eq. True) then 
    res2@gsnRightString             = ""
end if 
res2@gsnCenterString                = ""
res2@gsnLeftStringFontHeightF       = 0.015
; res2@gsnLeftStringParallelPosF      = -0.15
; res2@gsnLeftStringOrthogonalPosF    = 0.12
; res2@gsnRightStringFontHeightF      = 0.012
    
if ( all(ismissing(res2@my_trXMinMaxF_trYMinMaxF)) .eq. True) then 
    res2@trXMinF                    = -10
    res2@trXMaxF                    = 50
    res2@trYMinF                    = 0
    res2@trYMaxF                    = 5
else
    r_Temp                         := res2@my_trXMinMaxF_trYMinMaxF
    res2@trXMinF                    = r_Temp(0)
    res2@trXMaxF                    = r_Temp(1)
    res2@trYMinF                    = r_Temp(2)
    res2@trYMaxF                    = r_Temp(3)
end if 
    ;
res2@tmXTOn                         = False
res2@tmYROn                         = False
res2@tmXTBorderOn                   = False
res2@tmYRBorderOn                   = False

res2@vpHeightF                      = 0.2
res2@vpWidthF                       = 0.25
res2@tmXBLabelFontHeightF           = 0.010
res2@tmYLLabelFontHeightF           = 0.010
; res2@tmYLFormat                     = "f"            

if (ismissing(res2@vpYF) .eq. True) then 
    res2@vpYF                       = 0.8
end if 
if (ismissing(res2@vpXF) .eq. True) then 
    res2@vpXF                       = 0.15
end if 

plot                                = new(6, graphic)

x = r_PDF1&$r_PDF1!0$

;;; 约束前的PDF分布
r_TempColor                         = new((/2, 3/), float)
if ( all(ismissing(res2@PDF2LineColor)) .eq. True) then 
    r_TempColor(0, :)               = (/153, 0, 2/) / 255.
else 
    r_TempColor(0, :)               = res2@PDF1LineColor
end if 
if ( all(ismissing(res2@PDF1LineColor)) .eq. True) then 
    r_TempColor(1, :)               = (/0, 52, 102/) / 255.
else 
    r_TempColor(1, :)               = res2@PDF2LineColor
end if 
res2@xyLineColors                   = r_TempColor ; (/(/0, 0, 0/), (/112, 160, 205/)/) /255.
r_PDF2!0                            = "bins"
x                                   = r_PDF2&bins
i_dims                             := array_append_record(2, dimsizes(r_PDF1), 0)
r_DrawArr                           = conform_dims(i_dims, r_PDF1, (/1/))
copy_VarCoords(r_PDF2, r_DrawArr(0, :))
r_DrawArr(1, :)                     = r_PDF2 ; 
plot(0)                             = gsn_csm_xy(wks, x, r_DrawArr, res2)

res2@gsnXYBarChart                  = True
res2@gsnLeftString                  = ""
res2@gsnRightString                 = ""
res2@gsnCenterString                = ""
res2@gsnXYBarChartOutlineOnly       = True ; False
res2@gsnAboveYRefLineColor          = (/67, 147, 195/)/255.  ;
res2@gsnXYBarChartFillOpacityF      = 0.6
res2@gsnXYBarChartBarWidth          = 12
res2@xyLineColors                  := (/ "Transparent" /)
;
if (all(ismissing(res2@PDF2ShadingColor)) .eq. True) then 
    res2@gsnAboveYRefLineColor      = (/67, 147, 195/) / 255. ; (/146, 197, 222/) /255. 
else 
    res2@gsnAboveYRefLineColor      = res2@PDF2ShadingColor
end if 
dbin                                = r_PDF2&bins(1) - r_PDF2&bins(0)
plot(1)                             = gsn_csm_xy(wks, x+dbin/2, r_PDF2, res2)
overlay(plot(0), plot(1))
res2@gsnXYBarChartOutlineOnly       = False
x2                                 := x
if (all(ismissing(res2@PDF_Thredhold)) .eq. False) then 
    r_PDF2                          = where(r_PDF2&bins.lt.res2@PDF_Thredhold,\
                                            r_PDF2@_FillValue, r_PDF2)
    i_Ind                           = ind(r_PDF2&bins .ge. res2@PDF_Thredhold)
    x2                             := x2(i_Ind(0):)
    r_PDF2                         := r_PDF2(i_Ind(0):)
end if 
plot(4)                             = gsn_csm_xy(wks, x2+dbin/2, r_PDF2, res2)
; plot(4)                             = gsn_csm_xy(wks, x2, r_PDF2, res2)
overlay(plot(0), plot(4))


res2@gsnXYBarChartOutlineOnly       = True
if (all(ismissing(res2@PDF1ShadingColor)) .eq. True) then 
    res2@gsnAboveYRefLineColor     := (/252, 209, 197/) / 255.  ;
else 
    res2@gsnAboveYRefLineColor     := res2@PDF1ShadingColor
end if 
res2@gsnXYBarChartFillOpacityF      = 0.5
plot(2)                             = gsn_csm_xy(wks, x+dbin/2, r_PDF1, res2)
overlay(plot(0), plot(2))
res2@gsnXYBarChartOutlineOnly       = False
x2                                 := x
if (all(ismissing(res2@PDF_Thredhold)) .eq. False) then 
    r_PDF1                          = where(r_PDF1&bins.lt.res2@PDF_Thredhold,\
                                            r_PDF1@_FillValue, r_PDF1)
    i_Ind                           = ind(r_PDF1&bins .ge. res2@PDF_Thredhold)
    x2                             := x2(i_Ind(0):)
    r_PDF1                         := r_PDF1(i_Ind(0):)
end if 
plot(5)                             = gsn_csm_xy(wks, x2+dbin/2, r_PDF1, res2)
; plot(5)                             = gsn_csm_xy(wks, x2, r_PDF1, res2)
overlay(plot(0), plot(5))


if ( ismissing(res2@my_DrawLegend) .eq. True .or. res2@my_DrawLegend .eq. True ) then 
; Create legend
    lgres                           = True
    lgres@lgDashIndexes             = (/0, 0/) ;res2@xyDashPatterns
    lgres@lgMonoDashIndex           = False           ; Same line dash patterns 
    lgres@lgLineThicknessF          = (/4, 4/) 
    lgres@lgLabelFontHeightF        = .2 ; .14           ; set the legend label font thickness
    lgres@vpWidthF                  = 0.2          ; width of legend (NDC)
    lgres@vpHeightF                 = 0.06            ; height of legend
    lgres@lgPerimOn                 = False          ; draw the box perimeter in orange
    lgres@lgLineColors              = r_TempColor ; (/(/0, 0, 0/), (/153, 0, 2/) /)  / 255.
    ;res2@xyLineColors

    if ( all(ismissing(res2@my_LegendLabels)) .eq. True) then
        labels = "  " + (/ "PDF1", "PDF2" /)
    else 
        labels = "  " + res2@my_LegendLabels
    end if
    lbid = gsn_create_legend(wks, 2, labels, lgres)
;
    ; Add legend to bottom right corner of plot
    amres                           = True
    amres@amJust                    = "TopLeft"  ; Corner of plot for positioning
    if (all(ismissing(res2@my_LegendLabelsPos)) .eq. True) then 
        amres@amOrthogonalPosF      =  -0.4   ;  0.5 is the bottom edge of the plot.
        amres@amParallelPosF        =  -0.1   ;  0.5 is the right edge of the plot.
    else 
        r_Temp                     := res2@my_LegendLabelsPos
        amres@amOrthogonalPosF      = r_Temp(0)
        amres@amParallelPosF        = r_Temp(1)
    end if 
    annoid  = gsn_add_annotation(plot(0), lbid, amres)
end if 

; ;;;; PDF对应的面积和人口比例的阈值线
; if (all(ismissing(res2@PDF_Thredhold)) .eq. False) then
;     plres                       = True
;     plres@gsLineThicknessF      = 4
;     plres@gsLineDashPattern     = 1
;     ;
;     r_Line                      = (/ 0, res2@trYMaxF /)
;     r_Line_x                    = (/ res2@PDF_Thredhold, res2@PDF_Thredhold /)
;     plres@gsLineColor           = "black"
;     plot@$tostring("Thredhold")$= gsn_add_polyline(wks, plot(0), r_Line_x, \
;                                                    r_Line, plres)
; end if

;;;; PDF对应的面积和人口比例
s_After                         = sprintf("%2.0f", r_PDF1@my_percentage) + "%"
s_Before                        = sprintf("%2.0f", r_PDF2@my_percentage) + "%"
; s_After                         = s_After + "          "
r_TempColor                     = (/(/153, 0, 2/), (/0, 52, 102/) /)/ 255.
;; 通过xy的坐标添加
plotText                        = new(2, graphic)
txres                           = True
txres@txFontHeightF             = 0.015

r_TempLegend                    = (res2@trXMaxF - res2@trXMinF) / 2. ;+ res2@trXMinF
r_TempLegend2                   = abs(r_TempLegend) / 2.
if (all(ismissing(res2@PDF1LegendLocation)) .eq. True) then 
    ;; (x, y)
    xytxt1                      = (/ res2@trXMinF + r_TempLegend2, \
                                    res2@trYMaxF * 2. /) / 2.
else 
    xytxt1                      = res2@PDF1LegendLocation
end if 
if (all(ismissing(res2@PDF2LegendLocation)) .eq. True) then 
    ;; (x, y)
    xytxt2                      = xytxt1
    xytxt2(0)                   = res2@trXMinF + r_TempLegend2*3
else 
    xytxt2                      = res2@PDF2LegendLocation
end if 
txres@txFontColor               = r_TempColor(1, :)
plotText(1)                     = gsn_add_text(wks, plot(0), s_Before, \
                                               xytxt2(0), xytxt2(1), txres)
txres@txFontColor               = r_TempColor(0, :)
; plotText(1)                     = gsn_add_text(wks, plot(0), s_After, \
;                                                xytxt2(0), xytxt2(1), txres)
plotText(0)                     = gsn_add_text(wks, plot(0), s_After, \
                                               xytxt1(0), xytxt1(1), txres)

return(plot)
;
end
;
begin
;
; ===========================================================
;   1. 读取数据，注意先读取约束前和约束后的MME的空间PDF
s_4Models               = "_20Models"
s_ExceedingLocalSTD     = "_ExceedingLocalSTDFraction"

s_DirRead = "/WORK/chenzm/Work/OtherWork/2020/ProjectionInChina/Data/05.GSAT_scaled_RegionalTXx/06.SpatialPDF/"
s_FileRead              = s_DirRead + "SpatialPDF_TXx" + \
                            s_ExceedingLocalSTD + s_4Models + ".nc"
FileRead                = addfile(s_FileRead, "r")

s_FileRead              = s_DirRead + "PopuPDF_TXx" + \
                            s_ExceedingLocalSTD + s_4Models + ".nc"
FilePopu                = addfile(s_FileRead, "r")

s_Exp                   = "ssp245" ; "ssp585" ; 
s_ExpTitle              = "SSP2-4.5" ; "SSP5-8.5" ; 
s_VarName2              = "SpatialPDF_TXx_Constrained_MME_" + s_Exp
r_Constrained           = FileRead->$s_VarName2$
r_SigArea_Constrained   = r_Constrained@SigArea

s_VarName2              = "SpatialPDF_TXx_MME_" + s_Exp
r_Raw                   = FileRead->$s_VarName2$
r_SigArea_Raw           = r_Raw@SigArea
;
;; Population 
s_VarName2              = "PopuPDF_TXx_Constrained_MME_" + s_Exp
r_Constrained_Popu      = FilePopu->$s_VarName2$
r_SigPopu_Constrained   = r_Constrained_Popu@SigPopu 

s_VarName2              = "PopuPDF_TXx_MME_" + s_Exp
r_Raw_Popu              = FilePopu->$s_VarName2$
r_SigPopu_Raw           = r_Raw_Popu@SigPopu 

s_DirOut = "/WORK/chenzm/Work/OtherWork/2020/ProjectionInChina/Fig/05.GSAT_scaled_RegionalTXx/06.SpatialPDF/"
system("mkdir -p " + s_DirOut)
s_FigOut                = s_DirOut + "04.Vis_SpatialAndPopuPDF_GivenVar_" + \
                        s_Exp + s_ExceedingLocalSTD
wks                     = gsn_open_wks("X11", s_FigOut)

; ===========================================================
;   2. 参考之前的程序，进行绘图，绘制成近、中、远三个时期的结果，上中下结构
res                    := True
res@my_DrawLegend       = False

r_TrXMinMax_YMinMax     = (/ (/0,   2.5,   0, 25/), \
                             (/0.5, 3.5,   0, 20/), \
                             (/2,   7,     0, 15/)/) * 1.
r_TrXMinMax_YMinMax_Popu= (/ (/0,   2.5,   0, 30/), \
                             (/0.5, 3.5,   0, 25/), \
                             (/2,   7,     0, 20/)/) * 1.

if (s_Exp .eq. "ssp245") then
    r_TrXMinMax_YMinMax     = (/ (/0,   2.5,   0, 25/), \
                             (/0.5, 3.5,   0, 20/), \
                             (/  1,   4,     0, 24/)/) * 1.
    r_TrXMinMax_YMinMax_Popu= (/ (/0,   2.5,   0, 30/), \
                             (/0.5, 3.5,   0, 25/), \
                             (/  1,   4,     0, 28/)/) * 1.
end if

s_LeftString            = (/ "(a)", "(b)", "(c)" /)
s_LeftString2           = (/ "(d)", "(e)", "(f)" /)
s_Term                  = (/ "Near", "Mid", "Long" /) + "-term Projection in " \
                        + s_ExpTitle
res@tiXAxisString       = "TXx Changes (K)"

res@PDF_Thredhold       = 0

res@PDF1LegendLocation := (/7, 25/)
res@PDF2LegendLocation := (/8, 20/)

r_Position              = (/ (/-0.6, 0.05/), \
                             (/-0.6, 0.05/), \
                             (/-0.6, 0.05/) /)

r_Constrained@my_percentage = 0
r_Raw@my_percentage     = 0

r_vpYF                  = (/0.92, 0.6, 0.28/)

r_RunAveDot             = (/3, 3, 3/)
r_RunAveDot2            = (/3, 3, 3/)
plots                   = new((/3, 6/), graphic)
plotsPopu               = plots
do iTerm = 0, 2
    ;
    if (r_RunAveDot(iTerm) .ne. 0) then 
        r_Constrained(iTerm, :) = runave_n_Wrap(r_Constrained(iTerm, :), \
                                            r_RunAveDot(iTerm), 0, 0)
        r_Raw(iTerm, :)         = runave_n_Wrap(r_Raw(iTerm, :), \
                                            r_RunAveDot(iTerm), 0, 0)
    end if 
    ;
    res@my_DrawLegend   = True
    res@my_LegendLabels = (/ "Constrained: " + \
                                toint(r_SigArea_Constrained(iTerm)) + "%", \
                             "Unconstrained: " + toint(r_SigArea_Raw(iTerm)) + \
                             "%" /)
    res@my_LegendLabelsPos = r_Position(iTerm, :)

    res@gsnLeftString   = s_LeftString(iTerm) + " Area" ;+ s_Term(iTerm)
    res@tiYAxisString   = "Land Fraction (%)"
    res@vpYF            = r_vpYF(iTerm)
    res@my_trXMinMaxF_trYMinMaxF= r_TrXMinMax_YMinMax(iTerm, :)
    res@vpXF            = 0.1
    plots(iTerm, :)     = my_DoublePDF_HistogramAndShading(wks, \
                                r_Constrained(iTerm, :), r_Raw(iTerm, :), res)
    draw(plots(iTerm, 0))
    ;
    ;; Population
    if (r_RunAveDot2(iTerm) .ne. 0) then 
        r_Constrained_Popu(iTerm,:) = runave_n_Wrap(r_Constrained_Popu(iTerm,:),\
                                            r_RunAveDot2(iTerm), 0, 0)
        r_Raw_Popu(iTerm, :)        = runave_n_Wrap(r_Raw_Popu(iTerm, :), \
                                            r_RunAveDot2(iTerm), 0, 0)
    end if 
    res@vpXF            = 0.55
    res@tiYAxisString   = "Population Fraction (%)"
    res@gsnLeftString   = s_LeftString2(iTerm) + " Population"
    res@my_trXMinMaxF_trYMinMaxF= r_TrXMinMax_YMinMax_Popu(iTerm, :)
    ;
    res@my_DrawLegend   = True
    ; res@my_LegendLabels = (/ "Constrained", "Unconstrained" /)
    res@my_LegendLabels = (/ "Constrained: " + \
                                toint(r_SigPopu_Constrained(iTerm)) + "%", \
                             "Unconstrained: " + toint(r_SigPopu_Raw(iTerm)) + \
                             "%" /)
    res@my_LegendLabelsPos = r_Position(iTerm, :)
    ;
    plotsPopu(iTerm,:)  = my_DoublePDF_HistogramAndShading(wks, \
                        r_Constrained_Popu(iTerm, :), r_Raw_Popu(iTerm, :), res)
    ; res@my_DrawLegend   = False
    ;
    ;; main string 
    txres               = True         ; Text resources
    txres@txFontHeightF = 0.015
    txres@txJust        = "CenterLeft"
    labels              = s_Term(iTerm)
    gsn_text_ndc(wks, labels, 0.32, r_vpYF(iTerm) + 0.04, txres)
    ;
    draw(plotsPopu(iTerm, 0))
end do

frame(wks)
print(s_FigOut + "")
;
end