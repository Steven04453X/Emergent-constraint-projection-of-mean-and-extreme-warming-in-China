; ========================================================
; 本程序计算指定变量下，中国区域区域平均下预估变化序列，同时输出近、中、远期的空间分布
;   1. 指定模式、变量名、试验、输出文件信息
;   2. 读取CN05.1数据，作为中国陆地区域的标记矩阵
;   3. 逐个模式处理，首先读取数据，然后扣除1995-2014年的平均
;   4. 提取中国区域，然后进行区域平均
;   5. 输出结果，注意区域平均和三个分布都输出
; Usage: ncl s_VarName=\"tas\" 01.China_map.ncl
; ========================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code/my_FunctionInPNNL.ncl"
;
begin
;
; ========================================================
;   1. 指定模式、变量名、试验、输出文件信息
s_ModelSelect             = (/"CNRM-CM6-1", "EC-Earth3", "GFDL-CM4", "MIROC6" /)
; s_ModelSelect             = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
;                        "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
;                        "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
;                         "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
;                         "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
;                         "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
s_4Models                 = "_4Models" ; "_18Models"

s_Exp                     = (/ "ssp585", "ssp126", "ssp245", "ssp370" /)

if (isdefined("s_VarName") .eq. False) then
    s_VarName             = "TXx" ; "pr" ; 
end if 
i_Month                   = ispan(1, 12, 1) ; (/12, 1, 2/) ; (/6, 7, 8/) ;
s_Month                   = "AnMean" ; "DJF" ; "JJA" ;

s_DirOut = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/05.GSAT_scaled_RegionalTXx/06.SpatialPDF/"
system("mkdir -p " + s_DirOut)
s_FileOut                 = s_DirOut + s_VarName + "_ProjectionInChina_" \
                            + s_Month + s_4Models + ".nc"
ls                        = systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileOut               = addfile(s_FileOut, "c")
    GAttr                 = "True"
    GAttr@history         = systemfunc("date") + "  " + systemfunc("pwd") + \
                            "/00.China_map_ForExtremeIndex.ncl"
    GAttr@Month           = i_Month
    fileattdef(FileOut, GAttr)
else 
    FileOut               = addfile(s_FileOut, "w")
end if
;
s_FileOut2                = s_DirOut + s_VarName + "_ProjectionInChinaPattern_"\
                             + s_Month + s_4Models + ".nc"
ls                        = systemfunc("ls " + s_FileOut2)
if (ismissing(ls(0))) then
    FileOutPattern        = addfile(s_FileOut2, "c")
    GAttr                 = "True"
    GAttr@history         = systemfunc("date") + "  " + systemfunc("pwd") + \
                            "/00.China_map_ForExtremeIndex.ncl"
    GAttr@Month           = i_Month
    fileattdef(FileOutPattern, GAttr)
else 
    FileOutPattern        = addfile(s_FileOut2, "w")
end if
; ========================================================
;   2. 读取CN05.1数据，作为中国陆地区域的标记矩阵
s_DirRead                 = "$WACCEM/my_data/Obs/CN05.1/"
s_FileRead                = s_DirRead + "CN05.1_Pre_1961_2015_month_025x025.nc"
FileMask                  = addfile(s_FileRead, "r")
r_Mask                    = FileMask->pre(0, :, :)
; printVarSummary(r_Mask)
; ;; 台湾经纬度范围
; r_lonTaiwan               = (/119, 125/)
; r_latTaiwan               = (/20, 26/)
; r_Mask({r_latTaiwan(0):r_latTaiwan(1)}, {r_lonTaiwan(0):r_lonTaiwan(1)}) = 1
; r_Mask                    = my_mask([/r_Mask/], 0)

i_Mask                    = where(ismissing(r_Mask) .eq. False, 1, 0)
copy_VarCoords(r_Mask, i_Mask)

; ========================================================
;   3. 逐个模式处理，首先读取数据，然后扣除1995-2014年的平均
;
;; Historical 
s_DirRead                := "$WACCEM/my_data/CMIP6/tasmax/POST/TXx/AllReal/"
; s_FileRead               := s_DirRead + s_VarName + "_" + s_Month + \
;                             "_historical_CMIP6.nc"
; FileHist                 := addfile(s_FileRead, "r")
do iModel = 0, dimsizes(s_ModelSelect)-1
    ;
    print(iModel + ": " + s_ModelSelect(iModel))
    s_FileRead           := s_DirRead + s_VarName + "_" + s_Month + \
                            "_historical_" + s_ModelSelect(iModel) + \
                            "_CMIP6_AllReal.nc"
    FileHist             := addfile(s_FileRead, "r")
    ; s_VarName2           := s_VarName + "_" + s_ModelSelect(iModel) + \
    ;                         "_historical"
    ; print(s_VarName2 + "")
    ; r_VarHist            := FileHist->$s_VarName2$({1995:2014}, \
    ;                         {-90:90}, {0:360})
    ; r_VarHist            := dim_avg_n_Wrap(r_VarHist, 0)
    ; ;
    ; s_VarName2           := s_VarName + "_" + s_ModelSelect(iModel) + \
    ;                             "_historical"
    ; FileOutPattern->$s_VarName2$ = r_VarHist
    ;
    do iExp = 0, dimsizes(s_Exp)-1
        print(iExp + ": " + s_Exp(iExp))
        if (s_ModelSelect(iModel) .eq. "ACCESS-ESM1-5" .and. \
            s_Exp(iExp) .eq. "ssp126") then
            print(s_ModelSelect(iModel) + ", " + s_Exp(iExp))
            continue
        end if
        ;
        s_FileRead       := s_DirRead + s_VarName + "_" + s_Month + "_" + \
                            s_Exp(iExp) + "_" + s_ModelSelect(iModel) + \
                            "_CMIP6_AllReal.nc"
        FileSSP          := addfile(s_FileRead, "r")
        ;
        l_list           := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "day")
        s_Real           := l_list[1]
        l_list_          := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "day")
        s_Real_          := l_list_[1]
        if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real)) then
            s_Real       := s_Real_
        end if
        ;
        if (all(ismissing(s_Real)) ) then
            print(s_ModelSelect(iModel) + ", " + s_Exp + " Missing!")
            continue
        end if
        ;
        s_VarNameAll     := getfilevarnames(FileSSP)
        do iReal = 0, dimsizes(s_Real)-1
            s_VarName2   := s_VarName + "_" + s_ModelSelect(iModel) + "_" + \
                            s_Exp(iExp) + "_" + s_Real(iReal)
            i_ind        := ind(s_VarName2 .eq. s_VarNameAll)
            if (ismissing(i_ind(0))) then
                print(s_VarName2 + " Missing!")
                continue 
            end if
            print(s_VarName2 + "")
            r_VarSSP     := FileSSP->$s_VarName2$(:, {-90:90}, {0:360})
            ;
            s_VarName2   := s_VarName + "_" + s_ModelSelect(iModel) + \
                            "_historical_" + s_Real(iReal)
            if (my_VarNameExist(FileHist, s_VarName2) .eq. False) then
                print(s_VarName2 + " Missing!")
                continue
            end if
            r_VarHist    := FileHist->$s_VarName2$({1995:2014}, {-90:90}, \
                                                   {0:360})
            r_VarHist    := dim_avg_n_Wrap(r_VarHist, 0)
        ;
            ;;; 扣气候态
            r_Temp       := conform(r_VarSSP, r_VarHist, (/1, 2/))
            copy_VarCoords(r_VarSSP, r_Temp)
            r_VarSSP      = r_VarSSP - r_Temp
            ; printVarSummary(r_VarSSP)
        ;
            ; ==================================================
            ;   4. 提取中国区域，然后进行区域平均
            ;; 首先提取中国区域
            r_latChina   := (/ min(r_Mask&lat), max(r_Mask&lat) /)
            r_lonChina   := (/ min(r_Mask&lon), max(r_Mask&lon) /)
            r_VarChina   := r_VarSSP(:, {r_latChina(0):r_latChina(1)}, \
                                        {r_lonChina(0):r_lonChina(1)})
            r_latInter   := r_VarChina&lat
            r_lonInter   := r_VarChina&lon
        ;
        ; if (isdefined("i_Mask2_") .eq. False) then
            i_Mask2_     := i_Mask({r_latInter}, {r_lonInter})
            i_Mask2_     := conform(r_VarChina, i_Mask2_, (/1,2/))
            copy_VarCoords(r_VarChina, i_Mask2_)
        ; end if 
        ;
            r_VarChina_  := where(i_Mask2_ .ne. 0, r_VarChina, \
                                  r_VarChina@_FillValue)
            copy_VarCoords(r_VarChina, r_VarChina_)
        ;
            ;; 台湾经纬度范围
            r_lonTaiwan   = (/119, 125/)
            r_latTaiwan   = (/20, 26/)
            r_Temp       := r_VarChina
            r_Temp        = my_mask([/r_Temp/], 0)
            r_VarChina_(:, {r_latTaiwan(0):r_latTaiwan(1)}, \
                        {r_lonTaiwan(0):r_lonTaiwan(1)}) = r_Temp(:, \
                                            {r_latTaiwan(0):r_latTaiwan(1)}, \
                                            {r_lonTaiwan(0):r_lonTaiwan(1)})
        ;; 台湾经纬度范围
        ;
        ;; Hist Clim 
            r_VarHistChina   := r_VarHist({r_latChina(0):r_latChina(1)}, \
                                        {r_lonChina(0):r_lonChina(1)})
            r_Temp           := r_VarHistChina
            r_VarHistChina    = where(i_Mask2_(0, :, :) .ne. 0, r_VarHistChina,\
                                  r_VarHistChina@_FillValue)
            r_Temp            = my_mask([/r_Temp/], 0)
            r_VarHistChina({r_latTaiwan(0):r_latTaiwan(1)}, \
                           {r_lonTaiwan(0):r_lonTaiwan(1)}) = \
                                        r_Temp({r_latTaiwan(0):r_latTaiwan(1)},\
                                        {r_lonTaiwan(0):r_lonTaiwan(1)})
            r_VarHistChina   := my_AreaMean2([/r_VarHistChina/], \
                                    r_VarHistChina&lat, r_VarHistChina&lon, 0)
        ;
        ;;; 中国区域平均
            s_VarName2       := s_VarName + "_ChinaAreaCheck_" + \
                            s_ModelSelect(iModel)+ "_" + s_Exp(iExp) + "_" + \
                            s_Real(iReal)
            FileOut->$s_VarName2$ = r_VarChina_(0, :, :)

            r_VarChina_      := my_AreaMean2([/r_VarChina_/], r_VarChina_&lat, \
                                         r_VarChina_&lon, 0)
            s_VarName2       := s_VarName + "_ChinaArea_" + \
                            s_ModelSelect(iModel) + "_" + s_Exp(iExp) + "_" + \
                            s_Real(iReal)
            r_VarChina_@Clim  = r_VarHistChina
            r_VarChina_!0     = r_VarChina_!0 + "_" + s_Exp(iExp)
            FileOut->$s_VarName2$ = r_VarChina_
        ;
        ;;; 全球分布
            r_VarChina_      := r_VarSSP(:2, :, :)
            r_VarChina_(0,:,:) = dim_avg_n_Wrap(r_VarSSP({2021:2040}, :, :), 0)
            r_VarChina_(1,:,:) = dim_avg_n_Wrap(r_VarSSP({2041:2060}, :, :), 0)
            r_VarChina_(2,:,:) = dim_avg_n_Wrap(r_VarSSP({2080:2099}, :, :), 0)
        ;
            r_VarChina_!0     = "Term"
            r_VarChina_&Term  = ispan(0, 2, 1)
            s_VarName2       := s_VarName + "_" + s_ModelSelect(iModel) + "_" +\
                                s_Exp(iExp) + "_" + s_Real(iReal)
            FileOutPattern->$s_VarName2$= r_VarChina_
                ;
            if (isdefined("r_VarAllReal_ForEachModel") .eq. False .or. \
                iReal .eq. 0) then 
                i_dims          := array_append_record(dimsizes(s_Real), \
                                                       dimsizes(r_VarChina_), 0)
                r_VarAllReal_ForEachModel  := new(i_dims, float)
                copy_VarCoords(r_VarChina_, r_VarAllReal_ForEachModel(0,:,:,:))
                r_VarAllReal_ForEachModel!0 = "Real"
                r_VarAllReal_ForEachModel&Real = s_Real
            end if
            r_VarAllReal_ForEachModel(iReal, :, :, :) = r_VarChina_
        end do 
        ;; For warming STD 
        if (isdefined("r_VarAll") .eq. False) then
            i_dims          := array_append_record(dimsizes(s_Exp), \
                                                   dimsizes(r_VarChina_), 0)
            i_dims          := array_append_record(dimsizes(s_ModelSelect), \
                                                   i_dims, 0)
            r_VarAll        := new(i_dims, float)
                copy_VarCoords(r_VarChina_, r_VarAll(0, 0, :, :, :))
            r_VarAll!0       = "Models"
            r_VarAll&Models  = s_ModelSelect
            r_VarAll!1       = "Exp"
            r_VarAll&Exp     = s_Exp
            r_VarAll@Units   = "K"
        end if 
        r_VarAll(iModel, iExp, :, :, :) = dim_avg_n_Wrap(\
                                                r_VarAllReal_ForEachModel, 0) 
        delete([/r_VarAllReal_ForEachModel/])
        ;
    end do
    print(s_ModelSelect(iModel) + " Done!")
end do
;
;; Warming STD 
r_VarAll                    := dim_stddev_n_Wrap(r_VarAll, 0)
s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/05.GSAT_scaled_RegionalTXx/06.SpatialPDF/"
s_FileOut                    = s_DirOut + "TXx_ProjectedInterModelSTD_" + \
                            s_Month + s_4Models + ".nc"
ls                           = systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileSTD                  = addfile(s_FileOut, "c")
    GAttr                    = True
    GAttr@history            = systemfunc("date") + "  " + systemfunc("pwd") + \
                            "/00.China_map_ForExtremeIndex.ncl"
    GAttr@Des                = "inter-model standard deviation in the projection"
    fileattdef(FileSTD, GAttr)
else 
    FileSTD                  = addfile(s_FileOut, "w")
end if

do iExp = 0, dimsizes(s_Exp)-1
    s_VarName2               = "STD_" + s_VarName + "_" + s_Exp(iExp)
    r_Temp                  := r_VarAll(iExp, :, :, :)
    r_Temp@long_name         = "inter-model STD of " + s_VarName\
                                + " increase in " + s_Exp(iExp)
    FileSTD->$s_VarName2$    = r_Temp
end do

;
end






