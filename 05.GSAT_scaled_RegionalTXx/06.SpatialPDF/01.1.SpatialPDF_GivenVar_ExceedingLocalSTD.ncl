; ==================================================================
; 本程序计算输入变量，中国区域内的空间PDF分布，并计算超过局地升温标准差的面积占比
;   1. 指定模式、试验、变量名、三个时期的名称
;   2. 读取三个预估时期的升温标准差
;   3. 逐个模式、逐个试验进行，首先读取对应的异常场，然后挑选出中国区域
;   4. 首先计算整个中国的面积
;   5. 分别在三个预估时期，计算其空间的PDF，随后挑选出局地超过升温标准差的格点，并计算其面积占比
; 注意上述过程，同时对约束后的结果进行
; 此外，挑选出升温幅度超过局地升温标准差的面积占比出来
; ==================================================================
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

begin
;
; ==================================================================
;   1. 指定模式、试验、变量名、三个时期的名称
s_ModelSelected = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
; s_ModelSelected = (/"CNRM-CM6-1", "EC-Earth3", "GFDL-CM4", "MIROC6" /)
s_Exp           = (/"ssp126", "ssp245", "ssp370", "ssp585"/)
s_VarName       = "TXx"
s_Term          = (/"Near", "Mid", "Long"/)
s_Month         = "AnMean"
s_4Models       = "_18Models" ; "_4Models" ; 

i_STD_Times     = 2

s_DirRead       = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/05.GSAT_scaled_RegionalTXx/06.SpatialPDF/"
s_FileRead      = s_DirRead + "TXx_ProjectionInChinaPattern_AnMean" + s_4Models\
                  + ".nc"
FileRaw         = addfile(s_FileRead, "r")
;
;; Output 
s_DirOut = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/05.GSAT_scaled_RegionalTXx/06.SpatialPDF/"
system("mkdir -p " + s_DirOut)
s_FileOut       = s_DirOut + "SpatialPDF_" + s_VarName + \
                "_ExceedingLocalSTDFraction" + s_4Models + ".nc"
system("rm -rf " + s_FileOut)
ls              = systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileOut     = addfile(s_FileOut, "c")
    GAttr       = True
    GAttr@history = systemfunc("date") + "  " + systemfunc("pwd") + \
                    "/01.1.SpatialPDF_GivenVar_ExceedingLocalSTD.ncl"
    fileattdef(FileOut, GAttr)
else 
    FileOut     = addfile(s_FileOut, "w")
end if
;
; ==================================================================
;   2. 读取三个预估时期的升温标准差
s_DirRead       = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/05.GSAT_scaled_RegionalTXx/06.SpatialPDF/"
s_FileRead      = s_DirRead + "TXx_ProjectedInterModelSTD_AnMean" + s_4Models +\
                  ".nc"
FileSTD         = addfile(s_FileRead, "r")
;
; ==================================================================
;   3. 逐个模式、逐个试验进行，首先读取对应的异常场，然后挑选出中国区域
do iExp = 0, dimsizes(s_Exp)-1
    ;
    ;; STD 
    s_VarName2  = "STD_TXx_" + s_Exp(iExp)
    r_STD      := FileSTD->$s_VarName2$(:, {10:60}, {65:140})
    r_STD      := my_SelChina([/r_STD/])
    if (isdefined("i_STD_Times")) then
        r_STD   = r_STD * i_STD_Times
    end if
    ;
    ; ==================================================================
    ;   4. 首先计算整个中国的面积
    r_Area     := my_GridArea(r_STD&lat, r_STD&lon)
    r_Area     := conform(r_STD, r_Area, (/1, 2/))
    copy_VarCoords(r_STD, r_Area)
    r_Area      = where(ismissing(r_STD), r_Area@_FillValue, r_Area)
    r_AreaTotal:= dim_sum_n_Wrap(r_Area, (/1, 2/))
    print("Total Area: " + r_AreaTotal(0) + " m2")
    print("")
    ;
    do iModel = 0, dimsizes(s_ModelSelected)-1
        ;
        if (s_ModelSelected(iModel) .eq. "ACCESS-ESM1-5" .and. \
            s_Exp(iExp) .eq. "ssp126") then
            print(s_ModelSelected(iModel) + ", " + s_Exp(iExp) + ", continue")
            continue
        end if
        ;
        l_list         := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp(iExp), "day")
        s_Real         := l_list[1]
        l_list_        := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp(iExp), "day")
        s_Real_        := l_list_[1]
        if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real)) then
            s_Real     := s_Real_
        end if
        ;
        if (all(ismissing(s_Real)) ) then
            print(s_ModelSelected(iModel) + ", " + s_Exp + " Missing!")
            continue
        end if
        ;
        do iReal = 0, dimsizes(s_Real)-1
            s_VarName2      = "TXx_" + s_ModelSelected(iModel) + "_" + \
                                s_Exp(iExp) + "_" + s_Real(iReal)
            L_Test          = my_VarNameExist(FileRaw, s_VarName2)
            if (L_Test .eq. False) then
                print(s_VarName2 + "  Missing!")
                continue
            end if
            r_Var          := FileRaw->$s_VarName2$(:, {10:60}, {65:140})
        ;
            ;; China Area 
            r_Var          := my_SelChina([/r_Var/])
        ;;
            ; ==================================================================
            ;   3. 分别在三个预估时期，计算其空间的PDF
            i_Nbin          = 50
            res             = True
            res@bin_min     = -1.
            res@bin_max     = 7.
            do iTerm = 0, 2
                l_list     := my_SpatialPDF([/r_Var(iTerm, :, :)/], i_Nbin, res)
                r_Spatial  := l_list[0]
                ;
                if (iTerm .eq. 0) then
                    i_dims                := (/3, dimsizes(r_Spatial)/)
                    r_SpatialPDF           = new(i_dims, float)
                    copy_VarCoords(r_Spatial, r_SpatialPDF(0, :))
                    copy_VarCoords(r_Var(:, 0, 0), r_SpatialPDF(:, 0))
                    r_SpatialPDF@long_name = "Spatial PDF of " + \
                                s_ModelSelected(iModel) + " in " + s_Exp(iExp)
                end if
                r_SpatialPDF(iTerm, :)     = r_Spatial
            end do
            ;; The Area Fraction where warming exceed local STD (* i_STD_Times)
            r_AreaSig      := r_Area
            r_AreaSig       = where(abs(r_Var) .ge. r_STD, r_AreaSig, \
                                r_AreaSig@_FillValue)
            r_AreaSig      := dim_sum_n_Wrap(r_AreaSig, (/1, 2/))
            r_AreaSig       = r_AreaSig / r_AreaTotal * 100.
        ; print(r_AreaSig + "")
        ; print("")
        ;
            r_SpatialPDF@SigArea      = r_AreaSig
            r_SpatialPDF@SigAreaUnits = "%"
            r_SpatialPDF@SigDefinition= "one time local STD "
            if (isdefined("i_STD_Times")) then
                r_SpatialPDF@SigDefinition= i_STD_Times + " time local STD "
            end if
            ;;;;; Output 
            s_VarName2      = "SpatialPDF_" + s_VarName + "_" + \
                                s_ModelSelected(iModel) + "_" + s_Exp(iExp) + \
                                "_" + s_Real(iReal)
            FileOut->$s_VarName2$          = r_SpatialPDF
            ;;;;; Output 
            ;
            if (iReal .eq. 0 .or. \
                isdefined("r_SpatialPDF_AllReal") .eq. False) then
                i_dims     := array_append_record(dimsizes(s_Real), \
                                                  dimsizes(r_SpatialPDF), 0)
                r_SpatialPDF_AllReal  := new(i_dims, float)
                copy_VarCoords(r_SpatialPDF, r_SpatialPDF_AllReal(0, :, :))
                r_SpatialPDF_AllReal!0 = "Real"
                r_SpatialPDF_AllReal&Real = s_Real
                ;
                i_dims     := array_append_record(dimsizes(s_Real), \
                                                  dimsizes(r_Var), 0)
                r_VarAllReal := new(i_dims, float)
                copy_VarCoords(r_Var, r_VarAllReal(0, :, :, :))
                copy_VarCoords(r_SpatialPDF_AllReal(:, 0, 0), \
                               r_VarAllReal(:, 0, 0, 0))
            end if
            r_SpatialPDF_AllReal(iReal, :, :) = r_SpatialPDF
            r_VarAllReal(iReal, :, :, :)      = r_Var
        end do 
        ;
        if (iModel .eq. 0 .or. isdefined("r_VarAll") .eq. False) then
            i_dims                    := array_append_record(\
                                dimsizes(s_ModelSelected), dimsizes(r_Var), 0)
            r_VarAll                  := new(i_dims, float)
            copy_VarCoords(r_Var, r_VarAll(0, :, :, :))
            r_VarAll!0                 = "Models"
            r_VarAll&Models            = s_ModelSelected
        end if
        r_VarAll(iModel, :, :, :)      = dim_avg_n_Wrap(r_VarAllReal, 0)
        delete([/r_VarAllReal, r_SpatialPDF_AllReal, r_Var/])
        ;
        print(s_ModelSelected(iModel) + "  Done!")
    end do 
    ;
    r_VarAll           := dim_avg_n_Wrap(r_VarAll, 0)
    ;
    s_DirRead = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/05.GSAT_scaled_RegionalTXx/02.ConstrainAllReal/02.ConstraintByUsingLocalTrend/"
    s_FileRead          = s_DirRead + "FurtherCorrectConstrainedTXxWarming_" + \
                        s_Exp(iExp) + s_4Models + ".nc"
    FileConstrained    := addfile(s_FileRead, "r")
    ;; Constrained MME 
    s_VarName2          = "ConstrainedTXxWarming_" + s_Exp(iExp)
    r_Constrained      := FileConstrained->$s_VarName2$(:, {10:60}, {65:140})
    ;; China Area 
    r_Constrained      := my_SelChina([/r_Constrained/])
    ;
    do iTerm = 0, 2
        l_list         := my_SpatialPDF([/r_VarAll(iTerm, :, :)/], i_Nbin, res)
        r_Spatial      := l_list[0]
        ;
        l_list         := my_SpatialPDF([/r_Constrained(iTerm, :, :)/], i_Nbin,\
                                         res)
        r_SpatialConstrained = l_list[0]
        ;
        if (iTerm .eq. 0) then
            i_dims                := (/3, dimsizes(r_Spatial)/)
            r_SpatialPDF           = new(i_dims, float)
            copy_VarCoords(r_Spatial, r_SpatialPDF(0, :))
            copy_VarCoords(r_VarAll(:, 0, 0), r_SpatialPDF(:, 0))
            r_SpatialPDF_Constrained = r_SpatialPDF
            r_SpatialPDF@long_name = "Spatial PDF of raw MME in " + s_Exp(iExp)
            r_SpatialPDF_Constrained@long_name = \
                            "Spatial PDF of constrained MME in " + s_Exp(iExp)
        end if
        r_SpatialPDF(iTerm, :)     = r_Spatial
        r_SpatialPDF_Constrained(iTerm, :) = r_SpatialConstrained
    end do
    ;
    ;; The Area Fraction where warming exceed local STD (* i_STD_Times)
    r_AreaSig          := r_Area
    r_AreaSig           = where(abs(r_VarAll) .ge. r_STD, r_AreaSig, \
                                r_AreaSig@_FillValue)
    r_AreaSig          := dim_sum_n_Wrap(r_AreaSig, (/1, 2/))
    r_AreaSig           = r_AreaSig / r_AreaTotal * 100.
    r_SpatialPDF@SigArea= r_AreaSig
    r_SpatialPDF@SigAreaUnits = "%"
    r_SpatialPDF@SigDefinition= "one time local STD "
    if (isdefined("i_STD_Times")) then
        r_SpatialPDF@SigDefinition= i_STD_Times + " time local STD "
    end if
    ;
    r_AreaSig_Constrained           := r_Area
    r_AreaSig_Constrained            = where(abs(r_Constrained) .ge. r_STD, \
                        r_AreaSig_Constrained, r_AreaSig_Constrained@_FillValue)
    r_AreaSig_Constrained           := dim_sum_n_Wrap(r_AreaSig_Constrained, \
                                                      (/1, 2/))
    r_AreaSig_Constrained            = r_AreaSig_Constrained / r_AreaTotal*100.
    r_SpatialPDF_Constrained@SigArea = r_AreaSig_Constrained
    r_SpatialPDF_Constrained@SigAreaUnits = "%"
    r_SpatialPDF_Constrained@SigDefinition= "one time local STD "
    if (isdefined("i_STD_Times")) then
        r_SpatialPDF_Constrained@SigDefinition= i_STD_Times + " time local STD "
    end if
    ;
    print(s_Exp(iExp) + "")
    print("Unconstrained & Constrained Area Fraction ")
    print(r_AreaSig + ", " + r_AreaSig_Constrained)
    print("")
    ;
    ;;;;; Output 
    s_VarName2          = "SpatialPDF_" + s_VarName + "_MME_" + s_Exp(iExp)
    FileOut->$s_VarName2$          = r_SpatialPDF
    s_VarName2          = "SpatialPDF_" + s_VarName + "_Constrained_MME_" + \
                            s_Exp(iExp)
    FileOut->$s_VarName2$          = r_SpatialPDF_Constrained
    ;;;;; Output 

    print(s_Exp(iExp) + " Done!")
    print("")
    print("")
end do
;
end





