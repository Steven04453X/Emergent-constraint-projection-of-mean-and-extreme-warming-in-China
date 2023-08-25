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
; =====================================================================
;   1. indicate the model, exp, periods and their name 
; s_ModelSelect           = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
;                             "BCC-CSM2-MR", "CAMS-CSM1-0", "CESM2-WACCM", \
;                             "CNRM-CM6-1", "CNRM-ESM2-1", "CanESM5", \
;                             "CanESM5-CanOE", "EC-Earth3", "EC-Earth3-Veg", \
;                             "FGOALS-f3-L", "FGOALS-g3", "GFDL-CM4", \
;                             "GFDL-ESM4", "GISS-E2-1-G", "HadGEM3-GC31-LL", \
;                             "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", \
;                             "MCM-UA-1-0", "MIROC6", "MIROC-ES2L", \
;                             "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", \
;                             "NESM3", "NorESM2-LM", "UKESM1-0-LL"/)
; Rerun:  
s_ModelSelect           = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                            "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                            "CNRM-ESM2-1", "CanESM5", "CanESM5-CanOE",\
                            "EC-Earth3", "EC-Earth3-Veg", "FGOALS-f3-L", \
                            "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
                            "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8",\
                            "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR", \
                            "MRI-ESM2-0", "UKESM1-0-LL"/)

s_4Models               = "_23Models"
;
s_Exp                   = (/"ssp126", "ssp245", "ssp370", "ssp585"/)
i_yrBaseState           = (/1995, 2014/)
;
s_VarName               = "tas"
;
; =====================================================================
s_DirOut = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/01.ChinaMean/"
; s_DirOut  = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/ProjectionOfAllReal/"
system("mkdir -p " + s_DirOut)
s_FileOut               = s_DirOut + s_VarName + "_ProjectionInChina_AnMean.nc"
ls                     := systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileOut_ProjChanges = addfile(s_FileOut, "c")
    GAttr               = True
    GAttr@history       = systemfunc("date") + "  " + systemfunc("pwd") + \
                        "/01.AreaMeanWarmingInChina.ncl"
    fileattdef(FileOut_ProjChanges, GAttr)
else 
    FileOut_ProjChanges = addfile(s_FileOut, "w")
end if
;
s_FileOut               = s_DirOut + s_VarName + "_HistoricalInChina_AnMean.nc"
ls                     := systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileOut_Hist        = addfile(s_FileOut, "c")
    GAttr               = True
    GAttr@history       = systemfunc("date") + "  " + systemfunc("pwd") + \
                        "/01.AreaMeanWarmingInChina.ncl"
    fileattdef(FileOut_Hist, GAttr)
else 
    FileOut_Hist        = addfile(s_FileOut, "w")
end if
;
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/2.5x2.5/AllReal/"
do iExp = 0, dimsizes(s_Exp)-1
    do iModel = 0, dimsizes(s_ModelSelect)-1
        l_list         := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "Amon")
        ; s_ExpModel_Path:= l_list[0]
        s_Real         := l_list[1]
        l_list_        := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "Amon")
        ; s_ExpModel_Path_:= l_list_[0]
        s_Real_         := l_list_[1]
        if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real) ) then
            s_Real      := s_Real_
            ; s_ExpModel_Path := s_ExpModel_Path_
        end if
        ;
        print("Models No.: " + iModel + ", " + s_ModelSelect(iModel) + \
              "; Exp: " + s_Exp(iExp))
        ;
        s_FileRead     := s_DirRead + s_VarName + "_" + s_Exp(iExp) + "_" + \
                            s_ModelSelect(iModel) + "_CMIP6_AllReal.nc"
        ls             := systemfunc("ls " + s_FileRead)
        if (ismissing(ls(0))) then
            print(s_FileRead + "  Missing!")
            continue
        end if
        FileProj       := addfile(s_FileRead, "r")
        s_FileRead     := s_DirRead + s_VarName + "_historical_" + \
                        s_ModelSelect(iModel) + "_CMIP6_AllReal.nc"
        FileHist       := addfile(s_FileRead, "r")
        ;
        do iReal = 0, dimsizes(s_Real)-1
            s_VarName2 := s_VarName + "_" + s_ModelSelect(iModel) + "_" + \
                            s_Exp(iExp) + "_" + s_Real(iReal)
            if (my_VarNameExist(FileProj, s_VarName2) .eq. False) then
                print(s_VarName2 + "  Missing!")
                continue
            end if
            r_Proj     := FileProj->$s_VarName2$(:, :, {-90:90}, {0:360})
            r_Proj     := dim_avg_n_Wrap(r_Proj, 1)
            ;
            s_VarName2 := s_VarName + "_" + s_ModelSelect(iModel) + "_" + \
                            "historical_" + s_Real(iReal)
            if (my_VarNameExist(FileHist, s_VarName2) .eq. False) then
                print(s_VarName2 + "  Missing!")
                continue
            end if
            r_Hist     := FileHist->$s_VarName2$(:, :, {-90:90}, {0:360})
            r_Hist     := dim_avg_n_Wrap(r_Hist, 1)
            ;
            r_BasicState  := dim_avg_n_Wrap(r_Hist(\
                            {i_yrBaseState(0):i_yrBaseState(1)}, :, :), 0)
            r_BasicState0 := conform(r_Hist, r_BasicState, (/1, 2/))
            copy_VarCoords(r_Hist, r_BasicState0)
            r_Hist         = r_Hist - r_BasicState0
            ;
            r_BasicState0 := conform(r_Proj, r_BasicState, (/1, 2/))
            copy_VarCoords(r_Proj, r_BasicState0)
            r_Proj         = r_Proj - r_BasicState0
            ;
            ;; Area-Mean 
            r_Hist        := my_SelChina([/r_Hist/])
            r_Proj        := my_SelChina([/r_Proj/])
            r_Hist        := my_AreaMean2([/r_Hist/], r_Hist&lat, r_Hist&lon, 0)
            r_Proj        := my_AreaMean2([/r_Proj/], r_Proj&lat, r_Proj&lon, 0)
            ;
            r_BasicState  := my_SelChina([/r_BasicState/])
            r_BasicState  := my_AreaMean2([/r_BasicState/], r_BasicState&lat, \
                                          r_BasicState&lon, 0)
            r_Hist@BasicState = r_BasicState
            r_Proj@BasicState = r_BasicState
            ;
            ;; Output 
            s_VarName2    := s_VarName + "_ChinaArea_" + s_ModelSelect(iModel) \
                            + "_" + s_Exp(iExp) + "_" + s_Real(iReal)
            r_Proj@long_name = s_VarName2
            FileOut_ProjChanges->$s_VarName2$ = r_Proj
            ;
            s_VarName2    := s_VarName + "_ChinaArea_" + s_ModelSelect(iModel) \
                            + "_historical_" +s_Real(iReal)
            r_Hist@long_name = s_VarName2
            FileOut_Hist->$s_VarName2$ = r_Hist
        end do 
        print(s_ModelSelect(iModel) + ", " + s_Exp(iExp) + " Done!")
    end do 
    print("")
end do 
;
end 

















