; =======================================================
; Author: Ziming Chen
; This code is for: 
;   - sel the data of variables over China 
;   - cal the clim, interannual var, and trend in given period 
;   1. indicate models, period, var, exp
;   2. for each realization of all models, read the var at first, 
; and then take the data over China in given period 
;   3. cal the clim, interannual var, and trend 
;   4. output 
; =======================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code/my_FunctionInPNNL.ncl"
;
undef("my_VarClimSTD_Trend")
function my_VarClimSTD_Trend(Var[*][*][*]:numeric)
local r_Var 
begin
;
r_Var               = Var 
i_dims              = array_append_record(3, dimsizes(r_Var(0, :, :)), 0)
r_Return            = new(i_dims, float)
copy_VarCoords(r_Var(0, :, :), r_Return(0, :, :))
r_Return!0          = "Clim_STD_Trend"
r_Return&Clim_STD_Trend = (/0, 1, 2/)
r_Return@units      = r_Var@units + ", " + r_Var@units + ", " + r_Var@units + \
                        "/" + dimsizes(r_Var(:, 0, 0)) + "yr"
;
;; Clim 
r_Return(0, :, :)   = dim_avg_n_Wrap(r_Var, 0)
;
;; STD 
r_Return(1, :, :)   = dim_stddev_n_Wrap(r_Var, 0)
;
;; Trend 
r_Temp             := trend_manken(r_Var, False, 0)
r_Return(2, :, :)   = r_Temp(1, :, :) * dimsizes(r_Var(:, 0, 0))
;
r_Return(2, :, :)   = where(ismissing(r_Return(0, :, :)), r_Return@_FillValue, \
                            r_Return(2, :, :))
;
return(r_Return)
;
end
;
begin
; =======================================================
;   1. indicate models, period, var, exp
s_ModelSelect           = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                            "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                            "CNRM-ESM2-1", "CanESM5", "CanESM5-CanOE", \
                            "EC-Earth3", "EC-Earth3-Veg", "FGOALS-f3-L", \
                            "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
                            "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8", \
                            "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR", \
                            "MRI-ESM2-0", "UKESM1-0-LL"/)
; s_ModelSelect        = (/"GFDL-CM4", "CNRM-CM6-1", "EC-Earth3"/)
s_4Models            = "_23Models"
i_yrRead             = (/1981, 2014/)
s_Exp                = "historical"
; s_VarName            = (/"pr", "rlus", "rsds", "huss", "hfls", "ts"/)
s_VarName            = (/"pr", "rlus", "rsds", "huss", "psl", "hfls", "ts", \
                         "tasmax"/)
;
; =======================================================
;   2. for each realization of all models, read the var at first, 
; and then take the data over China in given period 
do iModel = 0, dimsizes(s_ModelSelect)-1
    l_list         := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp, "Amon")
    s_ExpModel_Path:= l_list[0]
    s_Real         := l_list[1]
    l_list_        := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp, "Amon")
    s_ExpModel_Path_:= l_list_[0]
    s_Real_         := l_list_[1]
    if (all(ismissing(s_Real)) .or. dimsizes(s_Real_).gt.dimsizes(s_Real)) then
        s_Real      := s_Real_
        s_ExpModel_Path := s_ExpModel_Path_
    end if
        ;
    print("Models No.: " + iModel + ", " + s_ModelSelect(iModel) + \
              "; Exp: " + s_Exp)
    ;
    do ivar = 0, dimsizes(s_VarName)-1
        print(s_VarName(ivar) + "")
        s_DirRead   = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/" + \
                      s_VarName(ivar) + "/2.5x2.5/AllReal/"
        s_FileRead := s_DirRead + s_VarName(ivar) + "_" + s_Exp + "_" + \
                            s_ModelSelect(iModel) + "_CMIP6_AllReal.nc"
        ls         := systemfunc("ls " + s_FileRead)
        if (ismissing(ls(0))) then
            print(s_FileRead + "  Missing!")
            continue
        end if
        FileHist   := addfile(s_FileRead, "r")
        ;
        do iReal = 0, dimsizes(s_Real)-1
            s_VarName2 := s_VarName(ivar) + "_" + s_ModelSelect(iModel) + "_" +\
                            s_Exp + "_" + s_Real(iReal)
            if (my_VarNameExist(FileHist, s_VarName2) .eq. False) then
                print(s_VarName2 + "  Missing!")
                continue
            end if
            r_Hist     := FileHist->$s_VarName2$({i_yrRead(0):i_yrRead(1)}, :, \
                                                 {-90:90}, {0:360})
            r_Hist     := dim_avg_n_Wrap(r_Hist, 1)
            ;
            r_Hist     := my_SelChina([/r_Hist/])
            ; FileTest = "./Test.nc"
            ; system("rm -rf " + FileTest)
            ; FileTest:= addfile(FileTest, "c")
            ; FileTest->Hist = r_Hist 
            ; exit()
            ; =======================================================
            ;   3. cal the clim, interannual var, and trend 
            r_Clim_STD_Trend := my_VarClimSTD_Trend(r_Hist)
            ;
            ; =======================================================
            ;   4. output 
            s_DirOut    = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/08.RelatedToOtherVar/"
            system("mkdir -p " + s_DirOut)
            s_FileOut  := s_DirOut + s_VarName(ivar) + "_Clim_STD_Trend_" + \
                          s_Exp + "_" + i_yrRead(0) + "-" + i_yrRead(1) + \
                          s_4Models + ".nc"
            ls         := systemfunc("ls " + s_FileOut)
            if (ismissing(ls(0)) ) then
                FileOut:= addfile(s_FileOut, "c")
                GAttr  := True
                GAttr@history = systemfunc("date") + "  " + systemfunc("pwd") + \
                                "/01.VarClim_STD_TrendOverChina.ncl"
                GAttr@Exp     = s_Exp
                GAttr@Period  = i_yrRead(0) + " to " + i_yrRead(1)
                fileattdef(FileOut, GAttr)
            else 
                FileOut:= addfile(s_FileOut, "w")
            end if
            s_VarName2 := s_VarName(ivar) + "_" + s_ModelSelect(iModel) + "_" + \
                            s_Exp + "_" + s_Real(iReal)
            FileOut->$s_VarName2$ = r_Clim_STD_Trend
            print(s_Real(iReal) + " Done!")
            ;
        end do 
        print("")
    end do
    print("")
end do
;
end












