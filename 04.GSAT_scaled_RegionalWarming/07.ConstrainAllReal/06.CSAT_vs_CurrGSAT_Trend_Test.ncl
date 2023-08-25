; =========================================================
; Author: Ziming Chen 
; This code is for:
;   - vis the scatter plot between projected CSAT warming and present-day 
; GSAT trend 
;   1. indicate the models, exp, 
;   2. read the GSAT trend and CSAT warming 
;   3. cal the inter-model correlation 
;   4. vis the scatter plot 
; =========================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code/my_FunctionInPNNL.ncl"
;
begin
;
; =========================================================
;   1. indicate the models, exp, 
s_ModelSelected           = (/ "ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                         "BCC-CSM2-MR", "CAMS-CSM1-0", "CanESM5", \
                         "CanESM5-CanOE", "CNRM-CM6-1", "CNRM-ESM2-1", \
                         "EC-Earth3", "EC-Earth3-Veg", "FGOALS-f3-L", \
                         "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
                         "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8", \
                         "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR", \
                         "MRI-ESM2-0", "UKESM1-0-LL"/)
s_4Models                 = "_23Models"
s_Exp                     = "ssp585"
s_ExpTitle                = "SSP5-8.5"

s_DirOut = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Fig/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/"
system("mkdir -p " + s_DirOut)
s_FigOut                  = s_DirOut + "06.CSAT_vs_CurrGSAT_Trend_Test_" + s_Exp
wks                       = gsn_open_wks("X11", s_FigOut)
;
; =========================================================
;   2. read the GSAT trend and CSAT warming 
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/GSAT_AndRegionalTrend_AllReal/"
s_FileRead                = s_DirRead + "GSAT_AndRegionalTrend_1981-2014"\
                                + s_4Models + ".nc"
FileGSAT_PD               = addfile(s_FileRead, "r")
;
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/ProjectionOfAllReal/"
s_FileRead                = s_DirRead + "ProjectionInDifTerms_AllReal" + \
                                s_4Models + ".nc"
FileGSAT_Proj             = addfile(s_FileRead, "r")
;
r_GSAT_Trend              = new(dimsizes(s_ModelSelected), float)
r_ProjWarming             = new((/dimsizes(s_ModelSelected), 3/), float)
r_ProjWarming!0           = "Models"
r_ProjWarming&Models      = s_ModelSelected
r_ProjWarming!1           = "Terms"
r_ProjWarming&Terms       = (/ "Near", "Mid", "Long" /)
;
do iModel = 0, dimsizes(s_ModelSelected)-1
    l_list               := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp, "Amon")
    s_Real               := l_list[1]
    l_list_              := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp, "Amon")
    s_Real_              := l_list_[1]
    if (all(ismissing((s_Real))) .or. dimsizes(s_Real_) .gt. dimsizes(s_Real)) then
        s_Real           := s_Real_
    end if
    ;
    if (all(ismissing(s_Real))) then
        print(s_ModelSelected(iModel) + " Missing!")
        continue
    end if
    r_GSAT_Trend         := new(dimsizes(s_Real), float)
    r_CSAT               := new((/dimsizes(s_Real), 3/), float)
    do iReal = 0, dimsizes(s_Real)-1
        s_VarName2        = "GSAT_Trend_" + s_ModelSelected(iModel) + \
                            "_historical_" + s_Real(iReal)
        if (my_VarNameExist(FileGSAT_PD, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue
        end if  
        r_GSAT_Trend(iReal) = FileGSAT_PD->$s_VarName2$(1)
        ;
        s_VarName2        = "tas_inDifTerms_" + s_ModelSelected(iModel) + \
                            "_" + s_Exp + "_" + s_Real(iReal)
        if (my_VarNameExist(FileGSAT_Proj, s_VarName2) .eq. False) then
            print(s_VarName2 + "  Missing!")
            continue  
        end if
        r_Temp           := FileGSAT_Proj->$s_VarName2$
        ;
        r_Temp           := my_SelChina([/r_Temp/])
        r_Temp           := my_AreaMean2([/r_Temp/], r_Temp&lat, r_Temp&lon, 0)
        ;
        r_CSAT(iReal, :)  = r_Temp
    end do
    ;
    if (iModel .eq. 0 .or. isdefined("r_GSAT_TrendAll") .eq. False) then
        i_dims           := (/dimsizes(s_ModelSelected), 3/)
        r_CSAT_All       := new(i_dims, float)
        r_CSAT_All!0      = "Models"
        r_CSAT_All&Models = s_ModelSelected
        r_CSAT_All!1      = "Terms"
        r_CSAT_All&Terms  = (/0, 1, 2/)
        ;
        r_GSAT_TrendAll  := r_CSAT_All(:, 0)
    end if    
    r_GSAT_TrendAll(iModel) = avg(r_GSAT_Trend)
    r_CSAT_All(iModel, :) = dim_avg_n_Wrap(r_CSAT, 0)
end do
; =========================================================
;   3. cal the inter-model correlation 
r_CorCoef                 = escorc_n(r_GSAT_TrendAll, r_CSAT_All, 0, 0)
print("Cor coef in three terms")
print(r_CorCoef + "")
print("")
print(s_ModelSelected + ": " + sprintf("%4.2f", r_GSAT_TrendAll) + ", " + \
      sprintf("%4.2f", r_CSAT_All(:, 2)))
;
; =========================================================
;   4. vis the scatter plot 
res                       = True
res@trXMinF               = min(r_GSAT_TrendAll) - 0.1
res@trXMaxF               = max(r_GSAT_TrendAll) + 0.1
res@trYMinF               = min(r_CSAT_All) - 0.2
res@trYMaxF               = max(r_CSAT_All) + 0.2
;
if (s_Exp .eq. "ssp126") then
    s_LeftString             := (/ "(a1) ", "(b1) ", "(c1) " /)
else if (s_Exp .eq. "ssp245") then
    s_LeftString             := (/ "(a2) ", "(b2) ", "(c2) " /)
else if (s_Exp .eq. "ssp370") then
    s_LeftString             := (/ "(a3) ", "(b3) ", "(c3) " /)
else if (s_Exp .eq. "ssp585") then
    s_LeftString             := (/ "(a4) ", "(b4) ", "(c4) " /)
end if
end if
end if
end if
;
do iSub = 0, 2
    res@tiMainString      = s_ExpTitle
    if (iSub .gt. 0) then
        res@tiMainString  = ""
    end if
    res@gsnLeftString     = s_LeftString(iSub) 
end do
;
end










