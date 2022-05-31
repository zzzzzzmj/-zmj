begin
filename1 = "/mnt/d/zmjuse/fortran/uc2.dat"
filename2 = "/mnt/d/zmjuse/fortran/vc2.dat"
filename3 = "/mnt/d/zmjuse/fortran/zc2.dat"

u = asciiread(filename1,(/16,20/), "float")
v = asciiread(filename2,(/16,20/), "float")
z = asciiread(filename3,(/16,20/), "float")
u!0 = "lat"
u!1 = "lon"
u&lat = fspan(28.5,73.5,16)
u&lon = fspan(89.5,146.5,20)
u&lat@units = "degree_N"
u&lon@units = "degree_E"
copy_VarCoords(u,v)
copy_VarCoords(u,z)

wks = gsn_open_wks("pdf", "/mnt/d/zmjuse/fortran/24_2")
gsn_define_colormap(wks, "BlueWhiteOrangeRed")

res = True


res@mpFillOn           = False
res@mpGridAndLimbOn    = True    ;设置经纬线
res@mpGridLineTicknessF = 0.5
res@mpGridLineDashPattern  = 2   ;等值线样式
res@GridSpacingF           = 10
res@mpMinLonF          = 89.5    ;地图区域 
res@mpMaxLonF          = 146.5
res@mpMinLatF          = 28.5    
res@mpMaxLatF          = 73.5


res@gsnDraw = False
res@gsnFrame = False
res@gsnAddCyclic = False
res@gsnCenterString = "wind-ght-24hours"

res@gsnLeftString    =""
res@gsnRightString    =""
res@pmTickMarkDisplayMode  = "Always"

resv   =  res
res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn     = True
res@cnLineLabelBackgroundColor = -1
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 5000
res@cnMaxLevelValF = 6000
res@cnLevelSpacingF = 50



resv@vcGlyphStyle = "LineArrow"
resv@vcMinDistanceF            = 0.03
resv@vcLineArrowHeadMaxSizeF = 0.005
resv@vcLineArrowHeadMinSizeF = 0.005

resv@vcRefAnnoOn              = True
resv@vcRefLengthF             = 0.03
resv@vcRefAnnoBackgroundColor = "white"
resv@vcRefAnnoPerimOn         = False
resv@vcRefAnnoFontHeightF     = 0.015
resv@vcRefAnnoString1On       = False
resv@vcRefAnnoString2On       = True
resv@vcRefMagnitudeF   = 20
resv@vcRefAnnoString2  = "20m/s"
resv@vcRefAnnoSide            = "top"
resv@vcRefAnnoOrthogonalPosF  = -0.05
resv@vcRefAnnoParallelPosF    = 1

plot_cn = gsn_csm_contour_map(wks, z, res)
plot_vc = gsn_csm_vector(wks, u, v, resv)



overlay(plot_cn,plot_vc)     ;合并

draw(plot_cn)
frame(wks)

end
