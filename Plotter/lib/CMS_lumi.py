import ROOT as rt

# CMS_lumi
#   Initiated by: Gautier Hamel de Monchenault (Saclay)
#   Translated in Python by: Joshua Hardenbrook (Princeton)
#   Updated by:   Dinko Ferencek (Rutgers)
#

cmsText     = "CMS";
cmsTextFont   = 61

writeExtraText = True
extraText   = "Preliminary"
#extraText   = "Private"
extraTextFont = 52

lumiTextSize     = 0.85
lumiTextOffset   = 0.2
lumiText_posX = 0.07
lumiText_posY = 0.015

cmsTextSize      = 0.75
cmsTextOffset    = 0.1
CMSText_posX = 0.1

relPosX    = 0.045
relPosY    = 0.035
relExtraDY = 1.2

extraOverCmsTextSize  = 0.76
extraText_posX = 0.07

lumi_14TeV = "3000 fb^{-1}"
lumi_13TeV = "35.9 fb^{-1}"
lumi_8TeV  = "19.7 fb^{-1}"
lumi_7TeV  = "5.1 fb^{-1}"
#lumi_sqrtS = "35.9 fb^{-1}"

drawLogo      = False
outOfFrame    = False

def CMS_lumi(pad,  iPeriod,  iPosX, year ):
    if year == '2016':
        lumi_sqrtS = "16.81 fb^{-1}"
    elif year == '-2016':
        lumi_sqrtS = "19.52 fb^{-1}"
    elif year == '2017':
        lumi_sqrtS = "41.48 fb^{-1}"
    elif year == '2018':
        lumi_sqrtS = "59.83 fb^{-1}"
    elif year == 'run2':
        lumi_sqrtS = "138 fb^{-1}"
    else:
        print("do not include at 2016/2017/2018")
        exit(0)
    global outOfFrame, relPosX
    if(iPosX/10==0 ): outOfFrame = True

    alignY_=3
    alignX_=2
    if( iPosX/10==0 ): alignX_=1
    if( iPosX==0    ): alignY_=1
    if( iPosX/10==1 ): alignX_=1
    if( iPosX/10==2 ): alignX_=2
    if( iPosX/10==3 ): alignX_=3
    align_ = 10*alignX_ + alignY_

    H = pad.GetWh()
    W = pad.GetWw()
    l = pad.GetLeftMargin()-0.03
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()
    e = 0.025

    pad.cd()

    lumiText = ""
    if( iPeriod==1 ):
        lumiText += lumi_7TeV
        lumiText += " (7 TeV)"
    elif ( iPeriod==2 ):
        lumiText += lumi_8TeV
        lumiText += " (8 TeV)"

    elif( iPeriod==3 ):
        lumiText = lumi_8TeV
        lumiText += " (8 TeV)"
        lumiText += " + "
        lumiText += lumi_7TeV
        lumiText += " (7 TeV)"
    elif ( iPeriod==4 ):
        # lumiText += lumi_13TeV
        lumiText += lumi_sqrtS
        lumiText += " (13 TeV)"
    elif ( iPeriod==7 ):
        if( outOfFrame ):lumiText += "#scale[0.85]{"
        lumiText += lumi_13TeV
        lumiText += " (13 TeV)"
        lumiText += " + "
        lumiText += lumi_8TeV
        lumiText += " (8 TeV)"
        lumiText += " + "
        lumiText += lumi_7TeV
        lumiText += " (7 TeV)"
        if( outOfFrame): lumiText += "}"
    elif ( iPeriod==12 ):
        lumiText += "8 TeV"
    elif ( iPeriod==13 ):
        if( outOfFrame ):lumiText += "#scale[0.90]{"
        lumiText += lumi_13TeV
        lumiText += " (13 TeV)"
        if( outOfFrame): lumiText += "}"
    elif ( iPeriod==14 ):
        if( outOfFrame ):lumiText += "#scale[0.90]{"
        lumiText += lumi_14TeV
        lumiText += " (14 TeV, 200 PU)"
        if( outOfFrame): lumiText += "}"
    elif ( iPeriod==0 ):
        lumiText += lumi_sqrtS
    elif ( iPeriod==5 ):
        lumiText += "13 TeV"

    print(lumiText)

    latex = rt.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(rt.kBlack)

    extraTextSize = extraOverCmsTextSize*cmsTextSize

    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(lumiTextSize*t)

    latex.DrawLatex(1-r-lumiText_posX,1-t+lumiTextOffset*t-lumiText_posY,lumiText)

    if( outOfFrame ):
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11)
        latex.SetTextSize(cmsTextSize*t)
        latex.DrawLatex(l-CMSText_posX,1-t+lumiTextOffset*t-lumiText_posY,cmsText)

    pad.cd()

    posX_ = 0
    if( iPosX%10<=1 ):
        posX_ =   l + relPosX*(1-l-r)   # left aligned
    elif( iPosX%10==2 ):
        posX_ =  l + 0.5*(1-l-r)        # centered
    elif( iPosX%10==3 ):
        posX_ =  1-r - relPosX*(1-l-r)  # right aligned

    posY_ = 1-t - relPosY*(1-t-b)

    if( not outOfFrame ):
        if( drawLogo ):
            posX_ =   l + 0.045*(1-l-r)*W/H
            posY_ = 1-t - 0.045*(1-t-b)
            xl_0 = posX_
            yl_0 = posY_ - 0.15
            xl_1 = posX_ + 0.15*H/W
            yl_1 = posY_
            CMS_logo = rt.TASImage("CMS-BW-label.png")
            pad_logo =  rt.TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 )
            pad_logo.Draw()
            pad_logo.cd()
            CMS_logo.Draw("X")
            pad_logo.Modified()
            pad.cd()
        else:
            latex.SetTextFont(cmsTextFont)
            latex.SetTextSize(cmsTextSize*t)
            latex.SetTextAlign(align_)
            latex.DrawLatex(posX_, posY_, cmsText)
            if( writeExtraText ) :
                latex.SetTextFont(extraTextFont)
                latex.SetTextAlign(align_)
                latex.SetTextSize(extraTextSize*t)
                latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText)
    elif( writeExtraText ):
        if( iPosX==0):
            posX_ =   l +  relPosX*(1-l-r)+extraText_posX
            posY_ =   1-t+lumiTextOffset*t-lumiText_posY

        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize*t)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_-CMSText_posX, posY_, extraText)

    pad.Update()
