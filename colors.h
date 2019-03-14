#ifndef PLOT_H
#define PLOT_H

#ifndef __CINT__
#include <iostream>
#include <vector>
#include <cassert>

#include "TColor.h"
#include "TH1.h"
#include "TPad.h"
#include "TList.h"
#include "TClass.h"

//======================================================================
TH1* getPadAxisHist(TPad* pad)
{
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  while (( obj=next() )) {
    // cout << obj->GetName() << endl;
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      // cout << "getPadAxisHist returning hist with name " << obj->GetName() << endl;
      return (TH1*)obj;
    }
  }
  return NULL;
}

//======================================================================
std::vector<TH1*> getPadHists(TPad* pad)
{  
  std::vector<TH1*> ret;
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  while (( obj=next() )) {
    //cout << obj->GetName() << endl;
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      //cout << "getPadAxisHist returning hist with name " << obj->GetName() << endl;
      ret.push_back((TH1*)obj);
    }
  }
  return ret;
}

//======================================================================
double getPadMax(TPad* pad)
{
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  Double_t runningMax=-9e99;//Hparam.ymax;
  while (( obj=next() )) {
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      TH1* curHist=(TH1*)obj;
      const double thisMax=curHist->GetBinContent(curHist->GetMaximumBin());
      if (thisMax > runningMax) {
        runningMax=thisMax;
      }
    }
  }
  return runningMax;
}

//======================================================================
void reMax(TPad* pad, double ymin=0, double headroom=1.1)
{
  TH1* firstHist=getPadAxisHist(pad);
  if(!firstHist){
    std::cerr << "reMax: No histogram in pad " << pad->GetName() << ". Can't reMax" << std::endl;
    return;
  }
  firstHist->GetYaxis()->SetRangeUser(ymin, headroom*getPadMax(pad));
  pad->Update();
}

//======================================================================
enum EColorPalettes
{
  kAlphabetPalette, kKellyPalette, k36Palette, kGlasbeyPalette, kBrewerSet1Palette, kBrewerDark2Palette
};

//======================================================================
const std::vector<int>& getColors(int whichColours)
{
  const int alphabetColors[]={
    TColor::GetColor(  25,  25,  25 ), //ebony
    TColor::GetColor( 240, 163, 255 ), //amethyst
    TColor::GetColor( 255,   0,  16 ), //red
    TColor::GetColor( 128, 128, 128 ), //iron
    TColor::GetColor(  25, 164,   5 ), //orpiment
    TColor::GetColor(   0, 153, 143 ), //turquoise
    TColor::GetColor(   0,  51, 128 ), //navy
    TColor::GetColor(  94, 241, 242 ), //sky
    TColor::GetColor( 157, 204,   0 ), //lime
    TColor::GetColor( 153,  63,   0 ), //caramel
    TColor::GetColor( 224, 255, 102 ), //uranium
    TColor::GetColor(  76,   0,  92 ), //damson
    TColor::GetColor( 255, 168, 187 ), //pink
    TColor::GetColor( 194,   0, 136 ), //mallow
    TColor::GetColor( 255, 204, 153 ), //honeydew
    TColor::GetColor(  43, 206,  72 ), //green
    TColor::GetColor(   0, 117, 220 ), //blue
    TColor::GetColor(   0,  92,  49 ), //forest
    TColor::GetColor( 255, 225,   0 ), //yellow
    TColor::GetColor( 153,   0,   0 ), //wine
    TColor::GetColor( 143, 124,   0 ), //khaki
    TColor::GetColor(  16,  10, 255 ), //violet
    TColor::GetColor( 255, 255, 128 ), //xanthin
    TColor::GetColor( 255,  80,   0 ), //zinnia
    TColor::GetColor( 148, 255, 181 ), //jade
    TColor::GetColor(  66, 102,   0 )  //quagmire
  };
  
  const int kellyColors[]={
    // TColor::GetColor("#f2f3f4"), // white,
    TColor::GetColor("#222222"), // black,
    TColor::GetColor("#f3c300"), // yellow,
    TColor::GetColor("#875692"), // purple,
    TColor::GetColor("#f38400"), // orange,
    TColor::GetColor("#a1caf1"), // lightblue,
    TColor::GetColor("#be0032"), // red,
    TColor::GetColor("#c2b280"), // buff,
    TColor::GetColor("#848482"), // gray,
    TColor::GetColor("#008856"), // green,
    TColor::GetColor("#e68fac"), // purplishpink,
    TColor::GetColor("#0067a5"), // blue,
    TColor::GetColor("#f99379"), // yellowishpink,
    TColor::GetColor("#604e97"), // violet,
    TColor::GetColor("#f6a600"), // orangeyellow,
    TColor::GetColor("#b3446c"), // purplishred,
    TColor::GetColor("#dcd300"), // greenishyellow,
    TColor::GetColor("#882d17"), // reddishbrown,
    TColor::GetColor("#8db600"), // yellowgreen,
    TColor::GetColor("#654522"), // yellowishbrown,
    TColor::GetColor("#e25822"), // reddishorange,
    TColor::GetColor("#2b3d26") // olivegreen
  };

  const int palette36Colors[]={
    TColor::GetColor("#5A5156"),
    // TColor::GetColor("#E4E1E3"), // much too light
    TColor::GetColor("#F6222E"),
    TColor::GetColor("#FE00FA"),
    TColor::GetColor("#16FF32"),
    TColor::GetColor("#3283FE"),
    TColor::GetColor("#FEAF16"),
    TColor::GetColor("#B00068"),
    TColor::GetColor("#1CFFCE"),
    TColor::GetColor("#90AD1C"),
    TColor::GetColor("#2ED9FF"),
    TColor::GetColor("#DEA0FD"),
    TColor::GetColor("#AA0DFE"),
    TColor::GetColor("#F8A19F"),
    TColor::GetColor("#325A9B"),
    TColor::GetColor("#C4451C"),
    TColor::GetColor("#1C8356"),
    TColor::GetColor("#85660D"),
    TColor::GetColor("#B10DA1"),
    TColor::GetColor("#FBE426"),
    TColor::GetColor("#1CBE4F"),
    TColor::GetColor("#FA0087"),
    TColor::GetColor("#FC1CBF"),
    TColor::GetColor("#F7E1A0"),
    TColor::GetColor("#C075A6"),
    TColor::GetColor("#782AB6"),
    TColor::GetColor("#AAF400"),
    TColor::GetColor("#BDCDFF"),
    TColor::GetColor("#822E1C"),
    TColor::GetColor("#B5EFB5"),
    TColor::GetColor("#7ED7D1"),
    TColor::GetColor("#1C7F93"),
    TColor::GetColor("#D85FF7"),
    TColor::GetColor("#683B79"),
    TColor::GetColor("#66B0FF"),
    TColor::GetColor("#3B00FB"),
  };

  const int glasbeyColors[]={
    TColor::GetColor("#0000FF"),
    TColor::GetColor("#FF0000"),
    TColor::GetColor("#00FF00"),
    TColor::GetColor("#000033"),
    TColor::GetColor("#FF00B6"),
    TColor::GetColor("#005300"),
    TColor::GetColor("#FFD300"),
    TColor::GetColor("#009FFF"),
    TColor::GetColor("#9A4D42"),
    TColor::GetColor("#00FFBE"),
    TColor::GetColor("#783FC1"),
    TColor::GetColor("#1F9698"),
    TColor::GetColor("#FFACFD"),
    TColor::GetColor("#B1CC71"),
    TColor::GetColor("#F1085C"),
    TColor::GetColor("#FE8F42"),
    TColor::GetColor("#DD00FF"),
    TColor::GetColor("#201A01"),
    TColor::GetColor("#720055"),
    TColor::GetColor("#766C95"),
    TColor::GetColor("#02AD24"),
    TColor::GetColor("#C8FF00"),
    TColor::GetColor("#886C00"),
    TColor::GetColor("#FFB79F"),
    TColor::GetColor("#858567"),
    TColor::GetColor("#A10300"),
    TColor::GetColor("#14F9FF"),
    TColor::GetColor("#00479E"),
    TColor::GetColor("#DC5E93"),
    TColor::GetColor("#93D4FF"),
    TColor::GetColor("#004CFF")
  };

  const int brewerSet1Colors[]={
    TColor::GetColor("#e41a1c"),
    TColor::GetColor("#377eb8"),
    TColor::GetColor("#4daf4a"),
    TColor::GetColor("#984ea3"),
    TColor::GetColor("#ff7f00"),
    TColor::GetColor("#ffff33"),
    TColor::GetColor("#a65628"),
    TColor::GetColor("#f781bf")
  };

  const int brewerDark2Colors[]={
    TColor::GetColor( 27, 158, 119),
    TColor::GetColor(217,  95,   2),
    TColor::GetColor(117, 112, 179),
    TColor::GetColor(231,  41, 138),
    TColor::GetColor(102, 166,  30),
    TColor::GetColor(230, 171,   2),
    TColor::GetColor(166, 118,  29),
    TColor::GetColor(102, 102, 102)
  };

  const int nPalettes=6;
  static std::vector<std::vector<int> > allPalettes(nPalettes);
  static bool firstCall=true;

  if(firstCall){
    for(int j=0; j<nPalettes; ++j){
      int nColors=0;
      const int* firstColor=NULL;

      switch(j){
      case 0:
        nColors=26;
        firstColor=alphabetColors;
        break;
      case 1:
        nColors=21;
        firstColor=kellyColors;
        break;
      case 2:
        nColors=35; // 36 because I commented one out
        firstColor=palette36Colors;
        break;
      case 3:
        nColors=30;
        firstColor=glasbeyColors;
        break;
      case 4:
        nColors=8;
        firstColor=brewerSet1Colors;
        break;
      case 5:
        nColors=8;
        firstColor=brewerDark2Colors;
        break;
      default:
        assert(0 && "whichColours must be 0-5");
      }

      for(int i=0; i<nColors; ++i){
        allPalettes[j].push_back(firstColor[i]);
      }
    }
    firstCall=false;
  }
  return allPalettes[whichColours];
}

//======================================================================
void autoColorHists(TPad* pad, int whichColours=kBrewerSet1Palette)
{
  const std::vector<int>& colours=getColors(whichColours);

  std::vector<TH1*> hists=getPadHists(pad);
  for(unsigned int i=0; i<hists.size(); ++i){
    hists[i]->SetLineColor(colours[i%colours.size()]);
  }
  pad->Draw();
}

#endif // #ifndef __CINT__
#endif // include guard
