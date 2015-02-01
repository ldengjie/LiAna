#include  <iostream>
//#include <stdio.h>
//#include  <stdlib.h>
#include  <dirent.h>//DIR
#include  <TCanvas.h>
#include  <TMath.h>
#include  <TLegend.h>
#include  <TH1F.h>
#include  <TFile.h>
#include  <fstream>
#include  <sstream>
#include  <TStyle.h>
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include  "RooSimultaneous.h"
#include  "RooAbsReal.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooChi2Var.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooMinuit.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include  "RooCategory.h"
#include  <string>
#include  <TChain.h>
#include  <TROOT.h>
#include  <vector>
#include "RooHistPdf.h"
#include "RooExtendPdf.h"
#include "RooDataHist.h"
#include "RooChebychev.h"
#include "RooDataSet.h"
#include  "RooAddition.h"
using namespace std;
using namespace RooFit;

TString nameStr,nameStr2;
RooRealVar* xe;
RooRealVar* xt;
RooRealVar* rateMu;

typedef struct isoItems
{
    string isoName;
    bool isFixed;
    double tau_val;
    double tau_val_low;
    double tau_val_high;
    double N_value;
    double N_val_low;
    double N_val_high;
    Color_t linecolor;
    double tauInSlice[3];
    double tauErrInSlice[3];
    double tNumInSlice[3];
    double tNumErrInSlice[3];
    double eNumInSlice[3];
    double eNumErrInSlice[3];
    double eNumErr;
    double realNumInSlice[3];
    double realNumErrInSlice[3];
    double totalRealNum;
    double totalRealNumErr2;
    double rate;
    double rateErr;
    RooRealVar* fitTau;
    RooRealVar* tFitNum;
    RooFormulaVar* fitLambda;
    RooExponential* fitExpPdf;
    RooRealVar* timeTCutCoe;
    RooRealVar* timeECutCoe;
    RooExtendPdf* fitExtendPdf;
    RooDataHist* ehd;//RooHistPdf neither owns or clone 'dhist' and the user must ensure the input histogram exists for the entire life span of this PDF.http://root.cern.ch/root/html/RooHistPdf.html#RooHistPdf:RooHistPdf@1
    RooHistPdf* fitHistPdf;
    RooRealVar* specECutCoe;
    RooRealVar* specTCutCoe;
    //RooRealVar* specECutCoe;
    RooFormulaVar* eFitNum;
    /*
       isoItems()
       {
       for( int i=0 ; i<5 ; i++ )
       {
       tauInSlice[i]=0.;
       tauErrInSlice[i]=0.;
       tNumInSlice[i]=0.;
       tNumErrInSlice[i]=0.;
       }

       }
       */
}isoItem;
typedef struct FitInf
{
    string mode;
    bool isNoRed;
    bool isbinned;
    bool doSimulFit;
    double timeElow;
    double timeEhigh;
    double timeTlow;
    double timeThigh;
    double specElow;
    double specEhigh;
    double specTlow;
    double specThigh;
    string timeCom[10];//component
    string specCom[10];//component
    string con;//contour
    double muonRate[4];
    double daqTime;
    double liveTime;
    map<string,isoItem*> timeComMap;
    map<string,isoItem*> specComMap;
    RooAddPdf* timeFitPdf;
    RooAddPdf* specFitPdf;
    double histPdfZero;
    RooDataSet* unBinnedData[3];
    RooDataHist* binnedData[6];
    int timeNdf[3];
    double timeChi[3];
    int specNdf[3];
    double specChi[3];
    string ifRed;
}fitInf;

isoItem isoB12={"B12",0,0.0291,0.011,0.8 ,3.e4,0,5.e6,632};
isoItem isoN12={"N12",0,0.0159,0.005,0.08,1.e3,0,5.e6,635};
isoItem isoC9 ={"C9" ,0,0.1825,0.051,0.5 ,5.e2,0,5.e6,800};
isoItem isoLi8={"Li8",0,1.21  ,0.511,5.0 ,2.e3,0,5.e6,601};
isoItem isoB8 ={"B8" ,0,1.11  ,0.511,5.0 ,1.e3,0,5.e6,909};
isoItem isoLi9={"Li9",0,0.1717,0.081,0.8 ,2.e3,0,5.e6,415};
isoItem isoHe8={"He8",0,0.2572,0.111,0.8 ,5.e2,0,5.e6,432};
isoItem isoBkg={"Bkg",0,0.    ,0.   ,50. ,5.e4,0,5.e6,416};

fitInf fitB12={"B12",1,1,1,5.5 ,20.0,0.002,0.100,5.5 ,20.0,0.002,0.100,{"B12","N12","C9","Bkg"},{"B12","N12","C9"},"N12"};
//fitInf fitB12={"B12",1,1,1,5.5 ,20.0,0.001,0.501,5.5 ,20.0,0.002,0.060,{"B12","N12","C9","He8","Li9","Bkg"},{"B12","N12","C9","He8","Li9"},"N12"};
fitInf fitN12={"N12",1,1,1,14.0,20.0,0.001,0.501,14.0,20.0,0.002,0.060,{"N12","C9","He8","Li9","Bkg"},{"N12","C9","He8","Li9"},"C9"};
fitInf fitLi8={"Li8",1,1,1,5.5 ,20.0,0.8  ,10.  ,5.5 ,20.0,0.6  ,4.0  ,{"Li8","B8","C9","Bkg"},{"Li8","B8","C9"},"B8"};
fitInf fitC9 ={"C9" ,1,1,1,12.0,20.0,0.15 ,2.0  ,12.0,20.0,0.2  ,0.6  ,{"C9" ,"B8","Li8","Bkg"},{"C9" ,"B8","Li8"},"Li8"};
map<string,isoItem*> iso;
map<string,fitInf*> fit;

void calDaqTime(fitInf* fitinf,string dataVer,string site)
{
    //livetime
    string runnum;
    TH1F* h[5];//the fifth is for daqtime
    Double_t totalTime[5]={0.};

    //B12 N12
    double tlivetime=0.;
    TH1F* showermuonNum[4]; 
    //TH1F* time2lastshowermuon[3];
    //TH1F* B12Result[5];
    //TH1F* hh[3];

    //===>get histograms from .root file for analyse
    string filename;
    filename.assign(dataVer);
    filename+="/";
    filename+=site;
    filename+="totalHisto_";
    filename+=dataVer;
    filename+=".root";
    TFile* f=new TFile(filename.c_str());

    //livetime
    for( int i=0 ; i<5 ; i++ )
    {

        stringstream hname;
        if( i==4 )
        {
            hname<<"LiveTime/DaqTime";
        } else
        {
            hname<<"LiveTime/AD";
            hname<<i+1;
            hname<<"LiveTime";	
        }
        h[i] = (TH1F*)f->Get(hname.str().c_str());
        if( !(h[i]) )
        {
            cout<<"ERROR : Can not get Hist : "<<hname.str()<<" from "<<site<<"/"<<runnum<<" ."<<endl;
            //return true;
            continue;
        }
    }

    //B12 N12,muon numbers in each slice.
    //for( int i=0 ; i<3 ; i++ )
    //{
    //TString hnameLi;
    //hnameLi="lidj/showermuonNumNoRed";
    ////hnameLi="lidj/showermuonNum";
    //hnameLi+=i+1;
    //showermuonNum[i]=(TH1F*)f->Get(hnameLi);
    //}
        showermuonNum[0]=(TH1F*)f->Get("lidj/muonEnergyI16");
        showermuonNum[1]=(TH1F*)f->Get("lidj/showermuonNumNoRed4");
        showermuonNum[2]=(TH1F*)f->Get("lidj/showermuonNumNoRed5");
        showermuonNum[3]=(TH1F*)f->Get("lidj/showermuonNumNoRed6");

    //===>calculate livetime and muonrate 
    TAxis *xaxis = h[0]->GetXaxis();
    int binNum = xaxis->GetNbins();
    //livetime
    int BinBound=h[0]->FindBin(1346428800);//2012.9.1 0:0:0
    double daqtDel=0.;
    if( BinBound!=0 )
    {
        if( site=="EH2" )
        {
            for(int j=1  ; j<=BinBound ; j++ )
            {
                h[1]->SetBinContent(j,0);
                daqtDel+=h[4]->GetBinContent(j);
            }
            for( int j=1 ; j<binNum ; j++ )
            {
                h[2]->SetBinContent(j,0);
                h[3]->SetBinContent(j,0);
            }
        }else if( site=="EH3" )
        {
            for(int j=1  ; j<=BinBound ; j++ )
            {
                h[3]->SetBinContent(j,0);
                daqtDel+=h[4]->GetBinContent(j);
            }
        }else
        {
            for( int j=1 ; j<binNum ; j++ )
            {
                h[2]->SetBinContent(j,0);
                h[3]->SetBinContent(j,0);
            }
        }
    }
    for( int i=0 ; i<5 ; i++ )
    {
        totalTime[i]=0.;
        for( int j=1 ; j<=binNum ; j++ )
        {
            totalTime[i]+=h[i]->GetBinContent(j);
        }

    }
    for( int i=0 ; i<4 ; i++ )
    {
        tlivetime+=totalTime[i];
    }
    fitinf->liveTime=tlivetime;
    double totalDaqtime=0.;
    if( site=="EH2" )
    {
        totalDaqtime=totalTime[4]*2-daqtDel;
    }else if( site=="EH3" )
    {
        totalDaqtime=totalTime[4]*4-daqtDel;
    } else
    {
        totalDaqtime=totalTime[4]*2;
    }
    fitinf->daqTime=totalDaqtime;
    //muon rate
    double NumMuon[3]={0.};
    double RateMuon[3]={0.};
    NumMuon[0]=showermuonNum[0]->Integral(1,showermuonNum[0]->FindBin(1500));
    fitinf->muonRate[0]=NumMuon[0]/totalDaqtime;
    for( int j=1 ; j<4 ; j++ )
    {
        for( int jbin=0 ; jbin<binNum ; jbin++ )
        {
            NumMuon[j]+=showermuonNum[j]->GetBinContent(jbin);
        }
        RateMuon[j]=NumMuon[j]/totalDaqtime;
        fitinf->muonRate[j]=RateMuon[j];
    }
}
void genIsoPdf(isoItem* myIso,fitInf* myfit)
{
    std::cout<<Form(">>> >>> prepare %-3s information ",myIso->isoName.c_str())<<endl;
    nameStr=Form("tau%s",myIso->isoName.c_str());
    if( myIso->isFixed )
    {
        myIso->fitTau = new RooRealVar(nameStr, nameStr, -myIso->tau_val,-myIso->tau_val_high,-myIso->tau_val_low, "s");
    }else
    {
        myIso->fitTau = new RooRealVar(nameStr, nameStr, -myIso->tau_val, "s");
    }
    nameStr=Form("N_%s",myIso->isoName.c_str());
    nameStr2=Form("total number of %s",myIso->isoName.c_str());
    myIso->tFitNum=new RooRealVar(nameStr,nameStr2,myIso->N_value,myIso->N_val_low,myIso->N_val_high);
    nameStr=Form("lambda%s",myIso->isoName.c_str());
    if( myIso->isoName=="Bkg" )
    {
        myIso->fitLambda=new RooFormulaVar(nameStr,nameStr,"@0 + @1",RooArgList(*myIso->fitTau, *rateMu));
    }else
    {
        myIso->fitLambda=new RooFormulaVar(nameStr,nameStr,"1/@0 + @1",RooArgList(*myIso->fitTau, *rateMu));
    }
    nameStr=Form("exp%s",myIso->isoName.c_str());
    nameStr2=Form("%s distribution",myIso->isoName.c_str());
    myIso->fitExpPdf=new RooExponential(nameStr,nameStr2, *xt, *myIso->fitLambda);

    nameStr=Form("%stimeTCutCoe",myIso->isoName.c_str());
    myIso->timeTCutCoe=new RooRealVar(nameStr,nameStr,1);

    nameStr=Form("%sspecTCutCoe",myIso->isoName.c_str());
    myIso->specTCutCoe=new RooRealVar(nameStr,nameStr,1);

    nameStr=Form("ExtendPdf%s",myIso->isoName.c_str());
    myIso->fitExtendPdf=new RooExtendPdf(nameStr,nameStr,*myIso->fitExpPdf,*(myIso->tFitNum)) ;
    if( myIso->isoName!="Bkg" )
    {
        nameStr=Form("%sspecHistogramVsp.d.f",myIso->isoName.c_str());
        TCanvas* ce = new TCanvas(nameStr,nameStr,800,400) ;
        ce->Divide(2) ;
        TFile* f=new TFile("IsoTheoreticalSpec.root","read");
        nameStr=Form("%sSpectraAfterCor",myIso->isoName.c_str());
        TH1F* h=(TH1F*)f->Get(nameStr);
        //h->Scale(10000000);
        ce->cd(1) ;  h->Draw() ;
        myIso->ehd=new RooDataHist("hd","hd",*xe,h);
        nameStr=Form("%shpdf",myIso->isoName.c_str());
        myIso->fitHistPdf=new RooHistPdf(nameStr,nameStr,*xe,*(myIso->ehd),2) ;
        RooPlot* framee = xe->frame(Title("spec histogram")) ;
        myIso->ehd->plotOn(framee,LineColor(kRed),DataError(RooAbsData::None));
        //ce->cd(2) ;  framee->Draw() ;

        double timeECoe=h->Integral(h->FindBin(myfit->timeElow),h->FindBin(myfit->timeEhigh))/h->Integral(1,h->GetNbinsX());
        cout<<"timeECoe  : "<<timeECoe<<endl;
        nameStr=Form("%stimeECutCoe",myIso->isoName.c_str());
        myIso->timeECutCoe=new RooRealVar(nameStr,nameStr,timeECoe);

        double specECoe=h->Integral(h->FindBin(myfit->specElow),h->FindBin(myfit->specEhigh))/h->Integral(1,h->GetNbinsX());
        cout<<"specECoe  : "<<specECoe<<endl;
        nameStr=Form("%sspecECutCoe",myIso->isoName.c_str());
        myIso->specECutCoe=new RooRealVar(nameStr,nameStr,specECoe);

        nameStr=Form("%seFitNum",myIso->isoName.c_str());
        myIso->eFitNum=new RooFormulaVar(nameStr,nameStr,"@0*@1*@2/(@3*@4)",RooArgList(*(myIso->tFitNum),*(myIso->specECutCoe),*(myIso->specTCutCoe),*(myIso->timeTCutCoe),*(myIso->timeECutCoe)));
        RooPlot* framee1 = xe->frame(Title("spec p.d.f")) ;
        myIso->fitHistPdf->plotOn(framee1,LineColor(kBlue)) ;
        ce->cd(2) ;  framee1->Draw() ;
        nameStr=Form("P4AB/dataEps/%sspecHistogramVsp.d.f.eps",myIso->isoName.c_str());
        ce->SaveAs(nameStr);
        cout<<"nameStr  : "<<nameStr<<endl;

    }
}
void genFitPdf(fitInf* fitinf)
{
    RooArgList timeFitComList;
    RooArgList specFitComList;
    RooArgList specFitNumList;
    for( int i=0 ; i<10 ; i++ )
    {
        if( fitinf->timeCom[i]!="" )
        {
            timeFitComList.add(*(iso[fitinf->timeCom[i]]->fitExtendPdf));
            //std::cout<<"iso["<<fitinf->timeCom[i].c_str()<<"]  : "<<iso[fitinf->timeCom[i].c_str()]<<endl;
            //std::cout<<"iso["<<fitinf->timeCom[i].c_str()<<"]->isoName  : "<<iso[fitinf->timeCom[i].c_str()]->isoName<<endl;
            fitinf->timeComMap.insert(map<string,isoItem*>::value_type(fitinf->timeCom[i].c_str(),iso[fitinf->timeCom[i].c_str()]));
            //std::cout<<"timeComMap["<<fitinf->timeCom[i].c_str() <<"]->isoName  : "<<fitinf->timeComMap[fitinf->timeCom[i].c_str()]->isoName<<endl;
        }
        if( fitinf->specCom[i]!="" )
        {
                specFitComList.add(*(iso[fitinf->specCom[i].c_str()]->fitHistPdf));
                specFitNumList.add(*(iso[fitinf->specCom[i].c_str()]->eFitNum));
            fitinf->specComMap.insert(map<string,isoItem*>::value_type(fitinf->specCom[i].c_str(),iso[fitinf->specCom[i].c_str()]));
            //std::cout<<"specComMap["<<fitinf->specCom[i].c_str() <<"]->isoName  : "<<fitinf->specComMap[fitinf->specCom[i].c_str()]->isoName<<endl;
        }
    }
    std::cout<<Form(">>> >>> prepare %-3s fit TimeFitPdf ",fitinf->mode.c_str())<<endl;
    fitinf->timeFitPdf=new RooAddPdf("timeModel","timeModel",timeFitComList) ;
    std::cout<<Form(">>> >>> prepare %-3s fit SpecFitPdf ",fitinf->mode.c_str())<<endl;
    fitinf->specFitPdf=new RooAddPdf("specMode","specMode",specFitComList,specFitNumList);
    int checkNum=10000;
    double checkWidth=(fitinf->specEhigh-fitinf->specElow)/checkNum;
    for( int i=1 ; i<checkNum; i++ )
    {
        if(fitinf->specFitPdf->getVal(*xe=fitinf->specElow+i*checkWidth)==0  )
        {
            std::cout<<"Find specFitPdf zero  at : "<<fitinf->specElow+i*checkWidth<<endl;
            fitinf->histPdfZero=fitinf->specElow+(i-1)*checkWidth;
            break;
        }

    }
}

void genData(fitInf* fitinf, string dataVer,string site)
{
    //TODO 3 slice for time fit
    std::cout<<Form(">>> >>> prepare %-3s fit data ",fitinf->mode.c_str())<<endl;
    nameStr=Form("%s/%siso_%s.root",dataVer.c_str(),site.c_str(),dataVer.c_str());
    TFile* f=new TFile(nameStr,"read");
    TH1F* ht[3];
    TH1F* he[3];
    TTree* t[7];
    if( fitinf->isNoRed )
    {
        fitinf->ifRed="NoRed";
    } else
    {
        fitinf->ifRed="";
    }
    for( int i=0 ; i<3 ; i++ )
    {
        if( fitinf->isbinned )
        {
            nameStr=Form("time2lastshowermuon%s%i_%0.1f_%0.1f",fitinf->ifRed.c_str(),i+3,fitinf->specElow,20.);
            if( i==0 )
            {
                nameStr=Form("I1time2lastshowermuon123_%0.1f_%0.1f",fitinf->specElow,20.); 
            }
            ht[i]=(TH1F*)f->Get(nameStr);
            if( !ht[i] )
            {
                std::cout<<"ERROR: Can't find histogram  : "<<nameStr<<endl;
                exit(0);
            }
            //cout<<"ht["<<i<<"]  : "<<ht[i]->GetEntries()<<endl;
            //cout<<"nameStr  : "<<nameStr<<endl;
            nameStr2=Form("%sSpec%sSlice%i_%0.1f_%0.1f",fitinf->mode.c_str(),fitinf->ifRed.c_str(),i+3,fitinf->specElow,20.);
            if( i==0 )
            {
                nameStr2=Form("%sI1SpecSlice123_%0.1f_%0.1f",fitinf->mode.c_str(),fitinf->specElow,20.);
            }
            he[i]=(TH1F*)f->Get(nameStr2);
            //cout<<"he["<<i<<"]  : "<<he[i]->GetEntries()<<endl;
            //cout<<"nameStr2  : "<<nameStr2<<endl;
            if( !he[i] )
            {
                std::cout<<"ERROR: Can't find histogram  : "<<nameStr2<<endl;
                exit(0);
            }

            //make sure that bincontent>=10 ,pdf in  xe range !=0 for MML fit;
            int hemax=he[i]->FindBin(fitinf->histPdfZero)-1;
            int hemin=he[i]->FindBin(fitinf->specElow);
            if( he[i]->GetBinLowEdge(hemin)<fitinf->specElow ) hemin+=1;
            xe->setRange(he[i]->GetBinLowEdge(hemin),he[i]->GetBinLowEdge(hemax+1));
            float heContent[1000]={0.};
            float heError[1000]={0.};
            float heEdge[1001]={0.};
            heEdge[0]=he[i]->GetBinLowEdge(hemin);
            double hebinCon=0.;
            double hebinErr2=0.;
            int henewBinNum=0;
            for( int bi=hemin ; bi<=hemax ; bi++ )
            {
                hebinCon+=he[i]->GetBinContent(bi);
                hebinErr2+=he[i]->GetBinError(bi)*he[i]->GetBinError(bi);
                if( hebinCon>=10. )
                {
                    heContent[henewBinNum]=hebinCon;
                    heError[henewBinNum]=sqrt(hebinErr2);
                    heEdge[++henewBinNum]=he[i]->GetBinLowEdge(bi+1);
                    cout<<" "<<henewBinNum<<"  : "<<heContent[henewBinNum-1]<<"+-"<<heError[henewBinNum-1]<<" "<<heEdge[henewBinNum]<<endl;
                    hebinCon=0.;
                    hebinErr2=0.;
                    continue;
                    
                }
                if( bi==hemax )
                {
                    do
                    {
                        henewBinNum--;
                        heContent[henewBinNum]=heContent[henewBinNum]+hebinCon;
                        hebinCon=heContent[henewBinNum];
                        heError[henewBinNum]=sqrt(heError[henewBinNum]*heError[henewBinNum]+hebinErr2);
                        hebinErr2=heError[henewBinNum]*heError[henewBinNum];
                    }while(hebinCon<10);
                    heEdge[++henewBinNum]=he[i]->GetBinLowEdge(hemax+1);
                }
            }
            cout<<"spec original bin number  : "<<hemax-hemin+1<<endl;
            cout<<"spec new bin number  : "<<henewBinNum<<endl;
            if(henewBinNum>=1000)
            {
                cout<<"ERROR : henewBinNum>=1000 "<<endl;
                exit(0);
            }
            TH1F* he4Fit= new TH1F("he4Fit","he4Fit",henewBinNum,heEdge);
            for( int bi=1 ; bi<=henewBinNum ; bi++ )
            {
                he4Fit->SetBinContent(bi,heContent[bi-1]);
                he4Fit->SetBinError(bi,heError[bi-1]);
            }
            fitinf->specNdf[i]=henewBinNum-fitinf->specComMap.size()-1;
            cout<<"he4Fit entries  : "<<he4Fit->Integral(1,henewBinNum)<<endl;
            cout<<"he4Fit entries  : "<<he4Fit->GetEntries()<<endl;

            //ht[i]->SetOption("E1");
            //ht[i]->Rebin(5);
            int htmax=ht[i]->FindBin(fitinf->timeThigh)-1;
            int htmin=ht[i]->FindBin(fitinf->timeTlow);
            cout<<"htmin  : "<<htmin<<endl;
            cout<<"ht[i]->GetBinLowEdge(htmin)  : "<<ht[i]->GetBinLowEdge(htmin)<<endl;
            cout<<"fitinf->timeElow  : "<<fitinf->timeElow<<endl;
            if( ht[i]->GetBinLowEdge(htmin)<fitinf->timeTlow ) htmin+=1;
            cout<<"htmin  : "<<htmin<<endl;
            xt->setRange(ht[i]->GetBinLowEdge(htmin),ht[i]->GetBinLowEdge(htmax+1));
            std::cout<<"h xt max  : "<<xt->getMax()<<endl;
            std::cout<<"h xt min  : "<<xt->getMin()<<endl;
            float htContent[1000]={0.};
            float htError[1000]={0.};
            float htEdge[1001]={0.};
            htEdge[0]=ht[i]->GetBinLowEdge(htmin);
            double htbinCon=0.;
            double htbinErr2=0.;
            int htnewBinNum=0;
            for( int bi=htmin ; bi<=htmax ; bi++ )
            {
                htbinCon+=ht[i]->GetBinContent(bi);
                htbinErr2+=ht[i]->GetBinError(bi)*ht[i]->GetBinError(bi);
                if( htbinCon>=10. )
                {
                    htContent[htnewBinNum]=htbinCon;
                    htError[htnewBinNum]=sqrt(htbinErr2);
                    htEdge[++htnewBinNum]=ht[i]->GetBinLowEdge(bi+1);
                    cout<<" "<<htnewBinNum<<"  : "<<htContent[htnewBinNum-1]<<"+-"<<htError[htnewBinNum-1]<<" "<<htEdge[htnewBinNum]<<endl;
                    htbinCon=0.;
                    htbinErr2=0.;
                    continue;
                    
                }
                if( bi==htmax )
                {
                    do
                    {
                        htnewBinNum--;
                        htContent[htnewBinNum]=htContent[htnewBinNum]+htbinCon;
                        htbinCon=htContent[htnewBinNum];
                        htError[htnewBinNum]=sqrt(htError[htnewBinNum]*htError[htnewBinNum]+htbinErr2);
                        htbinErr2=htError[htnewBinNum]*htError[htnewBinNum];
                    }while(htbinCon<10);
                    htEdge[++htnewBinNum]=ht[i]->GetBinLowEdge(htmax+1);
                }
            }
            cout<<"time original bin number  : "<<htmax-htmin+1<<endl;
            cout<<"time new bin number  : "<<htnewBinNum<<endl;
            if(htnewBinNum>=1000)
            {
                cout<<"ERROR : htnewBinNum>=1000 "<<endl;
                exit(0);
            }
            TH1F* ht4Fit= new TH1F("ht4Fit","ht4Fit",htnewBinNum,htEdge);
            for( int bi=1 ; bi<=htnewBinNum ; bi++ )
            {
                ht4Fit->SetBinContent(bi,htContent[bi-1]);
                ht4Fit->SetBinError(bi,htError[bi-1]);
            }
            fitinf->timeNdf[i]=htnewBinNum-fitinf->timeComMap.size()-1;



            fitinf->binnedData[i]=new RooDataHist(Form("%stime2lastmuonBinned",fitinf->mode.c_str()),"time2lastmuon binned data",*xt,ht4Fit);
            fitinf->binnedData[i+3]=new RooDataHist(Form("%sspecBinned",fitinf->mode.c_str()),"spec binned data",*xe,he4Fit);
            cout<<"numEntries  : "<<fitinf->binnedData[i+3]->numEntries()<<endl;
            //RooPlot* frame2 = xe->frame(Title("Time fit")) ;
            //fitinf->binnedData[i+3]->plotOn(frame2,DataError(RooAbsData::SumW2),MarkerStyle(9),DrawOption("Z")) ;
            //TCanvas* c = new TCanvas("c","c",400,300) ;
            //frame2->Draw() ;
            //nameStr2=Form("test%i.eps",i+1);
            //c->SaveAs(nameStr2);
            //delete c;
        }else
        {
//TODO:
            cout<<"Next will add unbinned fit  "<<endl;
            exit(0);
            nameStr2=Form("slice%s%i_%0.1f_%0.1f",fitinf->ifRed.c_str(),i+1,fitinf->specElow,20.);
            if( i==5 )
            {
                nameStr2=Form("slice%s6_%0.1f_%0.1f",fitinf->ifRed.c_str(),fitinf->specElow,20.);
            }
            t[i]=(TTree*)f->Get(nameStr2);
            fitinf->unBinnedData[i]=new RooDataSet(Form("%stime2lastmuon",fitinf->mode.c_str()),"time2lastmuon unbinned data",t[i],*xt);
        }
    }

    //f->Close();
}
void prepareInf( string dataVer,string site,string fitmode,bool doSimFit)
{
    //fit information
    //isotopes information
    iso.insert(map<string,isoItem*>::value_type("B12",&isoB12));
    iso.insert(map<string,isoItem*>::value_type("N12",&isoN12));
    iso.insert(map<string,isoItem*>::value_type("Li8",&isoLi8));
    iso.insert(map<string,isoItem*>::value_type("B8" ,&isoB8));
    iso.insert(map<string,isoItem*>::value_type("C9" ,&isoC9));
    iso.insert(map<string,isoItem*>::value_type("Li9",&isoLi9));
    iso.insert(map<string,isoItem*>::value_type("He8",&isoHe8));
    iso.insert(map<string,isoItem*>::value_type("Bkg",&isoBkg));
    std::cout<<"finish inserting iso information "<<endl;
    //isNoRed,isbinned
    fit.insert(map<string,fitInf*>::value_type("B12",&fitB12));
    fit.insert(map<string,fitInf*>::value_type("N12",&fitN12));
    fit.insert(map<string,fitInf*>::value_type("Li8",&fitLi8));
    fit.insert(map<string,fitInf*>::value_type("C9",&fitC9));
    std::cout<<"finish inserting fit information "<<endl;
    std::cout<<"fitmode  : "<<fitmode<<endl;

    xt=new RooRealVar("xt","x for time fit",(fit[fitmode]->timeThigh-fit[fitmode]->timeTlow)/2,fit[fitmode]->timeTlow,fit[fitmode]->timeThigh);
    xe=new RooRealVar("xe","x for spec fit",(fit[fitmode]->specEhigh-fit[fitmode]->specElow)/2,fit[fitmode]->specElow,fit[fitmode]->specEhigh);
    rateMu=new RooRealVar("rateMu","rateMu",0);
    std::cout<<"initialize xt setRange "<<endl;
    std::cout<<"xt : "<<xt->getVal()<<endl;
    std::cout<<"initialize xe setRange "<<endl;
    std::cout<<"xe : "<<xe->getVal()<<endl;

    fit[fitmode]->doSimulFit=doSimFit;
    calDaqTime(fit[fitmode],dataVer,site);
    for( int i=0 ; i<10 ; i++ )
    {
        if( fit[fitmode]->timeCom[i]!="" )
        {
            genIsoPdf(iso[fit[fitmode]->timeCom[i]],fit[fitmode]); //calculate kinds of coe.
        }
    }
    genFitPdf(fit[fitmode]);genData(fit[fitmode],dataVer,site);//genFitPdf should be first ,genData() will use a value form it.
    cout<<">>> >>> fit range : "<<endl;
    std::cout<<"fit["<<fitmode<<"] time fit Tlow  : "<<xt->getMin()<<endl;
    std::cout<<"fit["<<fitmode<<"] time fit Thigh  : "<<xt->getMax()<<endl;
    std::cout<<"fit["<<fitmode<<"] spec fit Elow  : "<<xe->getMin()<<endl;
    std::cout<<"fit["<<fitmode<<"] spec fit Ehigh  : "<<xe->getMax()<<endl;
}

int doFit(int siteNum,string dataVer,string fitMode,bool doSimFit,double fitLowRange,double fitHighRange)
{

    //===>initialize variable
    //bool anaIso=1;//N12/B12 bg
    //bool anaE=0;
    int ADNumOfSite[3]={0};
    double timeIsolationRatio[3]={0.0152,0.0446,0.8134};
    if( dataVer.find("13")!=string::npos || dataVer.find("14")!=string::npos || dataVer.find("15")!=string::npos || dataVer.find("16")!=string::npos)
    {
        ADNumOfSite[0]=2;
        ADNumOfSite[1]=2;
        ADNumOfSite[2]=4;
    } else
    {
        if( dataVer.find("11")!=string::npos || dataVer.find("12")!=string::npos )
        {
            ADNumOfSite[0]=2;
            ADNumOfSite[1]=1;
            ADNumOfSite[2]=3;
        }
    }
    string site;
    if( siteNum==1 )
    {
        site="EH1";
    } else if(siteNum==2)
    {
        site="EH2";
    }else if (siteNum==3)
    {
        site="EH3";
    }else
    {
        site="EH1";
    }

    //===>initialize histograms

    //===> fit iso
    std::cout<<"begin to analyse isotopes "<<endl;
    prepareInf(dataVer,site,fitMode,doSimFit);   
    if( fitLowRange!=0 && fitHighRange!=0 )
    {
        fit[fitMode]->timeTlow=fitLowRange;
        fit[fitMode]->timeThigh=fitHighRange;
    }
    std::cout<<"xt range  : "<<fit[fitMode]->timeTlow<<" ~ "<<fit[fitMode]->timeThigh <<endl;
    //return 1;
    if(!fit[fitMode]->doSimulFit)
    {
        bool draw6slices=1;
        TString xTitle[3]={"20Mev~1.5GeV","1.5~2.5GeV",">2.5GeV"};
        RooFitResult* fitres;

        TCanvas* c;
        if( draw6slices )
        {
            c = new TCanvas("c","c",1200,300) ; 
            c->Divide(3,1);
            gStyle->SetEndErrorSize(5.0);
            gStyle->SetMarkerSize(0.1);
        }
        for( int ihist=0 ; ihist<1 ; ihist++ )
        {
            for( int j=0 ; j<3 ; j++ )
            {
                std::cout<<"now is  : "<<j+1<<endl;
                rateMu->setVal(-fit[fitMode]->muonRate[j]);
                xt->setRange("timeTRange",fit[fitMode]->timeTlow,fit[fitMode]->timeThigh);
                xt->setRange("specTRange",fit[fitMode]->specTlow,fit[fitMode]->specThigh);
                xt->setRange(0,1.e6);
                for( map<string,isoItem*>::iterator it=fit[fitMode]->timeComMap.begin() ; it!=fit[fitMode]->timeComMap.end() ; it++ )
                {
                    RooAbsReal* timeTTmpCoe =it->second->fitExpPdf->createIntegral(*xt,NormSet(*xt),Range("timeTRange"));
                    std::cout<<"isoname  : "<<it->second->isoName<<endl;
                    std::cout<<"timeTTmpCoe  : "<<timeTTmpCoe->getVal()<<endl;
                    it->second->timeTCutCoe->setVal(timeTTmpCoe->getVal());
                    //it->second->timeTCutCoe=new RooRealVar(nameStr,nameStr,timeTTmpCoe->getVal());
                    delete timeTTmpCoe;
                    RooAbsReal* specTTmpCoe =it->second->fitExpPdf->createIntegral(*xt,NormSet(*xt),Range("specTRange"));
                    std::cout<<"specTTmpCoe  : "<<specTTmpCoe->getVal()<<endl;
                    nameStr=Form("%sspecTCutCoe",it->second->isoName.c_str());
                    it->second->specTCutCoe->setVal(specTTmpCoe->getVal());
                    //it->second->specTCutCoe=new RooRealVar(nameStr,nameStr,specTTmpCoe->getVal());
                    delete specTTmpCoe;
                }
                xt->setRange(fit[fitMode]->timeTlow,fit[fitMode]->timeThigh);
                //Num_tot[j]=hh[j]->Integral(1,hh[j]->FindBin(fit[fitMode]->timeThigh));
                if( fit[fitMode]->isbinned  )
                {
                    //fit[fitMode]->timeFitPdf->fitTo(*(fit[fitMode]->binnedData[j]),Save(),PrintLevel(-1));
                    //RooAbsReal* nll = fit[fitMode]->timeFitPdf->createNLL(*(fit[fitMode]->binnedData[j])) ;
                    //RooMinuit m(*nll) ; 
                    //m.migrad() ; 
                    //m.hesse() ;
                    //m.contour(*(iso[fit[fitMode]->con]->tFitNum),*(iso[fitMode]->tFitNum),1,2,3) ; 
                }
                if( fit[fitMode]->isbinned )
                {
                    fitres = fit[fitMode]->timeFitPdf->fitTo(*(fit[fitMode]->binnedData[j]),Save(),PrintLevel(-1));
                    RooChi2Var* timeChi2=new RooChi2Var("timeChi2", "timeChi2", *(fit[fitMode]->timeFitPdf), *(fit[fitMode]->binnedData[j]));
                    fit[fitMode]->timeChi[j]=timeChi2->getVal();
                }else
                {
                    fitres = fit[fitMode]->timeFitPdf->fitTo(*(fit[fitMode]->unBinnedData[j]),Save(),PrintLevel(-1),NumCPU(20));
                }
                fitres->Print();

                for( map<string,isoItem*>::iterator it=fit[fitMode]->timeComMap.begin() ; it!=fit[fitMode]->timeComMap.end() ; it++ )
                {
                    it->second->tNumInSlice[j]=it->second->tFitNum->getValV(); 
                    it->second->tNumErrInSlice[j]=it->second->tFitNum->getError(); 
                    cout<<it->second->isoName<<" tNumInSlice["<<j+1<<"]  : "<<it->second->tNumInSlice[j]<<endl;
                    //it->second->tFitNum->setError(0); 
                    if( it->second->isoName!="Bkg" )
                    {
                        it->second->realNumInSlice[j]=it->second->tFitNum->getValV()/it->second->timeECutCoe->getValV()/it->second->timeTCutCoe->getValV(); 
                        it->second->realNumErrInSlice[j]=it->second->tFitNum->getError()/it->second->timeECutCoe->getValV()/it->second->timeTCutCoe->getValV(); 
                    }
                    it->second->tauInSlice[j]=it->second->fitTau->getValV();
                    it->second->tauErrInSlice[j]=it->second->fitTau->getError();
                    if( j==0)
                    {
                        it->second->totalRealNum=it->second->realNumInSlice[j]/timeIsolationRatio[siteNum-1];
                        it->second->totalRealNumErr2=it->second->realNumErrInSlice[j]*it->second->realNumErrInSlice[j]/(timeIsolationRatio[siteNum-1]*timeIsolationRatio[siteNum-1]);
                    }else
                    {
                        it->second->totalRealNum+=it->second->realNumInSlice[j];
                        it->second->totalRealNumErr2+=it->second->realNumErrInSlice[j]*it->second->realNumErrInSlice[j];
                        
                    }
                }

                if( draw6slices )
                {
                    c->cd(j+1);
                    RooPlot* mesframe = xt->frame() ;
                    xTitle[j]+="        time since last muon (s)";
                    mesframe->GetXaxis()->SetTitle(xTitle[j]);
                    mesframe->GetYaxis()->SetTitle("Entries");
                    TLegend* leg = new TLegend(0.3,0.6,0.95,0.94);
                    //leg->SetTextSize(0.05);

                    if( fit[fitMode]->isbinned )
                    {
                        fit[fitMode]->binnedData[j]->plotOn(mesframe) ;
                    } else
                    {
                        fit[fitMode]->unBinnedData[j]->plotOn(mesframe) ;
                    }
                    fit[fitMode]->timeFitPdf->plotOn(mesframe,Name("sum"));
                    leg->SetHeader(Form("---%s Fit---",fitMode.c_str()));
                    for(map<string,isoItem*>::iterator it=fit[fitMode]->timeComMap.begin();it!=fit[fitMode]->timeComMap.end() ; it++)
                    {
                        fit[fitMode]->timeFitPdf->plotOn(mesframe,Components(*(it->second->fitExpPdf)),Name(Form("%s",it->second->isoName.c_str())),LineStyle(kDashed),LineColor(it->second->linecolor)) ;
                        leg->AddEntry(mesframe->findObject(Form("%s",it->second->isoName.c_str())),Form("%-3s ~ %0.2f",it->second->isoName.c_str(),it->second->tNumInSlice[j]),"l");
                    }
                    if( fit[fitMode]->isbinned )
                    {
                        TString timeChiStr = Form("Chi/ndf  =  %2.2f / %i = %2.2f", fit[fitMode]->timeChi[j],fit[fitMode]->timeNdf[j],fit[fitMode]->timeChi[j]/fit[fitMode]->timeNdf[j]);
                        leg->AddEntry((TObject*)0,timeChiStr,"");
                    }
                    leg->SetFillColor(10);
                    mesframe->Draw();
                    leg->Draw();
                    //gPad->SetLogy();
                }
            }
            nameStr2=Form("%s/simFitEps/%s%stimeFit%s.eps",dataVer.c_str(),site.c_str(),fitMode.c_str(),fit[fitMode]->ifRed.c_str());
            c->SaveAs(nameStr2);
            if( ihist==0 )
            {
                for( int i=0 ; i<ADNumOfSite[siteNum-1] ; i++ )
                {
                    //tlivetime+=totalTime[i];
                }
            }else
            {
                //tlivetime=totalTime[ihist-1];
            }
            std::cout<<" "<<endl;
            std::cout<<"livetime  : "<<fit[fitMode]->liveTime/(24*3600)<<endl;
            cout<<" "<<endl;
            std::cout<<" Slice |";
            //for(  map<string,isoItem*>::iterator it=fit[fitMode]->timeComMap.begin() ; it!=fit[fitMode]->timeComMap.end() ; it++ )
            for( int k=0;k<(int)(fit[fitMode]->timeComMap.size());k++) 
            {
                std::cout<<Form("    %-6s |",iso[fit[fitMode]->timeCom[k]]->isoName.c_str());
            }
            std::cout<<" MuonRate |";
            if( fit[fitMode]->isbinned )
            {
                std::cout<<"    Chi/ndf ";
            }
            std::cout<<endl;

            for( int j=0 ; j<3 ; j++ )
            {
                std::cout<<Form(" %5i |",j+1);//"│" | | ｜｜
                for( int k=0;k<(int)(fit[fitMode]->timeComMap.size());k++) 
                {
                    std::cout<<Form(" %9.1f |",fit[fitMode]->timeComMap[fit[fitMode]->timeCom[k]]->tNumInSlice[j]);//can't use "iso[fit[fitMode]->timeCom[k]]",not is reference type.
                }
                std::cout<<Form(" %8.5f |",fit[fitMode]->muonRate[j]);
                if( fit[fitMode]->isbinned )
                {
                    std::cout<<Form(" %8.2f/%i= %5.2f",fit[fitMode]->timeChi[j],fit[fitMode]->timeNdf[j],fit[fitMode]->timeChi[j]/fit[fitMode]->timeNdf[j]); 
                }
                std::cout<<endl;

            }

        }

    }else
    {
        //do simultaneous fit
        for( int j=0 ; j<3 ; j++ )
        {
            rateMu->setVal(-fit[fitMode]->muonRate[j]);
            xt->setRange(0,1.e6);
            xt->setRange("timeTRange",fit[fitMode]->timeTlow,fit[fitMode]->timeThigh);
            xt->setRange("specTRange",fit[fitMode]->specTlow,fit[fitMode]->specThigh);
            for( map<string,isoItem*>::iterator it=fit[fitMode]->timeComMap.begin() ; it!=fit[fitMode]->timeComMap.end() ; it++ )
            {
                RooAbsReal* timeTTmpCoe =it->second->fitExpPdf->createIntegral(*xt,NormSet(*xt),Range("timeTRange"));
                std::cout<<"isoname  : "<<it->second->isoName<<endl;
                std::cout<<"timeTTmpCoe  : "<<timeTTmpCoe->getVal()<<endl;
                it->second->timeTCutCoe->setVal(timeTTmpCoe->getVal());
                //it->second->timeTCutCoe=new RooRealVar(nameStr,nameStr,timeTTmpCoe->getVal());
                delete timeTTmpCoe;
                RooAbsReal* specTTmpCoe =it->second->fitExpPdf->createIntegral(*xt,NormSet(*xt),Range("specTRange"));
                std::cout<<"specTTmpCoe  : "<<specTTmpCoe->getVal()<<endl;
                nameStr=Form("%sspecTCutCoe",it->second->isoName.c_str());
                it->second->specTCutCoe->setVal(specTTmpCoe->getVal());
                //it->second->specTCutCoe=new RooRealVar(nameStr,nameStr,specTTmpCoe->getVal());
                delete specTTmpCoe;
            }
            xt->setRange(fit[fitMode]->timeTlow,fit[fitMode]->timeThigh);
            RooCategory sample("sample","sample") ;
            sample.defineType("spec");
            sample.defineType("time");
            RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
            simPdf.addPdf(*(fit[fitMode]->specFitPdf),"spec") ;
            simPdf.addPdf(*(fit[fitMode]->timeFitPdf),"time") ;

            std::cout<<"xt max  : "<<xt->getMax()<<endl;
            std::cout<<"xt min  : "<<xt->getMin()<<endl;
            std::cout<<"xe max  : "<<xe->getMax()<<endl;
            std::cout<<"xe min  : "<<xe->getMin()<<endl;

            //xt->setBins(fit[fitMode]->binnedData[j]->numEntries());
            //xe->setBins(fit[fitMode]->binnedData[j+3]->numEntries());
            xt->setBins((xt->getMax()-xt->getMin())/0.001);
            xe->setBins((xe->getMax()-xe->getMin())/0.25);
            RooDataHist combData("combData","combined data",RooArgSet(*xe,*xt),Index(sample),Import("spec",*(fit[fitMode]->binnedData[j+3])),Import("time",*(fit[fitMode]->binnedData[j]))) ;

            //1.
            //RooAbsReal* nll = simPdf.createNLL(combData) ;
            //RooMinuit m(*nll) ; 
            //m.migrad() ; 
            //m.hesse() ;
            //m.contour(*(iso[fit[fitMode]->con]->tFitNum),*(iso[fitMode]->tFitNum),1,2,3) ; 
            //2.
            //simPdf.fitTo(combData,SumW2Error(kTRUE),PrintEvalErrors(10)) ;
            //3.
            RooAbsReal* nlltime=fit[fitMode]->timeFitPdf->createNLL(*fit[fitMode]->binnedData[j]);
            RooAbsReal* nllspec=fit[fitMode]->specFitPdf->createNLL(*fit[fitMode]->binnedData[j+3]);
            //RooAbsReal* nlltime=fit[fitMode]->timeFitPdf->createChi2(*fit[fitMode]->binnedData[j],DataError(RooAbsData::SumW2));
            //RooAbsReal* nllspec=fit[fitMode]->specFitPdf->createChi2(*fit[fitMode]->binnedData[j+3],DataError(RooAbsData::SumW2));
            //RooAbsReal* nlltime=fit[fitMode]->timeFitPdf->createChi2(*fit[fitMode]->binnedData[j]);
            //RooAbsReal* nllspec=fit[fitMode]->specFitPdf->createChi2(*fit[fitMode]->binnedData[j+3]);
            RooAddition* nllTotal=new RooAddition("timeAndSpec","timeAndSpec",RooArgSet(*nlltime,*nllspec));
            RooMinuit mTotal(*nllTotal) ; 
            ////mTotal.optimizeConst(kTRUE);
            ////mTotal.setProfile(kTRUE);
            ////mTotal.setVerbose(kFALSE);
            mTotal.migrad() ; 
            mTotal.hesse() ;
            //mTotal.contour(*(iso[fit[fitMode]->con]->tFitNum),*(iso[fitMode]->tFitNum),1,2,3) ; 
            //mTotal.migrad() ; 
            //mTotal.hesse() ;

            std::cout<<"begin to get value "<<endl;
            std::cout<<"timeComMap.size  : "<<fit[fitMode]->timeComMap.size()<<endl;
            std::cout<<"specComMap.size  : "<<fit[fitMode]->specComMap.size()<<endl;
            cout<<"slice  "<<j+1<<" FITResult: "<<endl;
            for( map<string,isoItem*>::iterator it=fit[fitMode]->timeComMap.begin() ; it!=fit[fitMode]->timeComMap.end() ; it++ )
            {
                it->second->tNumInSlice[j]=it->second->tFitNum->getValV(); 
                cout<<it->second->isoName<<" tNumInSlice["<<j+1<<"]  : "<<it->second->tNumInSlice[j]<<" +- "<<it->second->tFitNum->getError()<<" "<<it->second->tFitNum->getError()/it->second->tNumInSlice[j]<<endl;
                it->second->tNumErrInSlice[j]=it->second->tFitNum->getError(); 
                it->second->tauInSlice[j]=it->second->fitTau->getValV();
                it->second->tauErrInSlice[j]=it->second->fitTau->getError();
                if( it->second->isoName!="Bkg" )
                {
                    it->second->realNumInSlice[j]=it->second->tFitNum->getValV()/it->second->timeECutCoe->getValV()/it->second->timeTCutCoe->getValV(); 
                    cout<<"timefit-realNumInSlice["<<j+1<<"]  : "<<it->second->realNumInSlice[j]<<endl;
                    it->second->realNumErrInSlice[j]=it->second->tFitNum->getError()/it->second->timeECutCoe->getValV()/it->second->timeTCutCoe->getValV(); 
                }
                if( j==0)
                {
                    it->second->totalRealNum=it->second->realNumInSlice[j]/timeIsolationRatio[siteNum-1];
                    cout<<"totalRealNum  : "<<it->second->totalRealNum<<endl;
                    it->second->totalRealNumErr2=it->second->realNumErrInSlice[j]*it->second->realNumErrInSlice[j]/(timeIsolationRatio[siteNum-1]*timeIsolationRatio[siteNum-1]);
                }else
                {
                    it->second->totalRealNum+=it->second->realNumInSlice[j];
                    it->second->totalRealNumErr2+=it->second->realNumErrInSlice[j]*it->second->realNumErrInSlice[j];

                }
            }
            for( map<string,isoItem*>::iterator it=fit[fitMode]->specComMap.begin() ; it!=fit[fitMode]->specComMap.end() ; it++ )
            {
                    it->second->eNumInSlice[j]=it->second->eFitNum->getValV(); 
                    //it->second->eNumErrInSlice[j]=it->second->eFitNum->getError(); 
                    cout<<"specfit-realNumInSlice["<<j<<"]  : "<<it->second->eFitNum->getValV()/it->second->specECutCoe->getValV()/it->second->specTCutCoe->getValV()<<endl;
            }
            RooChi2Var* timeChi2=new RooChi2Var("timeChi2", "timeChi2", *(fit[fitMode]->timeFitPdf), *(fit[fitMode]->binnedData[j]));
            fit[fitMode]->timeChi[j]=timeChi2->getVal();
            delete timeChi2;
            RooChi2Var* specChi2=new RooChi2Var("specChi2", "specChi2", *(fit[fitMode]->specFitPdf), *(fit[fitMode]->binnedData[j+3]));

            fit[fitMode]->specChi[j]=specChi2->getVal();
            cout<<">>> time fit chi/ndf  : "<<fit[fitMode]->timeChi[j]<<"/"<<fit[fitMode]->timeNdf[j]<<endl;
            cout<<">>> spec fit chi/ndf  : "<<fit[fitMode]->specChi[j]<<"/"<<fit[fitMode]->specNdf[j]<<endl;
            std::cout<<"begin to plot simultaneous fit "<<endl;
            delete specChi2;
            // Plot all data tagged as time sample
            RooPlot* frame1 = xt->frame(Title("Time fit")) ;
            combData.plotOn(frame1,Cut("sample==sample::time"),DataError(RooAbsData::SumW2),MarkerStyle(9),DrawOption("Z")) ;
            simPdf.plotOn(frame1,Slice(sample,"time"),ProjWData(sample,combData),LineWidth(1)) ;
            TLegend* leg1 = new TLegend(0.5,0.6,0.89,0.89);
            for(map<string,isoItem*>::iterator it=fit[fitMode]->timeComMap.begin();it!=fit[fitMode]->timeComMap.end() ; it++)
            {
                simPdf.plotOn(frame1,Slice(sample,"time"),Components(*(it->second->fitExpPdf)),ProjWData(sample,combData),Name(Form("%s",it->second->isoName.c_str())),LineStyle(kDotted),LineColor(it->second->linecolor),LineWidth(1)) ;
                leg1->AddEntry(frame1->findObject(Form("%s",it->second->isoName.c_str())),Form("%-3s ~ %0.0f +- %0.0f",it->second->isoName.c_str(),it->second->tNumInSlice[j],it->second->tNumErrInSlice[j]),"l");
            }
            TString timeChiStr;
            timeChiStr = Form("Chi2/ndf   %2.0f/%i", fit[fitMode]->timeChi[j],fit[fitMode]->timeNdf[j]);
            leg1->AddEntry((TObject*)0,timeChiStr,"");
            leg1->SetFillColor(10);
            leg1->SetBorderSize(0);
            //gPad->SetLogy();
            // The same plot for the spec sample slice
            RooPlot* frame2 = xe->frame(Title("Spectrum fit")) ;
            combData.plotOn(frame2,Cut("sample==sample::spec"),DataError(RooAbsData::SumW2),MarkerStyle(9),DrawOption("Z")) ;
            simPdf.plotOn(frame2,Slice(sample,"spec"),ProjWData(sample,combData),LineWidth(1)) ;
            TLegend* leg2 = new TLegend(0.5,0.6,0.89,0.89);
            for(map<string,isoItem*>::iterator it=fit[fitMode]->specComMap.begin();it!=fit[fitMode]->specComMap.end() ; it++)
            {
                simPdf.plotOn(frame2,Slice(sample,"spec"),Components(*(it->second->fitHistPdf)),ProjWData(sample,combData),Name(Form("%sSpec",it->second->isoName.c_str())),LineStyle(kDotted),LineColor(it->second->linecolor),LineWidth(1)) ;
                //leg2->AddEntry(frame1->findObject(Form("%s",it->second->isoName.c_str())),Form("%-3s ~ %0.0f ",it->second->isoName.c_str(),it->second->eNumInSlice[j]),"l");
            }
            TString specChiStr;
            specChiStr = Form("Chi2/ndf   %2.0f/%i", fit[fitMode]->specChi[j],fit[fitMode]->specNdf[j]);
            leg2->AddEntry((TObject*)0,specChiStr,"");
            //leg2->SetHeader(Form("---%s Fit---",fitMode.c_str()));
            leg2->SetFillColor(10);
            leg2->SetBorderSize(0);
            nameStr=Form("%s%s%sSimultaneousFit%sSlice%i",site.c_str(),dataVer.c_str(),fitMode.c_str(),fit[fitMode]->ifRed.c_str(),j+1);
            TCanvas* c = new TCanvas(nameStr,nameStr,800,400) ;
            c->Divide(2) ;
            c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;leg1->Draw();
            c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;leg2->Draw();
            nameStr2=Form("%s/simFitEps/%s%sSimultaneousFit%sSlice%i.eps",dataVer.c_str(),site.c_str(),fitMode.c_str(),fit[fitMode]->ifRed.c_str(),j+1);
            c->SaveAs(nameStr2);
            c->Draw();
            //if(j==4) c->Draw();
            //delete frame1;
            //delete leg1;
            //delete frame2;
            //delete leg2;
            //delete c;

        }
    }
    cout<<" "<<endl;
    std::cout<<"      |       Total number      |     Rate(/day/AD) "<<endl;
    double volumeEff=1.;
    for( int k=0;k<(int)(fit[fitMode]->timeComMap.size());k++) 
    {
        fit[fitMode]->timeComMap[fit[fitMode]->timeCom[k]]->rate=fit[fitMode]->timeComMap[fit[fitMode]->timeCom[k]]->totalRealNum/(fit[fitMode]->liveTime/(24*3600))/volumeEff;
        fit[fitMode]->timeComMap[fit[fitMode]->timeCom[k]]->rateErr=sqrt(fit[fitMode]->timeComMap[fit[fitMode]->timeCom[k]]->totalRealNumErr2)/(fit[fitMode]->liveTime/(24*3600))/volumeEff;
        double numFit=fit[fitMode]->timeComMap[fit[fitMode]->timeCom[k]]->totalRealNum;
        //double numErrFit=sqrt(numFit);
        double numErrFit=sqrt(fit[fitMode]->timeComMap[fit[fitMode]->timeCom[k]]->totalRealNumErr2);
        double rateFit=numFit/(fit[fitMode]->liveTime/(24*3600))/volumeEff;
        double rateErrFit=numErrFit/(fit[fitMode]->liveTime/(24*3600))/volumeEff;
        std::cout<<Form(" %-3s  | %10.1f +- %8.1f  | %7.1f +- %7.1f",fit[fitMode]->timeComMap[fit[fitMode]->timeCom[k]]->isoName.c_str(),numFit,numErrFit,rateFit,rateErrFit)<<endl;
    }

    cout<<" "<<endl;
    std::cout<<"      |       Yield(1E-7)       |     Rate(/day/AD/t) "<<endl;
    for( int k=0;k<(int)(fit[fitMode]->timeComMap.size());k++) 
    {
        double numFit=fit[fitMode]->timeComMap[fit[fitMode]->timeCom[k]]->totalRealNum;
        //double numErrFit=sqrt(numFit);
        double numErrFit=sqrt(fit[fitMode]->timeComMap[fit[fitMode]->timeCom[k]]->totalRealNumErr2);
        double yield=numFit/(0.860*258*fit[fitMode]->muonRate[3]*fit[fitMode]->daqTime)/volumeEff;
        double yieldErr=numErrFit/(0.860*258*fit[fitMode]->muonRate[3]*fit[fitMode]->daqTime)/volumeEff;
        double ratePerT=numFit/(fit[fitMode]->liveTime/(24*3600))/41./volumeEff;
        double ratePerTErr=numErrFit/(fit[fitMode]->liveTime/(24*3600))/41./volumeEff;
        std::cout<<Form(" %-3s  | %10.2f +- %8.2f  | %7.2f +- %7.2f",fit[fitMode]->timeComMap[fit[fitMode]->timeCom[k]]->isoName.c_str(),yield*1.e7,yieldErr*1.e7,ratePerT,ratePerTErr)<<endl;
    }
    cout<<" "<<endl;
    for( int j=0 ; j<3 ; j++ )
    {
        std::cout<<"muon rate  : "<<fit[fitMode]->muonRate[j]<<endl;
    }
    std::cout<<"original total muon rate  : "<<fit[fitMode]->muonRate[3]<<endl;

    //===>print live time result
    std::cout<<""<<endl;
    std::cout<<site <<"'s infomation : "<<endl;
    std::cout<<""<<endl;
    //for( int i=0 ; i<ADNumOfSite[siteNum-1] ; i++ )
    //{
    //std::cout<<"Total AD"<<i+1<<"LiveTime                         : "<<totalTime[i]/(24*3600)<<endl;
    //}
    //std::cout<<"Total DaqTime : "<<totalTime[4]/(24*3600)<<" day" <<endl;
    //std::cout<<""<<endl;

    //===>write into .root file
    //string rootname=site;
    //rootname+="FitResult_"+dataVer+".root";
    //TFile* file = new TFile(rootname.c_str(),"RECREATE");
    //file->cd();
    //for( int i=0 ; i<5 ; i++ )
    //{
    //h[i]->Write();
    //}
    //if( anaIso )
    //{
    //for( int i=0 ; i<5 ; i++ )
    //{
    //hh[i]->Write();
    //}
    ////B12Result[0]->Write();
    //}
    //file->Close();
    return 0;
}

vector<string> checkdata(string dataVer)
{
    string dataIsGood="1";
    string verNum[6]={"11","12","13","14","15","16"};
    vector<string> dataVerReal;
    for( int i=0 ; i<6 ; i++ )
    {
        int pos=0;
        string TmpVer;
        TmpVer=dataVer.substr(pos);
        while( pos!=-1 )
        {
            pos=TmpVer.find(verNum[i]);
            if( pos!=-1 )
            {
                string subdataVer;
                subdataVer="P"+verNum[i];
                subdataVer+=toupper(*(TmpVer.substr(pos+2,1).c_str()));
                for( vector<string>::iterator it=dataVerReal.begin() ; it!=dataVerReal.end() ; ++it  )
                {
                    if( *it==subdataVer )
                    {
                        dataVerReal.erase(it);
                        it--;
                    }
                }
                dataVerReal.push_back(subdataVer);
                TmpVer=TmpVer.substr(pos+3);
            }
        }
    }
    sort(dataVerReal.begin(),dataVerReal.end());

    std::cout<<" "<<endl;
    std::cout<<"Find following data version will be analysed : "<<endl;
    std::cout<<" "<<endl;
    string runlistSiteNum[3]={"EH1","EH2","EH3"};
    vector<string> runlistFileContent[3];
    for( int i=0 ; i<3 ; i++ )
    {
        runlistFileContent[i].clear();
    }
    string verSuf;
    for( vector<string>::iterator it=dataVerReal.begin() ; it!=dataVerReal.end() ; ++it )
    {
        std::cout<<*it<<endl;
        verSuf+="_"+*it;
        //check data path
        string datadir="/afs/ihep.ac.cn/users/l/lidj/largedata/IsotopesAna/";
        DIR *dir=NULL;
        datadir+=*it;
        dir = opendir(datadir.c_str());
        if( dir==NULL )
        {
            std::cout<<" !!! data dir doesn't exist : "<<datadir<<endl;
            dataIsGood="0";
        } else
        {
            std::cout<<" data dir exist : "<<datadir<<endl;
            closedir(dir);
        }


        std::cout<<" "<<endl;
    }
    //check merged file 
    for( int i=0 ; i<3 ; i++ )
    {
        string runlistName;
        runlistName.assign(verSuf,1,verSuf.size()-1);
        runlistName+="/";
        runlistName+=runlistSiteNum[i];
        runlistName+="totalHisto";
        runlistName+=verSuf;
        //runlistName+=*it;
        runlistName+=".root";
        TFile *ff = new TFile(runlistName.c_str());
        if( ff->IsZombie() )
        {
            std::cout<<" !!! merged file doesn't exist  : "<<runlistName<<endl;
            /*
               if( *it!="P12A" )
               {
               dataIsGood="0";
               }
               */
        } else
        {
            std::cout<<" merged file exist  : "<<runlistName<<endl;
        }
        ff->Close();
    }
    verSuf=verSuf.substr(1);
    std::cout<<"data version : "<<verSuf<<endl;

    dataVerReal.push_back(verSuf);
    dataVerReal.push_back(dataIsGood);

    return dataVerReal; 
}

//int IsoNumFit(string dataVer,int siteNum=0,string FitMode="B12",bool doSimFit=1,double FitLowRange=0.,double FitHighRange=0.)
int main(int argc, char *argv[])
{
    string dataVer;
    int siteNum=0;
    string FitMode;
    bool doSimFit=0;
    double FitLowRange=0.;
    double FitHighRange=0.;
    std::cout<<"argc  : "<<argc<<endl;
    dataVer=argv[1];
    std::cout<<"dataVer  : "<<dataVer<<endl;
    siteNum=atoi(argv[2]);
    std::cout<<"siteNum  : "<<siteNum<<endl;
    FitMode=argv[3];
    std::cout<<"FitMode  : "<<FitMode<<endl;
    doSimFit=atoi(argv[4]);
    std::cout<<"doSimFit  : "<<doSimFit<<endl;
    if( argc==7 )
    {
        FitLowRange=atof(argv[5]);
        std::cout<<"FitLowRange  : "<<FitLowRange<<endl;
        FitHighRange=atof(argv[6]);
        std::cout<<"FitHighRange  : "<<FitHighRange<<endl;
    }
    /*
    //vector<string> dataVerVec=checkdata(dataVer);
    std::cout<<" "<<endl;
    std::cout<<" "<<endl;
    std::cout<<" "<<endl;
    std::cout<<" "<<endl;
    //dataVer=dataVerVec[dataVerVec.size()-2];

    if( dataVerVec[dataVerVec.size()-1]=="0" )
    {
        std::cout<<"  The dataVerion["<<dataVer<<"] or data is wrong ,please check! "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        return 0;
    }
    */

    if( siteNum==0 )
    {
        for( int i=1 ; i<=3 ; i++ )
        {
            std::cout<<"====> begin to analyse EH"<<i<<"'s DaqTime, N12/B12 "<<endl;
            std::cout<<"dataVersion  : "<<dataVer<<endl;
            doFit(i,dataVer,FitMode,doSimFit,FitLowRange,FitHighRange);
            std::cout<<" "<<endl;
            std::cout<<" "<<endl;
            std::cout<<" "<<endl;
            std::cout<<" "<<endl;
        }

    }else if(siteNum>=1 && siteNum<=3)
    {
        std::cout<<"====> begin to analyse EH"<<siteNum<<"'s DaqTime,N12/B12 "<<endl;
        std::cout<<"dataVersion  : "<<dataVer<<endl;
        doFit(siteNum,dataVer,FitMode,doSimFit,FitLowRange,FitHighRange);
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
    }
    std::cout<<"ALL DONE!!! "<<endl;
    std::cout<<" "<<endl;
    std::cout<<" "<<endl;
    std::cout<<" "<<endl;
    return 0;
}

