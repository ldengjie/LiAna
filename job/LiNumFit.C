#include  <iostream>
#include  <TStyle.h>
#include  <dirent.h>//DIR
#include  <TCanvas.h>
#include  <TMath.h>
#include  <TLegend.h>
#include  <TH1F.h>
#include  <TFile.h>
#include  <fstream>
#include  <sstream>
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include    "RooPlot.h"
#include "RooFitResult.h"
#include  <string>
#include  <TChain.h>
#include  <TROOT.h>
#include  <vector>
using namespace RooFit;

int fitHisto(int siteNum,string dataVer)
{

//===>initialize variable
	bool anaLiHe=true;//He8/Li9 bg
    int ADNumOfSite[3]={0};
    int daqHistNum=5;
    if( dataVer.find("13")!=-1 || dataVer.find("14")!=-1 || dataVer.find("15")!=-1 || dataVer.find("16")!=-1)
    {
        ADNumOfSite[0]=2;
        ADNumOfSite[1]=2;
        ADNumOfSite[2]=4;
        daqHistNum=5;
    } else
    {
        if( dataVer.find("11")!=-1 || dataVer.find("12")!=-1 )
        {
            ADNumOfSite[0]=2;
            ADNumOfSite[1]=1;
            ADNumOfSite[2]=3;
            daqHistNum=4;
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
	//livetime
	string runnum;
	TH1F* h[daqHistNum];//the fifth is for daqtime
	Double_t totalTime[5]={0.};

	//Li9 He8
	double n98total=0.;
	//double in98total=0.;
    double ifitn98total=0.;
    double n98Rate=0.;
    double in98Rate=0.;
    double tlivetime=0.;
	TH1F* showermuonNum[6]; 
	TH1F* time2lastshowermuon[6];
	TH1F* hh[6];

//===>get histograms from .root file for analyse
	string filename;
    filename=dataVer;
    filename+="/";
    filename+=site;
    filename+="mergedHist_";
    filename+=dataVer;
    filename+=".root";
    TFile* f=new TFile(filename.c_str());

		//livetime
		for( int i=0 ; i<daqHistNum ; i++ )
		{

			stringstream hname;
			if( i==daqHistNum-1 )
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
				cout<<"Can not get Hist : "<<hname.str()<<" from "<<site<<"/"<<runnum<<" ."<<endl;
				//return true;
                continue;
			}
		}

		//Li9 He8
        if(anaLiHe)
		{
            for( int i=0 ; i<6 ; i++ )
            {
                TString hnameLi;
                hnameLi="lidj/showermuonNum";
                hnameLi+=i+1;
                //hnameLi+="NoRed";
                showermuonNum[i]=(TH1F*)f->Get(hnameLi);
                hnameLi="lidj/time2lastshowermuon";
                hnameLi+=i+1;
                hnameLi+="4Li";
                //hnameLi+="NoRed";
                time2lastshowermuon[i]=(TH1F*)f->Get(hnameLi);
            }
		}
	
//===>analyse 
	TAxis *xaxis = h[0]->GetXaxis();
	int binNum = xaxis->GetNbins();
    int BinBound=h[0]->FindBin(1346428800);//2012.9.1 0:0:0
	//livetime
    if( BinBound!=0 )
    {
        if( site=="EH2" )
        {
            for(int j=1  ; j<=BinBound ; j++ )
            {
                h[1]->SetBinContent(j,0);
            }
        }
        if( site=="EH3" )
        {
            for(int j=1  ; j<=BinBound ; j++ )
            {
                h[3]->SetBinContent(j,0);
            }
        }
    }
    double totalTimeSpecical=0.;
    for( int j=BinBound ; j<=binNum ; j++ )
    {
        totalTimeSpecical+=h[4]->GetBinContent(j);
    }

    for( int i=0 ; i<5 ; i++ )
    {
        totalTime[i]=0.;
        for( int j=1 ; j<=binNum ; j++ )
        {
            totalTime[i]+=h[i]->GetBinContent(j);
        }

    }
    double totalDaqtime=0.;
    if( site=="EH2" )
    {
        totalDaqtime=totalTime[4]+totalTimeSpecical;
    }else if( site=="EH3" )
    {
        totalDaqtime=totalTime[4]*3+totalTimeSpecical;
    } else
    {
        totalDaqtime=totalTime[4]*2;
    }

	//Li He
    std::cout<<"begin to analyse Li "<<endl;
	if(anaLiHe)
	{
        int hnum=time2lastshowermuon[0]->FindBin(100);
        //std::cout<<"0.001  : "<<time2lastshowermuon[0]->FindBin(0.001)<<endl;
        //std::cout<<"100  : "<<time2lastshowermuon[0]->FindBin(100)<<endl;
        //std::cout<<"hnum  : "<<hnum<<endl;
        for( int k=0 ; k<6 ; k++ )
        {
            TString hhname="slice";
            hhname+=k+1;
            hh[k]=new TH1F(hhname,"hh",hnum-2,0.001,100);
                for( int i=1 ; i<hnum-2 ; i++ )
                {
                    hh[k]->SetBinContent(i,time2lastshowermuon[k]->GetBinContent(i+1));
                }
                hh[k]->SetOption("E1");
        }
        double showerTh[7] = {0.02, 0.5, 1.5, 2.5, 3.5, 4.5, 5.0};
        TH1F* LiResult[5];
        LiResult[0]= new TH1F("LiResult", "All", 6, showerTh);
        LiResult[1]= new TH1F("LiResult", "AD1", 6, showerTh);
        LiResult[2]= new TH1F("LiResult", "AD2", 6, showerTh);
        LiResult[3]= new TH1F("LiResult", "AD3", 6, showerTh);
        LiResult[4]= new TH1F("LiResult", "AD4", 6, showerTh);
		double NumMuon[6]={0.};
		double RateMuon[6]={0.};
		double NumIbd[6]={0.};
        
		RooRealVar x("x","x",0.001,40., "s");
		RooRealVar tauLi9("tauLi9", "tauLi9", -0.257, "s");
		RooRealVar tauHe8("tauHe8", "tauHe8", -0.172, "s");
		RooRealVar rateMu("rateMu","rate of showermuon",-0.1,-10., 0.,"Hz");
		RooRealVar Li9Frac("Li9Frac","Li9's Frac", 0.0);//R
		RooRealVar N98("N98","total number of Li9 and He8",500.,0.,1.e5);
		//RooRealVar N98("N98","total number of Li9 and He8",3e2, 1e1, 2.e5);
		RooFormulaVar lambdaLi9("lambdaLi9","lambdaLi9","1/@0 + @1",RooArgList(tauLi9, rateMu));
		RooFormulaVar lambdaHe8("lambdaHe8","lambdaHe8","1/@0 + @1",RooArgList(tauHe8, rateMu));
		RooFormulaVar coeLi9("coeLi9","coeLi9","@0 * @1 ",RooArgList(N98,Li9Frac));
		RooFormulaVar coeHe8("coeHe8","coeHe8","@0 * ( 1 - @1 )",RooArgList(N98,Li9Frac));
		RooRealVar NIbd("NIbd","number of background",3e2, 1e1, 2.e6);
		RooExponential expLi9("expLi9", "Li9 distribution", x, lambdaLi9);
		RooExponential expHe8("expHe8", "He8 distribution", x, lambdaHe8);
        double Ibdrate[3]={0.0076,0.0067,0.00082};
        RooRealVar rateIbd("rateIbd","rateIbd",-Ibdrate[siteNum-1],-10,0.,"Hz");
        RooFormulaVar lambdaIbd("lambdaIbd","lambdaIbd","@0 + @1",RooArgList(rateIbd, rateMu));
		RooExponential expIbd("expIbd","Ibd distribution",x,rateMu);
		//RooExponential expIbd("expIbd","Ibd distribution",x,lambdaIbd);
		
		RooFitResult* fitres;
		RooAddPdf* sum;
		RooDataHist* data;
        RooPlot* mesframe[6];
		double n98[41]={0.};
		double in98[41]={0.};
		double nIbd[41]={0.};
		double inIbd[41]={0.};
		double minNLL[41]={0.};
        double rateMuValue[41]={0.};
        double rateMuErr[41]={0.};
		//double minNl[6]={0.};
		int minindex[6]={0};
		double n98fit[6]={0.};
		double in98fit[6]={0.};
        double rateMufit[6]={0.};
		double binwidth=0.;
        TH1F* lh[5][6];
        TString lname[5]={"","AD1","AD2","AD3","AD4"};
        int lcolor[5]={4,3,2,6,5};
        bool draw6slice=false;
        TCanvas* c;
        TString xTitle[6]={"slice 1 : 0.02~0.5GeV","slice 2 : 0.5~1.5GeV","slice 3 : 1.5~2.5GeV","slice 4 : 2.5~3.5GeV","slice 5 : 3.5~4.5GeV","slice 6 : >4.5GeV"};
        if( draw6slice )
        {
	        c = new TCanvas("c","c",2000,800) ; 
	        c->Divide(3,2);
	        gStyle->SetEndErrorSize(5.0);
	        gStyle->SetMarkerSize(0.1);
	        //gStyle->SetHistLineWidth(1);
	        if( 0 )
	        {
	            // it's strange ,if you can't draw following six figs in the big 3*2 pad,you can fist use these code,run once,don't exit the ROOT,then you can draw six figs!!!
			    data = new RooDataHist("data", "data", x, time2lastshowermuon[0]);
				sum  = new RooAddPdf("sum","sum",RooArgList(expLi9, expHe8, expIbd),RooArgList(coeLi9, coeHe8, NIbd));
		        RooPlot* mesframe = x.frame() ;
		        data->plotOn(mesframe) ;
		        c->cd(1);
		        mesframe->Draw();
	            return true;
	        }
        }
        for( int ihist=0 ; ihist<1 ; ihist++ )
        {
			double minNl[6]={0.};
			for( int j=0 ; j<6 ; j++ )
			{
	            std::cout<<"now is  "<< j+1<<endl;
	            NumIbd[j]=hh[j]->Integral(1,hh[j]->FindBin(40));
	            hh[j]->Rebin(100);
				data = new RooDataHist("data", "data", x, hh[j]);
				sum  = new RooAddPdf("sum","sum",RooArgList(expLi9, expHe8, expIbd),RooArgList(coeLi9, coeHe8, NIbd));
				for( int i=0 ; i<41 ; i++ )
				{
					Li9Frac.setVal(i*0.025);
					fitres = sum->fitTo((*data),Save(),PrintLevel(-1),Extended(kTRUE));
	                //fitres->Print();
					n98[i]=N98.getVal(0);
					in98[i]=N98.getError();
	                nIbd[i]=NIbd.getVal(0);
	                inIbd[i]=NIbd.getError();
	                rateMuValue[i]=rateMu.getVal(0);
	                rateMuErr[i]=rateMu.getError();
					minNLL[i] = fitres->minNll();
					if( minNLL[i]<minNl[j] )
					{
						minNl[j]=minNLL[i];
						minindex[j]=i;
					}
				}
				n98fit[j]=n98[minindex[j]];
				in98fit[j]=in98[minindex[j]];
	            rateMufit[j]=rateMuValue[minindex[j]];
	            if( draw6slice )
	            {
		            //for drawing fit figure
					Li9Frac.setVal(minindex[j]*0.025);
					fitres = sum->fitTo((*data),Save(),PrintLevel(-1),Extended(kTRUE));
		            //fitres->Print();
		            mesframe[j] = x.frame() ;
		            data->plotOn(mesframe[j]) ;
		            sum->plotOn(mesframe[j]);
		            sum->plotOn(mesframe[j],Components(expIbd),LineStyle(kDashed),LineColor(kGreen)) ;
		            sum->plotOn(mesframe[j],Components(RooArgSet(expLi9,expHe8)),LineStyle(kDashed),LineColor(kRed)) ;
		            xTitle[j]+="        time since last muon (s)";
		            mesframe[j]->GetXaxis()->SetTitle(xTitle[j]);
		            mesframe[j]->GetYaxis()->SetTitle("Entries");
		            c->cd(j+1);
		            mesframe[j]->Draw();
		            //gPad->SetLogy();
	            }
			}
	        
			n98total=(n98fit[0]+n98fit[1]+n98fit[2]+(n98fit[3]+n98fit[4]+n98fit[5])*exp(-1/0.257))/0.678;
			ifitn98total=sqrt(in98fit[0]*in98fit[0]+in98fit[1]*in98fit[1]+in98fit[2]*in98fit[2]+((in98fit[3]+in98fit[4]+in98fit[5])*exp(-1/0.257)*(in98fit[3]+in98fit[4]+in98fit[5])*exp(-1/0.257)))/0.678;
	        if( ihist==0 )
	        {
	            for( int i=0 ; i<ADNumOfSite[siteNum-1] ; i++ )
	            {
	                tlivetime+=totalTime[i];
	            }
	        }else
	        {
	            tlivetime=totalTime[ihist-1];
	        }
	        n98Rate=n98total/(tlivetime/(24*3600));
	        in98Rate=sqrt((ifitn98total/(tlivetime/(24*3600)))*(ifitn98total/(tlivetime/(24*3600)))+(n98total/(tlivetime/(24*3600))/2)*(n98total/(tlivetime/(24*3600))/2));
	        for( int j=0 ; j<6 ; j++ )
	        {
	            LiResult[ihist]->SetBinContent(j+1,n98fit[j]);
	            LiResult[ihist]->SetBinError(j+1,in98fit[j]);
			    for( int jbin=0 ; jbin<binNum ; jbin++ )
	            {
	                NumMuon[j]+=showermuonNum[j]->GetBinContent(jbin);
	            }
				RateMuon[j]=NumMuon[j]/totalDaqtime;
	            //std::cout<<"RateMuon[j]  : "<<RateMuon[j]<<endl;
	            //std::cout<<"RateMuon0[j]  : "<<RateMuon0[j]<<endl;
	            std::cout<<"  R "<<minindex[j]*0.025<<"   n98["<<minindex[j]<<"] "<<n98fit[j]<<"   in98["<<minindex[j]<<"] "<< in98fit[j]<<"   fitRateMu["<<minindex[j]<<"] "<<rateMufit[j]<<"  realRateMu "<<RateMuon[j]<<"   NumTol "<<NumIbd[j]<<endl;
	        }
	        std::cout<<"n98Num   : "<<n98total<<" +- "<<ifitn98total<<" Rate  : "<<n98Rate<<" +- "<< in98Rate<<endl;
        }
        
            /*
        //TCanvas* c4 = new TCanvas("c4", "c4", 600, 400);
        //gStyle->SetEndErrorSize(0.0);
        LiResult[0]->SetMarkerStyle(20);
		LiResult[0]->SetMarkerColor(lcolor[0]);
		LiResult[0]->SetLineColor(lcolor[0]);
		LiResult[0]->SetMarkerSize(1.0);
		LiResult[0]->SetMinimum(0);
		LiResult[0]->SetMaximum(LiResult[0]->GetMaximum()*1.5);
		LiResult[0]->GetXaxis()->SetTitle("Muon visible energy (GeV)");
		LiResult[0]->GetYaxis()->SetTitle("Fitted ^{9}Li events");
		LiResult[0]->SetTitle("");
        LiResult[0]->SetStats(kFALSE);
		LiResult[0]->Draw("EP");
        TLegend *legend=new TLegend(.6,.65,.79,.89);
        TString legendLabel=dataVer+" "+site+" "+"All";
        legend->AddEntry(LiResult[0],legendLabel,"lp");

        for( int ihist=1 ; ihist<ADNumOfSite[siteNum-1]+1; ihist++ )
        {
            LiResult[ihist]->SetMarkerStyle(20);
		    LiResult[ihist]->SetMarkerColor(lcolor[ihist]);
		    LiResult[ihist]->SetLineColor(lcolor[ihist]);
		    LiResult[ihist]->SetMarkerSize(1.0);
            LiResult[ihist]->SetStats(kFALSE);
		    LiResult[ihist]->Draw("EPsame");
            legend->AddEntry(LiResult[ihist],lname[ihist],"lp");
        }
        legend->SetFillColor(0);
        legend->Draw();
        */
	}	
	
    //f->Close();
//===>print result
	std::cout<<""<<endl;
	std::cout<<site <<"'s infomation : "<<endl;
	std::cout<<""<<endl;
	for( int i=0 ; i<ADNumOfSite[siteNum-1] ; i++ )
	{
		std::cout<<"Total AD"<<i+1<<"LiveTime                         : "<<totalTime[i]/(24*3600)<<endl;
	}
	std::cout<<"Total DaqTime : "<<totalTime[daqHistNum-1]/(24*3600)<<" day" <<endl;
	std::cout<<""<<endl;

	if(anaLiHe)
	{
	    std::cout<<" "<<endl;
        /*
	    std::cout<<"Li9/He8 "<<endl;
	    std::cout<<"n98total  : "<<n98total <<" +- "<<in98total <<" Rate:"<<n98total/totalTime[0]*24*3600 <<" +- "<<in98total/totalTime[0]*24*3600 <<endl;
	    std::cout<<"n98total0(without rpc)  : "<<n98total0 <<" +- "<<in98total0 <<" Rate:"<<n98total0/totalTime0[0]*24*3600 <<" +- "<<in98total0/totalTime0[0]*24*3600 <<endl;
        */
	}

//===>write into .root file
	string rootname=site;
	rootname+="FitResult_"+dataVer+".root";
	TFile* file = new TFile(rootname.c_str(),"RECREATE");
    file->cd();
    for( int i=0 ; i<daqHistNum ; i++ )
    {
        h[i]->Write();
    }
    
    if( anaLiHe )
    {
		//for( int i=0 ; i<ADNumOfSite[siteNum-1] ; i++ )
		for( int i=0 ; i<6 ; i++ )
		{
			hh[i]->Write();
		}
    }
	file->Close();
	f->Close();
	return 0;
}

int LiNumFit(string dataVer,int siteNum=0)
{
    std::cout<<" "<<endl;
    std::cout<<" "<<endl;
    std::cout<<" "<<endl;
    std::cout<<" "<<endl;

    if( siteNum==0 )
    {
        for( int i=1 ; i<=3 ; i++ )
        {
            
            std::cout<<"====> begin to analyse EH"<<i<<"'s DaqTime, He8/Li9 "<<endl;
            std::cout<<"dataVersion  : "<<dataVer<<endl;
            fitHisto(i,dataVer);
            std::cout<<" "<<endl;
            std::cout<<" "<<endl;
            std::cout<<" "<<endl;
            std::cout<<" "<<endl;
        }
        
    }else if(siteNum>=1 && siteNum<=3)
    {
    
        std::cout<<"====> begin to analyse EH"<<siteNum<<"'s DaqTime,He8/Li9 "<<endl;
        std::cout<<"dataVersion  : "<<dataVer<<endl;
        fitHisto(siteNum,dataVer);
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
    }
    std::cout<<"ALL DONE!!! "<<endl;
    return 0;
}
