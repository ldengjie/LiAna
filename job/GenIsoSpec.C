#include  <iostream>
#include  <TChain.h>
#include  <TROOT.h>
#include  <string>
using namespace std;
int anaTree(int siteNum,string dataVer)
{
    gROOT->ProcessLine(".L SingleTree.C+");
    TChain chain("Tree/LiTree");
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

    string chainFile;

    //if( dataVer=="P4AB" )
    //{
    //chainFile="/afs/ihep.ac.cn/users/l/lidj/largedata/LiAna/P14A/";
    //chainFile+=site;
    //chainFile+="/*LiAna.root";
    //chain.Add(chainFile.c_str());
    //chainFile="/afs/ihep.ac.cn/users/l/lidj/largedata/LiAna/P14B/";
    //chainFile+=site;
    //chainFile+="/*LiAna.root";
    //chain.Add(chainFile.c_str());
    //}else
    //{
    //chainFile="/afs/ihep.ac.cn/users/l/lidj/largedata/LiAna/";
    //chainFile+=dataVer;
    //chainFile+="/";
    //chainFile+=site;
    //chainFile+="/*LiAna.root";
    //std::cout<<"chainFile  : "<<chainFile<<endl;
    //chain.Add(chainFile.c_str());
    //}

    chainFile=dataVer;
    chainFile+="/";
    chainFile+=site;
    chainFile+="totalTree_";
    chainFile+=dataVer;
    chainFile+=".root";
    std::cout<<"chainFile  : "<<chainFile<<endl;
    chain.Add(chainFile.c_str());

    //chainFile="/afs/ihep.ac.cn/users/l/lidj/largedata/LiAna/P14A/EH1/run21344*LiAna.root";
    std::cout<<dataVer<<" : "<<chain.GetEntries()<<endl;
    if( dataVer=="P12E" )
    {
        chainFile="/afs/ihep.ac.cn/users/l/lidj/largedata/LiAna/";
        chainFile+="P12A";
        chainFile+="/";
        chainFile+=site;
        chainFile+="/*LiAna.root";
        chain.Add(chainFile.c_str());
        std::cout<<dataVer<<"+P12A : "<<chain.GetEntries()<<endl;
    }
    int evtnum=chain.GetEntries();
    TString  siteAndDataVer=site+dataVer+"Li9";
    siteAndDataVer+=evtnum;
    std::cout<<"site And DataVer  : "<<siteAndDataVer<<endl;
    std::cout<<"lidj total single number  : "<<evtnum<<endl;

    chain.Process("SingleTree",siteAndDataVer.Data());
    return 0;
}
int GenIsoSpec(string dataVer,int sitenum=0)
{
    if( !(dataVer=="P13A"||dataVer=="P12E"||dataVer=="P14A"||dataVer=="P4AB") )
    {
        std::cout<<"  The dataVerion is wrong ,please check! "<<endl;
    }
    if( sitenum==0 )
    {
        anaTree(1,dataVer);
        anaTree(2,dataVer);
        anaTree(3,dataVer);
    }
    if( sitenum==1 )
    {
        anaTree(1,dataVer);
    }
    if( sitenum==2 )
    {
        anaTree(2,dataVer);
    }
    if( sitenum==3 )
    {
        anaTree(3,dataVer);
    }
    return 0;
}
