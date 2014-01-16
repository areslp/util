#include <vector>
#include <set>
#include <fstream>
#include <list>
#include <map>
#include "util.h"

#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <pcl/common/common_headers.h>
#include <pcl/common/common_headers.h>
#include <pcl/features/normal_3d.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>



#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/if.hpp>

#include <boost/filesystem.hpp>
using namespace boost::filesystem;

using namespace boost::lambda;
using namespace std;

namespace po = boost::program_options;


int main(int argc, char const* argv[])
{
    std::string infile,orifile;

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("in", po::value<string>(&infile), "input file name")
        ("ori", po::value<string>(&orifile), "ori model file")
        ;
    po::positional_options_description pos;
    pos.add("in", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
                      options(desc).positional(pos).run(), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    
    if(!vm.count("ori")) {
        cout << "ori file must be set! exit...";
        return 1;
    } 

    //read points PointNormal type
    typedef pcl::PointCloud<pcl::PointNormal>::Ptr points_ptr;
    points_ptr points1(new pcl::PointCloud<pcl::PointNormal>);
    points_ptr points2(new pcl::PointCloud<pcl::PointNormal>);
    Util::read_pointnormals(infile, points1);
    Util::read_pointnormals(orifile, points2);
    assert(points1->size()==points2->size());

    int pNum=points1->size();

    int eNumber=0; //法矢估计错误较大的点数

    int nNumber=0; //法矢为空的点数
    
	int rNumber=0; //旋转过法矢的点数
    std::vector<int> eIdx(pNum,0);
    std::vector<int> nIdx(pNum,0);
    std::vector<int> rIdx(pNum,0);

    std::vector<float> eVals(pNum,0); //记录估计法矢与真实法矢之间的误差

	// for rendering
	for (int i = 0; i < points1->size(); i++) {
		pcl::PointNormal& p1=points1->points[i];
		pcl::PointNormal& p2=points2->points[i];
		Eigen::Vector3f n1=Eigen::Vector3f(p1.normal_x,p1.normal_y,p1.normal_z); 
		Eigen::Vector3f n2=Eigen::Vector3f(p2.normal_x,p2.normal_y,p2.normal_z); 
		if (n1.norm()==0) {
			continue;
		}
		n1.normalize();
		n2.normalize();
		// cout<<n1[0]<<" "<<n1[1]<<" "<<n1[2]<<endl;
		// cout<<n2[0]<<" "<<n2[1]<<" "<<n2[2]<<endl;
		float tmp=n1.dot(n2);

		//若n2为空，则跳过

		if (tmp!=tmp) {
			// std::cout<<"?"<<std::endl;
			continue;
		}
		if (tmp<0) {
			p1.normal_x=-p1.normal_x; 
			p1.normal_y=-p1.normal_y; 
			p1.normal_z=-p1.normal_z; 
		}
	}


    // for bias normals

    for (int i = 0; i < points1->size(); i++) {
        pcl::PointNormal& p1=points1->points[i];
        pcl::PointNormal& p2=points2->points[i];
        Eigen::Vector3f n1=Eigen::Vector3f(p1.normal_x,p1.normal_y,p1.normal_z); 
        Eigen::Vector3f n2=Eigen::Vector3f(p2.normal_x,p2.normal_y,p2.normal_z); 
        n1.normalize();
        n2.normalize();
        Eigen::Vector3f pn=n1.cross(n2);
        pn.normalize();
        //angel between sn and cn
        float tmp=n1.dot(n2);
        if (tmp>1) {
            tmp=1;
        }else if(tmp<-1){
            tmp=-1;
        }
        float agl=acos(tmp);
        // std::cout<<"angle:"<<agl<<std::endl;

        if (agl!=agl) {
            //if agl is nan, p1.normal() is 0
            // std::cout<<"agl is nan"<<std::endl;
            //空法矢

            nIdx[i]=1;
            eIdx[i]=1;
            nNumber++;
            eNumber++;
            eVals[i]=-1;
            continue;
        }

        // 若夹角与90相差10度，则认为是正好落在了另一个面上

        if (!(abs(agl-PI/2)<0.1745)) {
            eVals[i]=(agl>PI/2)?(PI-agl):agl;
            //大于30度认为是误差较大的点

            if (eVals[i]>PI/6) {
                // cout<<"agl:"<<agl<<endl;
                // cout<<"eVals[i]:"<<eVals[i]<<endl;
                // cout<<"eNumber++"<<endl;
                eNumber++;
                eIdx[i]=1;
            }
            continue;
        }

        rNumber++;
        rIdx[i]=1;
        eVals[i]=abs(agl-PI/2);
        // 旋转则说明误差夹角会小于10度，所以不用再考略eIdx

        //rotate matrix

        //Eigen::Matrix3f m;
        //float x=pn[0];
        //float y=pn[1];
        //float z=pn[2];
        //m << cos(agl)+( 1-cos(agl)*x*x ), (1-cos(agl))*x*y-sin(agl)*z, (1-cos(agl))*x*z+sin(agl)*y,
        //  (1-cos(agl))*y*x+sin(agl)*z, cos(agl)+(1-cos(agl))*(y*y), (1-cos(agl))*y*z-sin(agl)*x,
        //  (1-cos(agl))*z*x-sin(agl)*y, (1-cos(agl))*z*y+sin(agl)*x, cos(agl)+(1-cos(agl))*(z*z);
        //// std::cout<<m<<std::endl;
        //Eigen::Vector3f nn=m*n1;
        //// std::cout<<"nn[0]:"<<nn[0]<<"nn[1]:"<<nn[1]<<"nn[2]:"<<nn[2]<<std::endl;
        //p1.normal_x=nn[0];
        //p1.normal_y=nn[1];
        //p1.normal_z=nn[2];

    }

    
    float maxVal=*std::max_element(eVals.begin(),eVals.end());
    //将法矢为空的点误差值置为最大误差

    std::for_each( eVals.begin(), eVals.end(), 
            if_then(boost::lambda::_1 == -1, boost::lambda::_1=maxVal) );

    


    // 显示用模型

    std::string outfile="consist_"+infile;
    ofstream ofile(outfile.c_str());
    for(int i = 0; i < points1->size(); i++) {
        pcl::PointNormal p = points1->points[i];
        ofile << p.x << " " << p.y << " " << p.z << " " << p.normal_x << " " << p.normal_y << " " << p.normal_z << std::endl;
    }
    ofile.close();


    // Quality模型

    outfile=infile+"q";
    ofile.open(outfile.c_str());
    for(int i = 0; i < points1->size(); i++) {
        pcl::PointNormal p = points1->points[i];
        ofile << p.x << " " << p.y << " " << p.z << " " << p.normal_x << " " << p.normal_y << " " << p.normal_z << " " << eVals[i] << std::endl;
    }
    ofile.close();


    // 误差较大点列表

    outfile="elist_"+infile;
    ofile.open(outfile.c_str());
    for(int i = 0; i < eIdx.size(); i++) {
        int tmp=eIdx[i];
        if (tmp==1) {
            ofile <<(i+1)<< std::endl;
        }
    }
    ofile.close();

    // 旋转点列表

    outfile="rlist_"+infile;
    ofile.open(outfile.c_str());
    for(int i = 0; i < rIdx.size(); i++) {
        int tmp=rIdx[i];
        if (tmp==1) {
            ofile <<(i+1)<< std::endl;
        }
    }
    ofile.close();

    // 空法矢点列表

    outfile="nlist_"+infile;
    ofile.open(outfile.c_str());
    for(int i = 0; i < nIdx.size(); i++) {
        int tmp=nIdx[i];
        if (tmp==1) {
            ofile <<(i+1)<< std::endl;
        }
    }
    ofile.close();

    
    // std::cout<<"points size is:"<<pNum<<std::endl;
    // cout<<"eNumber:"<<eNumber<<endl;
    // cout<<"nNumber:"<<nNumber<<endl;
    // cout<<"rNumber:"<<rNumber<<endl;
    // cout<<"Null point percent:"<<(float)nNumber/(float)pNum<<endl;
    // 计算平均误差

    float aVals=0;
    std::for_each( eVals.begin(), eVals.end(), var(aVals)+=boost::lambda::_1);
    aVals=aVals/(float)pNum;


    //分析参数
    
    int idx=infile.find_last_of(".");
    std::string filename=infile.substr(0,idx);

    //vector<string> splitVars;
    //boost::split(splitVars,filename,boost::is_any_of("_")); 
    //map<string,string> kvmap;
    //for (int i = 0; i < splitVars.size(); i++) {
    //    string sv=splitVars[i];
    //    vector<string> svs;
    //    boost::split(svs,sv,boost::is_any_of(",")); 
    //    assert(svs.size()==2);
    //    kvmap.insert(make_pair(svs[0],svs[1]));
    //}
    //path path_ori(orifile);

    //cout<<path_ori.filename().c_str()<<" "<<pNum<<" "<<kvmap["l"]<<" "<<kvmap["tp"]<<" "<<kvmap["ng"]<<" "<<kvmap["i"]<<" "<<kvmap["li"]
    //    <<" "<<eNumber<<" "<<nNumber<<" "<<rNumber<<" "<<aVals<<endl;

	cout<<filename<<"point number:"<<pNum<<",error:"<<" "<<eNumber<<",null:"<<nNumber<<",rotate:"<<rNumber<<" "<<aVals<<endl;

    return 0;
}
