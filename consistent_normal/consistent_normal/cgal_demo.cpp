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
#include <pcl/kdtree/kdtree_flann.h>


#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/if.hpp>

#include <boost/filesystem.hpp>
using namespace boost::filesystem;

using namespace boost::lambda;
using namespace std;

namespace po = boost::program_options;
#define pi 3.1415926


void Rotate_Point3D(float theta, Eigen::Vector3f& axis, Eigen::Vector3f& input, Eigen::Vector3f& out)
{
    float nx=axis.x();
    float ny=axis.y();
    float nz=axis.z();
    // out[0] =  input[0] * (cosf(theta) + nx * nx * (1 - cosf(theta))) +    //transform by matrix
        // input[1] * (nx * ny * (1 - cosf(theta)) - nz * sinf(theta)) + 
        // input[2] * (nx * nz * (1 - cosf(theta) + ny * sinf(theta)));
    // out[1] = input[0] * (nx * ny * (1 - cosf(theta)) + nz * sinf(theta)) +  
      // input[1] * (ny * ny * (1 - cosf(theta)) + cosf(theta)) + 
      // input[2] * (ny * nz * (1 - cosf(theta)) - nx * sinf(theta));
    // out[2] = input[0] * (nx * nz * (1 - cosf(theta) - ny * sinf(theta))) + 
        // input[1] * (ny * nz * (1 -cosf(theta)) + nx * sinf(theta)) + 
        // input[2] * (nz * nz * (1 - cosf(theta)) + cosf(theta));
    // out=input*cos(theta)+axis.cross(input)*sin(theta)+axis*(axis.dot(input))*(1-cos(theta));
    float cf=cosf(theta);
    float sf=sinf(theta);
    Eigen::Matrix3f rotation_matrix;
    rotation_matrix << cf+nx*nx*(1-cf), nx*ny*(1-cf)-nz*sf, nx*nz*(1-cf)+ny*sf,
                    ny*nx*(1-cf)+nz*sf, cf+ny*ny*(1-cf), ny*nz*(1-cf)-nx*sf,
                    nz*nx*(1-cf)-ny*sf, nz*ny*(1-cf)+nx*sf, cf+nz*nz*(1-cf);
    out=rotation_matrix*input;
}

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

    int n=points1->size();
    pcl::KdTreeFLANN<pcl::PointNormal> kdtree;
    kdtree.setInputCloud(points2);
    int K = 1;
    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);
    std::vector<int> errorlist;
	// for rendering
	for (int i = 0; i < n; i++) {
		pcl::PointNormal& p1=points1->points[i];
        if ( kdtree.nearestKSearch (p1, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0 )
        {
            int idx=pointIdxNKNSearch[0];
            pcl::PointNormal& p2=points2->points[idx];

            Eigen::Vector3f n1=Eigen::Vector3f(p1.normal_x,p1.normal_y,p1.normal_z); 
            Eigen::Vector3f n2=Eigen::Vector3f(p2.normal_x,p2.normal_y,p2.normal_z); 
            if (n1.norm()==0) {
                continue;
            }
            n1.normalize();
            n2.normalize();
            float angle=acos(n1.dot(n2));
            //若n2为空，则跳过
            if (angle!=angle) {
                continue;
            }
            if (angle>pi/2.0) {
                n1*=-1;
                angle=pi-angle;
                p1.normal_x=n1.x(); 
                p1.normal_y=n1.y(); 
                p1.normal_z=n1.z(); 
            }
            // 如果角度大于45,就旋转90度
            if ( angle>pi/4.0 )
            {
                Eigen::Vector3f axis=n2.cross(n1); // 这样就是每次转-90度
                axis.normalize();                
                Eigen::Vector3f out;
                Rotate_Point3D(-pi/2, axis, n1, out);
                out.normalize();
                p1.normal_x=out.x(); 
                p1.normal_y=out.y(); 
                p1.normal_z=out.z(); 
            }
            Eigen::Vector3f out(p1.normal_x,p1.normal_y,p1.normal_z);
            if(acos(out.dot(n2))>=pi/4.0){
                errorlist.push_back(i);
            }
        }
	}

    std::string outfile="consist_"+infile;
    ofstream ofile(outfile.c_str());
    for(int i = 0; i < n; i++) {
        pcl::PointNormal p = points1->points[i];
        ofile << p.x << " " << p.y << " " << p.z << " " << p.normal_x << " " << p.normal_y << " " << p.normal_z << std::endl;
    }
    ofile.close();

    std::string errorfile="elist.txt";
    ofile.open(errorfile.c_str());
    for(int i=0;i<(int)errorlist.size();i++){
        ofile<<errorlist[i]<<endl;
    }
    ofile.close();
    return 0;
}
