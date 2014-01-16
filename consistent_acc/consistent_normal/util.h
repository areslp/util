/*
 * =====================================================================================
 *
 *       Filename:  util.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2011年12月13日 19时30分15秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  util_INC
#define  util_INC
/*
 * =====================================================================================
 *        Class:  Util
 *  Description:  
 * =====================================================================================
 */

#include	<boost/shared_ptr.hpp>
#include <pcl/common/common_headers.h>
#include <pcl/features/normal_3d.h>
#include <sstream>
#include <iomanip>
#include <Eigen/Dense>

#define PI 3.1415926

class Util
{
    public:

        /* ====================  LIFECYCLE     ======================================= */
        Util(){};
        ~Util (){};                            /* destructor       */

        /* ====================  ACCESSORS     ======================================= */

        /* ====================  MUTATORS      ======================================= */

        /* ====================  OPERATORS     ======================================= */
        static float angle_between_vectors(Eigen::Vector3f& v1, Eigen::Vector3f& v2){
            v1.normalize();
            v2.normalize();
            float dotv = v1.dot(v2);
            //fix numerical problem
            if(dotv > 1) {
                dotv = 1;
            } else if(dotv < -1) {
                dotv = -1;
            }
            float angle = acos(dotv);
            assert(angle >= 0 && angle <= PI + 0.001);
            if(angle > PI / 2) {
                angle = PI - angle;
            }
            return angle;
        }

        static pcl::PointCloud<pcl::PointXYZ>::Ptr subPointCloud(pcl::PointCloud<pcl::PointNormal>::Ptr points, std::vector<int> subIdx)
        {
            pcl::PointCloud<pcl::PointXYZ>::Ptr subPointCloud(new pcl::PointCloud<pcl::PointXYZ>());
            for(size_t i = 0; i < subIdx.size(); ++i) {
                int idx = subIdx.at(i);
                pcl::PointNormal p = points->points[idx];
                subPointCloud->push_back(pcl::PointXYZ(p.x, p.y, p.z));
            }
            return subPointCloud;
        }
     
        static void read_points(std::string infile, pcl::PointCloud<pcl::PointNormal>::Ptr points)
        {
            std::ifstream in(infile.c_str());
            std::string tmp;
            int line_count = 0;

            while(std::getline(in, tmp)) {
                line_count++;
            }
            in.clear();
            in.seekg(0);
            points->resize(line_count);
            int i = 0;
            double x, y, z, nx, ny, nz;
            while(in >> x >> y >> z >> nx >> ny >> nz) {
                pcl::PointNormal& p = points->points[i];
                p.x = x;
                p.y = y;
                p.z = z;
                p.normal_x = 0;
                p.normal_y = 0;
                p.normal_z = 0;
                i++;
            }
            in.close();
        }
        static void read_pointnormals(std::string infile, pcl::PointCloud<pcl::PointNormal>::Ptr points)
        {
            std::ifstream in(infile.c_str());
            std::string tmp;
            int line_count = 0;
            while(std::getline(in, tmp)) {
                line_count++;
            }
            in.clear();
            in.seekg(0);
            points->resize(line_count);
            int i = 0;
            double x, y, z, nx, ny, nz;
            while(in >> x >> y >> z >> nx >> ny >> nz) {
                pcl::PointNormal& p = points->points[i];
                p.x = x;
                p.y = y;
                p.z = z;
                p.normal_x = nx;
                p.normal_y = ny;
                p.normal_z = nz;
                i++;
            }
            in.close();
        }

        Util& operator = ( const Util &other ); /* assignment operator */

    protected:
        /* ====================  DATA MEMBERS  ======================================= */

    private:

}; /* -----  end of class Util  ----- */
#endif   /* ----- #ifndef util_INC  ----- */

