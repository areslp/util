#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <boost/lexical_cast.hpp>

#include <utility> // defines std::pair
#include <list>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;

int main(int argc,char** argv)
{
	int k;
	std::string fname;
	std::string oname;
	if (argc==4)
	{
		fname=argv[1];
		k=boost::lexical_cast<int>(argv[2]);
		oname=argv[3];
	}else{
		std::cout<<"argc:"<<argc<<",exiting..."<<std::endl;
		std::cout<<"normal_estimation infile k outfile"<<std::endl;
		return 0;
	}
	// Reads a .xyz point set file in points[].
	std::list<PointVectorPair> points;
	std::ifstream stream(fname.c_str());
	if (!stream ||
		!CGAL::read_xyz_points(stream,
		std::back_inserter(points),
		CGAL::First_of_pair_property_map<PointVectorPair>()))
	{
		std::cerr << "Error: cannot read data file" << std::endl;
		return EXIT_FAILURE;
	}

	// Estimates normals direction.
	// Note: pca_estimate_normals() requires an iterator over points
	// as well as property maps to access each point's position and normal.
	const int nb_neighbors = k; // K-nearest neighbors = 3 rings
	CGAL::pca_estimate_normals(points.begin(), points.end(),
		CGAL::First_of_pair_property_map<PointVectorPair>(),
		CGAL::Second_of_pair_property_map<PointVectorPair>(),
		nb_neighbors);


	std::ofstream pca("pca.xyzn");
	for (std::list<PointVectorPair>::iterator it=points.begin();it!=points.end();it++)
	{
		PointVectorPair pp=*it;
		Point point=pp.first;
		Vector normal=pp.second;
		pca<<point.x()<<" "<<point.y()<<" "<<point.z()<<" "<<normal.x()<<" "<<normal.y()<<" "<<normal.z()<<std::endl;
	}
	pca.close();

	// Orients normals.
	// Note: mst_orient_normals() requires an iterator over points
	// as well as property maps to access each point's position and normal.
	std::list<PointVectorPair>::iterator unoriented_points_begin =
		CGAL::mst_orient_normals(points.begin(), points.end(),
		CGAL::First_of_pair_property_map<PointVectorPair>(),
		CGAL::Second_of_pair_property_map<PointVectorPair>(),
		6);

	// Optional: delete points with an unoriented normal
	// if you plan to call a reconstruction algorithm that expects oriented normals.
	points.erase(unoriented_points_begin, points.end());

	// Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
	std::list<PointVectorPair>(points).swap(points);

	std::ofstream ot(oname.c_str());
	for (std::list<PointVectorPair>::iterator it=points.begin();it!=points.end();it++)
	{
		PointVectorPair pp=*it;
		Point point=pp.first;
		Vector normal=pp.second;
		ot<<point.x()<<" "<<point.y()<<" "<<point.z()<<" "<<normal.x()<<" "<<normal.y()<<" "<<normal.z()<<std::endl;
	}
	ot.close();

	return EXIT_SUCCESS;
}
