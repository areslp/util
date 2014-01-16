#include <iostream>
#include <vector>
#include <fstream>
#include <stack>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

#define pi 3.1415926

using namespace std;

bool is_in_the_same_clazz(MyMesh& mesh, int i, int j); //判断i，j两点是否属于一类
bool is_corner_point(MyMesh& mesh, int j); //判断j点是否是角点
void correct_normal_direction_as_the_first(MyMesh& mesh, int tid,std::vector<int> clazz, std::vector<OpenMesh::Vec3f>& normals); //将clazz里的点法矢调整为tid的方向
void Rotate_Point3D(float theta, OpenMesh::Vec3f& axis, OpenMesh::Vec3f& in, OpenMesh::Vec3f& out) // 将in绕axis旋转theta度，结果保存在out中
{
    float nx=axis[0];
    float ny=axis[1];
    float nz=axis[2];
    out[0] =  in[0] * (cosf(theta) + nx * nx * (1 - cosf(theta))) +    //transform by matrix
        in[1] * (nx * ny * (1 - cosf(theta)) - nz * sinf(theta)) + 
        in[2] * (nx * nz * (1 - cosf(theta) + ny * sinf(theta)));
    out[1] = in[0] * (nx * ny * (1 - cosf(theta)) + nz * sinf(theta)) +  
      in[1] * (ny * ny * (1 - cosf(theta)) + cosf(theta)) + 
      in[2] * (ny * nz * (1 - cosf(theta)) - nx * sinf(theta));
    out[2] = in[0] * (nx * nz * (1 - cosf(theta) - ny * sinf(theta))) + 
        in[1] * (ny * nz * (1 -cosf(theta)) + nx * sinf(theta)) + 
        in[2] * (nz * nz * (1 - cosf(theta)) + cosf(theta));
}


bool is_in_the_same_clazz(MyMesh& mesh, int i, int j){
    // 只有i点相邻面中有一个面与j点相邻面中的法矢非常近似时才认为i，j属于一类
    MyMesh::VertexHandle hi=mesh.vertex_handle(i);
    MyMesh::VertexHandle hj=mesh.vertex_handle(j);

    MyMesh::VertexFaceIter vf_it=mesh.vf_iter(hi);
    for (; vf_it; ++vf_it )
    {
        // 对i点的每一个面
        OpenMesh::Vec3f ni=mesh.normal(vf_it);
        ni.normalize();
        // 遍历j点的每一个面
        MyMesh::VertexFaceIter vf_it2=mesh.vf_iter(hj);
        for (; vf_it2; ++vf_it2 )
        {
            OpenMesh::Vec3f nj=mesh.normal(vf_it2);
            nj.normalize();
            float angle=acos(OpenMesh::dot(ni,nj));
            if ( angle<=5*pi/180 )
            {
                return true;
            }
        }
    }
    return false;
}
bool is_corner_point(MyMesh& mesh, int j){
    MyMesh::VertexHandle hj=mesh.vertex_handle(j);
    MyMesh::VertexIHalfedgeIter half_iter=mesh.vih_iter(hj);
    int count=0; //计算异面角比较大的半边的数量
    for ( ; half_iter; ++half_iter )
    {
        MyMesh::HalfedgeHandle handle = half_iter.handle();
        float angle = mesh.calc_dihedral_angle_fast(handle);
        if(fabs(angle) > 60.0 / 180.0 * pi) {
            ++count;
        }
    }
    if ( count>2 )
    {
        return true;
    }else{
        return false;
    }
}
void correct_normal_direction_as_the_first(MyMesh& mesh, int tid,std::vector<int> clazz, std::vector<OpenMesh::Vec3f>& normals){
    std::vector<int> mark(mesh.n_vertices(),0); //记录已访问过的点
    int clz=clazz[tid];
    std::stack<int> idxs; // 待处理的点堆栈
    idxs.push(tid);
    mark[tid]=1;
    while(true) {
        int pid=idxs.top(); // 得到最顶层元素值
        idxs.pop(); // 删除顶层元素
        OpenMesh::Vec3f nt=normals[pid];
        nt.normalize();
        MyMesh::VertexHandle hv = mesh.vertex_handle(pid);
        // 找相邻点，看相邻点是否为sharp point
        MyMesh::VertexVertexIter vv_it;
        for(vv_it = mesh.vv_iter(hv); vv_it; ++vv_it) {
            MyMesh::VertexHandle hvv = vv_it.handle();
            int idx = hvv.idx();
            if ( clazz[idx]==clz && mark[idx]==0 )
            {
                // 调整该点法矢
                OpenMesh::Vec3f ni=normals[idx];
                ni.normalize();
                float dotvar=OpenMesh::dot(nt,ni);
                if ( dotvar>1 ||dotvar<-1 )
                {
                    dotvar=min(dotvar,(float)1.0);
                    dotvar=max(dotvar,(float)-1.0);
                }
                float angle=acos(dotvar);
                if ( angle>pi/2 )
                {
                    ni*=-1;
                    angle=pi-angle;
                }
                if ( angle>45*pi/180 )
                {
                    // nt作为目标，将ni作为输入，ni绕ntxni旋转-90度，这一段有些问题，还没搞清楚，有时候会出错
                    // OpenMesh::Vec3f axis=OpenMesh::cross(nt,ni);
                    // OpenMesh::Vec3f out;
                    // Rotate_Point3D(-pi/2,axis,ni,out);
                    // out.normalize();
                    // normals[idx]=out;
                    
                    // METHOD2: 遍历idx点的相邻面，在里面找一个与nt足够相近的面法矢作为ni
                    MyMesh::VertexFaceIter vf_it;
                    for ( vf_it=mesh.vf_iter(mesh.vertex_handle(pid));vf_it ; ++vf_it )
                    {
                        OpenMesh::Vec3f nif=mesh.normal(vf_it);
                        nif.normalize();
                        dotvar=OpenMesh::dot(nt,nif);
                        if ( dotvar>1 ||dotvar<-1 )
                        {
                            dotvar=min(dotvar,(float)1.0);
                            dotvar=max(dotvar,(float)-1.0);
                        }
                        angle=acos(dotvar);
                        if ( angle<10*pi/180 )
                        {
                            normals[idx]=nif;
                        }
                    }
                }
                mark[idx]=1;
                idxs.push(idx);
            }
        }
        // 没有待处理元素则跳出while
        if ( idxs.size()==0 )
        {
            break;
        }
    }
}
int main(int argc, char **argv)
{
    MyMesh  mesh;

    // check command line options
    if(argc != 2) {
        std::cerr << "Usage:  " << argv[0] << "save_normal infile\n";
        return 1;
    }

    // read mesh from stdin
    if(! OpenMesh::IO::read_mesh(mesh, argv[1])) {
        std::cerr << "Error: Cannot read mesh from " << argv[1] << std::endl;
        return 1;
    }

    cout << "input point size:" << mesh.n_vertices() << endl;
    cout << "updating face normal..." << endl;
    mesh.request_face_normals();
    mesh.update_face_normals();
    cout << "updated face normal" << endl;
    vector<int> sharp_vertex(mesh.n_vertices(), 0);

    // 找到sharp points,并保存到sharp_vertex里
    for(int i = 0; i < (int)mesh.n_halfedges(); ++i) {
        MyMesh::HalfedgeHandle handle = mesh.halfedge_handle(i);
        float angle = mesh.calc_dihedral_angle_fast(handle);

        if(fabs(angle) > 25.0 / 180.0 * pi) {
            MyMesh::VertexHandle v1 = mesh.from_vertex_handle(handle);
            MyMesh::VertexHandle v2 = mesh.to_vertex_handle(handle);
            sharp_vertex[v1.idx()] = 1;
            sharp_vertex[v2.idx()] = 1;
        }
    }

    // 保存下sharp points，在meshlab中可视化看是否正确
    ofstream spf("sharp_points.txt");
    for(int i = 0; i < (int)sharp_vertex.size(); ++i) {
        spf << i << " " << sharp_vertex[i] << endl;
    }
    spf.close();

    // DEBUG 将角点保存出去可视化
    ofstream cpf("corner_points.txt");
    for(int i = 0; i < (int)mesh.n_vertices(); ++i) {
        if ( is_corner_point(mesh,i) )
        {
            cpf<<i<<" "<<1<<endl;
        }
    }
    cpf.close();


    // 将sharp points分解为几个连通的子分组
    int clazz_idx = 0; //保存当前的类型id
    std::vector<int> clazz(mesh.n_vertices(), 0); // 记录sharp点的类型

    for(int i = 0; i < (int)mesh.n_vertices(); ++i) {
        if(sharp_vertex[i] == 0 || clazz[i]!=0 || is_corner_point(mesh,i)) {
            continue;
        }
        // 找到了一个非角点的初始点
        ++clazz_idx;
        clazz[i]=clazz_idx;
        std::stack<int> idxs; // 待处理的点堆栈
        idxs.push(i);
        while(true) {
            int pid=idxs.top(); // 得到最顶层元素值
            idxs.pop(); // 删除顶层元素
            MyMesh::VertexHandle hv = mesh.vertex_handle(pid);
            // 找相邻点，看相邻点是否为sharp point
            MyMesh::VertexVertexIter vv_it;
            // 如果邻域包含角点，则所有点都只标记，而不加入堆栈
            bool hasCorner=false;
            std::vector<int> ips; //待加入堆栈的点
            for(vv_it = mesh.vv_iter(hv); vv_it; ++vv_it) {
                MyMesh::VertexHandle hvv = vv_it.handle();
                int idx = hvv.idx();
                assert(idx!=pid);
                // 判断：1.当前点是否已经被处理过
                // 2. 当前点是否为sharp points中的点
                // 3. 当前点的法矢与i点法矢差值是否在0,90,180附近
                if(sharp_vertex[idx] == 1 && clazz[idx]==0) {
                    if ( is_in_the_same_clazz(mesh,pid,idx))
                    {
                        clazz[idx]=clazz_idx;
                        // 如果新点是角点就不继续找了
                        if ( !is_corner_point(mesh,idx) )
                        {
                            ips.push_back(idx);
                        }else{
                            hasCorner=true;
                        }
                    }
                }
            }
            if ( !hasCorner )
            {
                for ( std::vector<int>::iterator it = ips.begin(); it!=ips.end(); ++it)
                {
                    idxs.push(*it);
                }
            }
            // 没有待处理元素则跳出while
            if ( idxs.size()==0 )
            {
                break;
            }
        }
    }

    // 保存分组可视化结果
    ofstream visf("vis.xyznq");
    for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        OpenMesh::Vec3f p = mesh.point(v_it);
        int idx=v_it.handle().idx();
        int q=clazz[idx];
        visf << p[0] << " " << p[1] << " " << p[2] << " " << 0 << " " << 0 << " " << 0 << " " << q << endl;
    }
    visf.close();

    // 单独保存一份分组结果
    ofstream gf("group.txt");
    for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        int idx=v_it.handle().idx();
        int q=clazz[idx];
        gf << q << endl;
    }
    gf.close();
    
    cout <<"Finish write file, face normal to vertex normal"<<endl;
    vector<OpenMesh::Vec3f> normals;
    // 给每一个点的法矢先赋一个值
    for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        // do something with *v_it, v_it->, or v_it.handle()
        MyMesh::VertexFaceIter  vf_it;
        for(vf_it = mesh.vf_iter(v_it); vf_it; ++vf_it) {
            //cout<<"iterating faces around the vertex"<<endl;
            OpenMesh::Vec3f n;
            n = mesh.normal(vf_it);
            // cout<<"normal:"<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
            if(n[0] == 0 && n[1] == 0 && n[2] == 0) {
                cout << "Warning: NULL normal!" << endl;
                cout << n[0] << " " << n[1] << " " << n[2] << endl;
            }
            n.normalize();
            normals.push_back(n);
            break;
        }
    }

    // TEST Rotate_Point3D
    // OpenMesh::Vec3f tn1=normals[27];
    // OpenMesh::Vec3f tn2=normals[28];
    // OpenMesh::Vec3f axis=OpenMesh::cross(tn2,tn1);
    // axis.normalize();
    // cout<<"axis"<<axis[0]<<" "<<axis[1]<<" "<<axis[2]<<endl;
    // OpenMesh::Vec3f tout;
    // Rotate_Point3D(-pi/2,axis,tn1,tout);
    // tout.normalize();
    // normals[27]=tout;
    // cout<<"dot:"<<OpenMesh::dot(tout,axis)<<endl;
    // cout<<"dot:"<<OpenMesh::dot(tn1,axis)<<endl;
    // cout<<"tout:"<<tout[0]<<" "<<tout[1]<<" "<<tout[2]<<endl;

    // cout<<"write face2vertex file"<<endl;
    // ofstream fvf("face2vertex.xyzn");
    // for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        // OpenMesh::Vec3f p = mesh.point(v_it);
        // OpenMesh::Vec3f n = normals[v_it.handle().idx()];
        // fvf<< p[0] << " " << p[1] << " " << p[2] << " " << n[0] << " " << n[1] << " " << n[2] << endl;
    // }
    // fvf.close();


    cout <<"Rotate vertex normal"<<endl;
    // 遍历clazz，调整法矢方向, clazz_idx为总的分类
    for ( int i = 1; i <= clazz_idx; ++i ) // i为一类
    {
        for ( int j = 0; j < (int)clazz.size(); ++j ) //在所有点中找
        {
            if ( clazz[j]==i && ! is_corner_point(mesh,j) )
            {
                correct_normal_direction_as_the_first(mesh,j,clazz,normals);
                break;
            }
        }
    }

    cout << "output" << endl;
    ofstream outfile("standard_normal.xyzn");
    int index = 0;

    for(MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        OpenMesh::Vec3f p = mesh.point(v_it);
        OpenMesh::Vec3f n = normals[index];
        outfile << p[0] << " " << p[1] << " " << p[2] << " " << n[0] << " " << n[1] << " " << n[2] << endl;
        ++index;
    }

    outfile.close();
    return 0;
}
