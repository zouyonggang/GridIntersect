//
// 文件名: main.C
// 软件包: JAUMIN application
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 136 $
// 修改  : $Date: 2011-07-22 08:41:30 +0800 (五, 2011-07-22) $
// 描述  : 主控程序(矩形非结构网格上, 求解线性对流对流问题).
//

#include <string>
#include <time.h>
#include <vector>
using namespace std;

#include "ArrayData.h"
#include "CompareData.h"
#include "EntityUtilities.h"
#include "GridGeometry.h"
#include "GridTopology.h"
#include "HierarchyTimeIntegrator.h"
#include "InputManager.h"
#include "JAUMINManager.h"
#include "JAUMIN_config.h"
#include "JaVisDataWriter.h"
#include "PatchHierarchy.h"
#include "RestartManager.h"
#include "TimerManager.h"
#include "Utilities.h"
#include "VariableDatabase.h"
using namespace JAUMIN;

#include "GridIntersect.h"
#include "LinAdv.h"
#include "LinAdvLevelIntegrator.h"

/*!
*************************************************************************
* @brief   在网格文件名前面追加input文件的路径
************************************************************************
*/
static void prefixInputDirName(const string& input_filename,
                               tbox::Pointer<tbox::Database> input_db) {
  string path_name = "";
  string::size_type slash_pos = input_filename.find_last_of('/');
  if (slash_pos != string::npos)
    path_name = input_filename.substr(0, slash_pos + 1);

  string mesh_file = input_db->getDatabase("GridGeometry")
                         ->getDatabase("MeshImportationParameter")
                         ->getString("file_name");

  slash_pos = mesh_file.find_first_of('/');
  if (slash_pos != 0) {
    input_db->getDatabase("GridGeometry")
        ->getDatabase("MeshImportationParameter")
        ->putString("file_name", path_name + mesh_file);
  }
}

/*!
*************************************************************************
*
* @brief 基于JAUMIN框架的非结构网格, 求解线性对流问题.
*
* 该程序分以下几个步骤:
* -# 预处理: 初始化MPI和JAUMIN环境, 解析输入文件, 读取主程序控制参数;
* -# 创建网格片层次结构和时间积分算法类对象, 主要包括:
*       - 网格几何       hier::GridGeometry<NDIM>
*       - 网格拓扑       hier::GridTopology<NDIM>
*       - 网格片层次结构 hier::PatchHierarchy<NDIM>
*    -# 网格片积分算法 LinAdv
*       - 应用级: 提供求解线性对流方程的数值计算子程序
*    -# 网格层积分算法 LinAdvLevelIntegrator
*       - 应用级: 提供基于网格层的线性对流方程求解流程
*    -# 网格层时间积分算法 algs::HierarchyTimeIntegrator<NDIM>
* -# 初始化网格片层次结构和物理量数据片;
* -# 循环: 时间步进;
* -# 后处理: 释放应用类对象, 释放JAUMIN和MPI内部资源.
*
************************************************************************
*/
int main(int argc, char* argv[]) {
  struct timespec start, end;
  // 初始化MPI和JAUMIN环境.
  tbox::MPI::init(&argc, &argv);
  tbox::JAUMINManager::startup();

  {
    /*******************************************************************************
     *                               预  处  理 *
     *******************************************************************************/
    // 解析命令行参数:
    string input_filename;  // 输入文件名.
    if (argc != 2) {
      tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
                 << "<restart dir> <restore number> [options]\n"
                 << "  options:\n"
                 << "  none at this time" << endl;
      tbox::MPI::abort();
      return (-1);
    } else {
      input_filename = argv[1];
    }

    /// 把信息输出到log文件
    tbox::plog << "input_filename = " << input_filename << endl;

    // 解析输入文件的计算参数到输入数据库, 称之为根数据库.
    tbox::Pointer<tbox::Database> input_db =
        new tbox::InputDatabase("input_db");
    tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

    // 在网格文件名前面追加input文件的路径
    prefixInputDirName(input_filename, input_db);

    // 从根数据库中获得名称为"Main"的子数据库.
    tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

    /*******************************************************************************
     *                     创建网格片层次结构和积分算法类对象 *
     *******************************************************************************/
    //(1) 创建非结构网格几何.
    tbox::Pointer<hier::GridGeometry<NDIM> > grid_geometry =
        new hier::GridGeometry<NDIM>("GridGeometry",
                                     input_db->getDatabase("GridGeometry"));

    //(2) 创建非结构网格拓扑.
    tbox::Pointer<hier::GridTopology<NDIM> > grid_topology =
        new hier::GridTopology<NDIM>("GridTopology",
                                     input_db->getDatabase("GridTopology"));

    //(3) 创建网格片层次结构.
    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
        new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry,
                                       grid_topology, true);

    //(4) 创建网格片时间积分算法类（应用级:
    //提供求解线性对流方程的数值计算子程序）.
    LinAdv* linadv_advection_model =
        new LinAdv("LinAdv", input_db->getDatabase("LinAdv"));

    //(5) 创建网格层时间积分算法类（应用级:
    //提供基于网格层的线性对流方程求解流程）.
    LinAdvLevelIntegrator* linadv_level_integrator = new LinAdvLevelIntegrator(
        "LinAdvLevelIntegrator", linadv_advection_model);

    //(6) 创建网格层时间积分算法类.
    tbox::Pointer<algs::HierarchyTimeIntegrator<NDIM> > time_integrator =
        new algs::HierarchyTimeIntegrator<NDIM>(
            "HierarchyTimeIntegrator",
            input_db->getDatabase("HierarchyTimeIntegrator"), patch_hierarchy,
            linadv_level_integrator, true);

    // 初始化网格片层次结构和物理量.
    time_integrator->initializeHierarchy();

    /************************************************************************************
     *                              网 格 相 交 测 试 数 据 *
     ************************************************************************************/
    {
      //获取网格层
      int number_levels = patch_hierarchy->getNumberOfLevels();
      tbox::pout << endl << "网格层总数：" << number_levels << endl;
      for (int level_id = 0; level_id < number_levels; level_id++) {
        tbox::pout << "当前网格层：" << level_id << endl;
        tbox::Pointer<hier::PatchLevel<NDIM> > patch_level =
            patch_hierarchy->getPatchLevel(level_id);
        tbox::Pointer<hier::Patch<NDIM> > source_patch;
        for (hier::PatchLevel<NDIM>::Iterator p(patch_level); p; p++) {
          tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++"
                     << endl;
          tbox::Pointer<hier::Patch<NDIM> > patch = patch_level->getPatch(p());
          if (patch->getIndex() == 0) source_patch = patch;
          tbox::pout << "patch index:" << patch->getIndex() << endl;

          int cell_number =
              patch->getNumberOfEntities(hier::EntityUtilities::CELL, 0);
          //   int cell_number_with_ghost =
          //       patch->getNumberOfEntities(hier::EntityUtilities::CELL, 1);
          tbox::pout << "Patch cell entity:" << cell_number << endl;
          int node_number =
              patch->getNumberOfEntities(hier::EntityUtilities::NODE, 0);
          //   int node_number_with_ghost =
          //       patch->getNumberOfEntities(hier::EntityUtilities::NODE, 1);
          tbox::Pointer<hier::PatchGeometry<NDIM> > geometry =
              patch->getPatchGeometry();
          tbox::Pointer<pdat::CellData<NDIM, double> > cell_coordinates =
              geometry->getCellCoordinates();
          double* cell_coordinate = cell_coordinates->getPointer();
          // double cell_coordinate[3] = {0, 0, 0};
          tbox::Pointer<pdat::NodeData<NDIM, double> > node_coordinates =
              geometry->getNodeCoordinates();
          const double* node_coordinate = node_coordinates->getPointer(0);

          //将x轴最下方的顶点赋值给boundary-point
          double x_lower_point[3] = {
              (*cell_coordinates)(0, 0, 0),
              (*cell_coordinates)(1, 0, 0),
              (*cell_coordinates)(2, 0, 0),
          };
          for (int i = 0; i < node_number; i++)
            if (node_coordinate[i * NDIM] > x_lower_point[0])
              for (int j = 0; j < NDIM; j++)
                x_lower_point[j] = node_coordinate[i * NDIM + j];

          //进行点与网格定位
          clock_gettime(CLOCK_MONOTONIC_RAW, &start);
          tbox::Pointer<GridIntersect> intersect = new GridIntersect(patch);
          clock_gettime(CLOCK_MONOTONIC_RAW, &end);
          cout << "构造函数: "
               << end.tv_sec + end.tv_nsec * 1e-9 - start.tv_sec -
                      start.tv_nsec * 1e-9
               << " s\n";

          //取网格单元坐标作为内部输入测试点
          clock_gettime(CLOCK_MONOTONIC_RAW, &start);
          vector<int> inside_result =
              intersect->pointInGrid(cell_coordinate, 10000);
          clock_gettime(CLOCK_MONOTONIC_RAW, &end);
          cout << "\n点定位内部(10000): "
               << end.tv_sec + end.tv_nsec * 1e-9 - start.tv_sec -
                      start.tv_nsec * 1e-9
               << " s\n";

          //取x方向上最大顶点作为边界输入测试点
          clock_gettime(CLOCK_MONOTONIC_RAW, &start);
          vector<int> boundary_result =
              intersect->pointInGrid(x_lower_point, 1);
          clock_gettime(CLOCK_MONOTONIC_RAW, &end);
          cout << "\n点定位边界(1): "
               << end.tv_sec + end.tv_nsec * 1e-9 - start.tv_sec -
                      start.tv_nsec * 1e-9
               << " s\n";

          // 取影像区网格单元中心坐标作为外部测试点输入
          //   double outside_point[3] = {x_lower_point[0] + 1,
          //   x_lower_point[1],
          //                              x_lower_point[2]};
          vector<int> outside_result;
          clock_gettime(CLOCK_MONOTONIC_RAW, &start);
          outside_result =
              intersect->pointInGrid(cell_coordinate + cell_number * 3, 100);
          clock_gettime(CLOCK_MONOTONIC_RAW, &end);
          cout << "\n点定位外部(100): "
               << end.tv_sec + end.tv_nsec * 1e-9 - start.tv_sec -
                      start.tv_nsec * 1e-9
               << " s\n";

          //进行射线网格相交测试,取上诉外部点作为起点，x轴负方向作为方向
          //   double points[3] = {outside_point[0], outside_point[1],
          //                       outside_point[2]};
          //   double direction[3] = {-1, 0, 0};
          //   double intersection_coordinates[2 * 3];
          double points[100000 * 3] = {0};
          double direction[100000 * 3];
          double intersection_coordinates[100000 * 3];
          for (int i = 0; i < 300000; i++) {
            points[i] = 0;
            direction[i] = i - 100000;
          }

          std::vector<int> ray_intersect_result;
          clock_gettime(CLOCK_MONOTONIC_RAW, &start);
          intersect->rayIntersectGrid(points, direction, 100,
                                      ray_intersect_result,
                                      intersection_coordinates);
          clock_gettime(CLOCK_MONOTONIC_RAW, &end);
          cout << "\n射线相交(100): "
               << end.tv_sec + end.tv_nsec * 1e-9 - start.tv_sec -
                      start.tv_nsec * 1e-9
               << " s\n";

          //
          // 网格相交测试
          //
          std::vector<int> grid_intersect_result;
          int intersect_number;
          clock_gettime(CLOCK_MONOTONIC_RAW, &start);
          intersect->gridIntersectGrid(source_patch, intersect_number,
                                       grid_intersect_result);
          clock_gettime(CLOCK_MONOTONIC_RAW, &end);
          cout << "\n网格相交: "
               << end.tv_sec + end.tv_nsec * 1e-9 - start.tv_sec -
                      start.tv_nsec * 1e-9
               << " s\n";

          //
          // 输出结果
          //
          tbox::pout << "inside相交网格编号为:";
          for (int i = 0; i < static_cast<int>(inside_result.size()); i++)
            tbox::pout << inside_result[i] << " ";
          tbox::pout << endl << "boundary相交网格编号为:";
          for (int i = 0; i < static_cast<int>(boundary_result.size()); i++)
            tbox::pout << boundary_result[i] << " ";
          tbox::pout << endl << "outside相交网格编号为:";
          for (int i = 0; i < static_cast<int>(outside_result.size()); i++)
            tbox::pout << outside_result[i] << " ";
          tbox::pout << endl << "射线相交:" << endl;
          for (int i = 0; i < static_cast<int>(ray_intersect_result.size());
               i++) {
            tbox::pout << " 交点网格单元 " << ray_intersect_result[i]
                       << " , 交点坐标 ";
            tbox::pout << intersection_coordinates[i * 3] << " "
                       << intersection_coordinates[i * 3 + 1] << " "
                       << intersection_coordinates[i * 3 + 2];
            tbox::pout << endl;
          }
          tbox::pout << endl << "网格相交:" << endl;
          tbox::pout << "相交单元数:" << intersect_number << endl
                     << "交点索引号：";
          for (int i = 0; i < static_cast<int>(grid_intersect_result.size());
               i++) {
            tbox::pout << grid_intersect_result[i] << ",";
          }
          tbox::pout << endl;
        }  // end patch
      }    // end level
    }      // end test
  }

  // 释放JAUMIN和MPI内部资源.
  tbox::JAUMINManager::shutdown();
  tbox::MPI::finalize();

  return (0);
}
