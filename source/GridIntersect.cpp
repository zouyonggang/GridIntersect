#ifndef GRID_INTERSECT_CPP
#define GRID_INTERSECT_CPP

#include "GridIntersect.h"

#include "Array.h"
#include "BoundingBox.h"
#include "DoubleVector.h"
#include "IntervalTree.h"
#include "NodeData.h"
#include "PatchGeometry.h"
#include "PatchTopology.h"

#include <assert.h>
#include <map>
#include <string>
#include <unordered_map>

#include <iostream>
using namespace std;

/**
 * @brief 描述一个3维结点，里面包含patch内节点编号，节点坐标
 *
 */
struct NodeDescriptor {
  int index;
  hier::DoubleVector<3> coord;
};

/**
 * @brief 描述一个面，包含patch内面编号，面周围节点编号、棱编号
 *
 */
struct FaceDescriptor {
  int index;
  std::vector<int> nodes;
  std::vector<int> cells;  // 用来获取网格外表面
};

/**
 * @brief 描述一个三维体，包含patch体编号，体周围节点编号、棱编号、面编号
 *
 */
struct CellDescriptor {
  int index;
  std::vector<int> nodes;
  std::vector<int> faces;
};

class GridIntersectImpl {
public:
  /**
   * @brief构造函数，从jaumin框架内得到点、边、面、体（非影像区），并保存到该类结构中
   *
   * @param source_patch 输入的patch网格片
   */
  GridIntersectImpl(tbox::Pointer<hier::Patch<3> >& source_patch) {
    delta_ = 1.0e-20;
    int node_number =
        source_patch->getNumberOfEntities(hier::EntityUtilities::NODE, 0);
    int face_number =
        source_patch->getNumberOfEntities(hier::EntityUtilities::FACE, 0);
    int cell_number =
        source_patch->getNumberOfEntities(hier::EntityUtilities::CELL, 0);
    tbox::Pointer<hier::PatchTopology<3> > topology =
        source_patch->getPatchTopology();

    // 获取节点,保存节点
    tbox::Pointer<hier::PatchGeometry<3> > patch_geometry =
        source_patch->getPatchGeometry();
    tbox::Pointer<pdat::NodeData<3, double> > nodes =
        patch_geometry->getNodeCoordinates();
#ifdef DEBUG_CHECK_ASSERTION
    assert(nodes->getDepth() == 1);
    assert(nodes->getGroup() == 3);
#endif
    const double* nodes_coordinates = nodes->getPointer(0);
    for (int i = 0; i < node_number; i++) {
      NodeDescriptor* node = new NodeDescriptor;
      node->index = i;
      node->coord = hier::DoubleVector<3>(nodes_coordinates[i * 3],
                                          nodes_coordinates[i * 3 + 1],
                                          nodes_coordinates[i * 3 + 2]);
      patch_nodes_[i] = node;
    }

    // 获取面
    tbox::Array<int> face_adj_nodes_extent;
    tbox::Array<int> face_adj_nodes_indices;
    topology->getFaceAdjacencyNodes(face_adj_nodes_extent,
                                    face_adj_nodes_indices);
    tbox::Array<int> face_adj_cells_extent;
    tbox::Array<int> face_adj_cells_indices;
    topology->getFaceAdjacencyCells(face_adj_cells_extent,
                                    face_adj_cells_indices);
    for (int i = 0; i < face_number; i++) {
      FaceDescriptor face;
      face.index = i;

      // 获取到的面点表中，有可能出现点重复的情况，所有要排除重复的点
      for (int j = face_adj_nodes_extent[i]; j < face_adj_nodes_extent[i + 1];
           j++) {
        std::vector<int>::iterator iter = std::find(
            face.nodes.begin(), face.nodes.end(), face_adj_nodes_indices[j]);
        if (iter == face.nodes.end()) {
          face.nodes.push_back(face_adj_nodes_indices[j]);
        }
      }
      // 如果组成该面的点小于3则抛弃该点
      if (static_cast<int>(face.nodes.size()) < 3) continue;

      // 获取面周围的体，并判断该面是否为外表面
      for (int j = face_adj_cells_extent[i]; j < face_adj_cells_extent[i + 1];
           j++) {
        if (face_adj_cells_extent[i + 1] - face_adj_cells_extent[i] == 1)
          patch_outer_faces_[i] = face_adj_cells_indices[j];
        face.cells.push_back(face_adj_cells_indices[j]);
      }
      patch_faces_[i] = face;
    }

    // 获取体
    tbox::Array<int> cell_adj_nodes_extent;
    tbox::Array<int> cell_adj_nodes_indices;
    topology->getCellAdjacencyNodes(cell_adj_nodes_extent,
                                    cell_adj_nodes_indices);
    tbox::Array<int> cell_adj_faces_extent;
    tbox::Array<int> cell_adj_faces_indices;
    topology->getCellAdjacencyFaces(cell_adj_faces_extent,
                                    cell_adj_faces_indices);
    for (int i = 0; i < cell_number; i++) {
      CellDescriptor cell;
      cell.index = i;
      for (int j = cell_adj_nodes_extent[i]; j < cell_adj_nodes_extent[i + 1];
           j++)
        cell.nodes.push_back(cell_adj_nodes_indices[j]);
      for (int j = cell_adj_faces_extent[i]; j < cell_adj_faces_extent[i + 1];
           j++) {
        // 如果该面在patch_faces_中不存在（构成该面的点小于3），则不添加该面
        std::unordered_map<int, FaceDescriptor>::iterator iter =
            patch_faces_.find(cell_adj_faces_indices[j]);
        if (iter != patch_faces_.end())
          cell.faces.push_back(cell_adj_faces_indices[j]);
      }

      patch_cells_[i] = cell;
    }

    // 构造区间树，将patch的所有网格单元bbox存入区间书的叶子节点
    interval_tree_ = new hier::IntervalTree<3>(patch_cells_.size());
    for (std::unordered_map<int, CellDescriptor>::iterator it =
             patch_cells_.begin();
         it != patch_cells_.end(); it++) {
      CellDescriptor& cell = it->second;
      // 区间树输入要求列主序
      double bbox[6] = {1.0e16, -1.0e16, 1.0e16, -1.0e16, 1.0e16, -1.0e16};
      for (int i = 0; i < static_cast<int>(cell.nodes.size()); i++) {
#ifdef DEBUG_CHECK_ASSERTION
        assert(cell.nodes[i] <= patch_nodes_.size());
#endif
        NodeDescriptor* node = patch_nodes_[cell.nodes[i]];
        for (int d = 0; d < 3; d++) {
          bbox[d * 2] = std::min(bbox[d * 2], node->coord[d]);
          bbox[d * 2 + 1] = std::max(bbox[d * 2 + 1], node->coord[d]);
        }
      }
      interval_tree_->addElement(it->first, bbox);
    }
    interval_tree_->constructTree();
  }

  ~GridIntersectImpl() {
    for (std::unordered_map<int, NodeDescriptor*>::iterator iter =
             patch_nodes_.begin();
         iter != patch_nodes_.end(); iter++)
      delete iter->second;
  }

  /**
   * @brief 判断点2，是否在点0，点1组成的线段的左边，只能适用于2维情况
   *
   * @param point0 点0，长度为2
   * @param point1 点1，长度为2
   * @param point2 点2，长度为2
   * @return int >0 p2 在通过p0,p1线段的左边
   *             =0 p2 在通过p0,p1线段上
   *             <0 p2 在通过p0,p1线段的右边
   */
  inline double isLeft(const double* point0, const double* point1,
                       const double* point2) {
#ifdef DEBUG_CHECK_ASSERTION
    assert(point1 != NULL && point2 != NULL && point0 != NULL);
#endif
    return ((point1[0] - point0[0]) * (point2[1] - point0[1]) -
            (point2[0] - point0[0]) * (point1[1] - point0[1]));
  }

  /**
   * @brief 利用缠绕算法，计算二维平面内，点是否在多边形内
   *
   * @param point 输入点
   * @param vertex 输入多边形的顶点
   * @param n 多边形顶点的个数
   * @return int !=0 点在多边形内部
   *             =0 点在多边形外部
   *
   * @note vertex要保证最后一个V[n] = V[0]
   */
  int wnPnpoly(const double* point, const double* vertex, int n) {
#ifdef DEBUG_CHECK_ASSERTION
    assert(point != NULL && vertex != NULL);
    assert(n > 2);
#endif

    int wn = 0;  // 记录winding number 数量

    // 遍历多边形的所有边
    for (int i = 0; i < n; i++) {  //  vertex[i]和vertex[i+1]组成的边
      const double V[2] = {vertex[i * 2], vertex[i * 2 + 1]};
      const double V_next[2] = {vertex[(i + 1) * 2], vertex[(i + 1) * 2 + 1]};

      if (isLeft(V, V_next, point) == 0 &&
          ((point[0] >= V[0] && point[0] <= V_next[0]) ||
           (point[0] >= V_next[0] && point[0] <= V[0])))
        ++wn;  // 如果该点在边上
      else {
        if ((V[1] - point[1]) <= 1.0e-8) {  // 输入点在线段下顶点上方
          if ((V_next[1] - point[1]) > delta_)      // 一个向上的缠绕
            if (isLeft(V, V_next, point) > delta_)  //  p在边的左边
              ++wn;
        } else {  // 输入点在线段下顶点上方
          if ((V_next[1] - point[1]) <= 1.0e-8)     // 一个向下的缠绕
            if (isLeft(V, V_next, point) < delta_)  //  p在边的右边
              --wn;
        }
      }
    }

    return wn;
  }

  /**
   * @brief 计算射线与平面的交点
   *
   * @param point 射线起点
   * @param direction 射线方向
   * @param face 平面（由三个点构成）行主序（存完第一个点，再存第二个点）
   * @param intersection 输出参数，交点坐标
   * @param dist 输出参数，射线起点到交点距离
   * @param tag 输出参数，标识状态，
   *              0：成功求出交点
   *             -1：射线与平面平行
   *             -2：射线所在直线与平面的交点在射线的副方向
   */
  void lineSurfaceIntersect(const double* point, const double* direction,
                            const double* face, double* intersection,
                            double& dist, int& tag) {
#ifdef DEBUG_CHECK_ASSERTION
    assert(point != NULL && direction != NULL && face != NULL);
#endif

    const double epson = 1.0e-12;
    dist = -1.0e-15;
    double x1, y1, z1, x2, y2, z2, AA, BB, CC, DD;
    double delta, tmp, t;

    // 平面上的两个切方向
    x1 = face[3] - face[0];
    y1 = face[4] - face[1];
    z1 = face[5] - face[2];
    x2 = face[6] - face[0];
    y2 = face[7] - face[1];
    z2 = face[8] - face[2];
    // 平面法向，平面方程 AA x + BB y + CC z + DD = 0
    AA = y1 * z2 - z1 * y2;
    BB = z1 * x2 - x1 * z2;
    CC = x1 * y2 - y1 * x2;
    DD = -1 * face[0] * AA - face[1] * BB - face[2] * CC;
    // 平面发向与射线切向的内积
    delta = AA * direction[0] + BB * direction[1] + CC * direction[2];
    // 射线方程(xm代表射线起点，xd代表射线方向) x = xm + xd t
    //                                    y = ym + yd t
    //                                    z = zm + zd t
    if (abs(delta) <= epson)
      tag = -1;  // 射线切向与平面平行
    else {
      // 射线与平面不平行，计算交点对应的t值
      tmp = -1 * AA * point[0] - BB * point[1] - CC * point[2] - DD;
      t = tmp / delta;
      double compare = 0;
      if (t >= compare)  //  t>0,解有意义，计算交点坐标
      {
        tag = 0;
        intersection[0] = point[0] + direction[0] * t;
        intersection[1] = point[1] + direction[1] * t;
        intersection[2] = point[2] + direction[2] * t;
        dist =
            t * sqrt(direction[0] * direction[0] + direction[1] * direction[1] +
                     direction[2] * direction[2]);
      } else
        //  t < 0，说明射线与平面的交点在射线的负方向，射线本上与平面无交点
        tag = -2;
    }
  }

  /**
   * @brief 将同一个面三维点降为二维点,去除法向量绝对值最大的分量维度
   *
   * @param points 输入三维点集合,行主序（存完第一个点，再存第二个点）
   * @param n 点个数
   * @param result 返回二维点集合
   */
  void reduceDim(const double* points, int n, double* result) {
    int dim_index = 0;  // 存储即将要去掉维度

    //  x,y,z存储法向量的各维度
    double x1 = points[0];
    double y1 = points[1];
    double z1 = points[2];
    double x2 = points[3];
    double y2 = points[4];
    double z2 = points[5];
    double x3 = points[6];
    double y3 = points[7];
    double z3 = points[8];
    double x = (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1);
    double y = (z2 - z1) * (x3 - x1) - (z3 - z1) * (x2 - x1);
    double z = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    double max = -1;
    if (abs(x) > max) {
      max = abs(x);
      dim_index = 0;
    }
    if (abs(y) > max) {
      max = abs(y);
      dim_index = 1;
    }
    if (abs(z) > max) {
      max = abs(z);
      dim_index = 2;
    }

    for (int i = 0; i < n; i++) {
      if (dim_index == 0) {
        result[i * 2] = points[i * 3 + 1];
        result[i * 2 + 1] = points[i * 3 + 2];
      } else if (dim_index == 1) {
        result[i * 2] = points[i * 3];
        result[i * 2 + 1] = points[i * 3 + 2];
      } else {
        result[i * 2] = points[i * 3];
        result[i * 2 + 1] = points[i * 3 + 1];
      }
    }
  }

  /**
   * @brief 计算射线与点组成的面相交，并判断交点是否在组成面的多边形中
   *
   * @param point 射线起点
   * @param direction 射线方向
   * @param face 组成面的点
   * @param n 组成面的点的个数
   * @param intersection 输出参数，交点坐标
   * @param dist 输出参数，射线起点到交点的距离
   * @param tag 输出参数，标识状态
   *              0：成功求出交点，并且交点在face面内部
   *             -1：射线与平面平行
   *             -2：射线所在直线与平面的交点在射线的负方向
   *             -3：射线与所在平面相交，但是交点在face外
   */
  void rayIntersectFace(const double* point, const double* direction,
                        const double* face, int n, double* intersection,
                        double& dist, int& tag) {
#ifdef DEBUG_CHECK_ASSERTION
    assert(n >= 3);
#endif
    // 取前三个点组成的面与射线相交，获取交点与距离
    double sub_face[9];
    for (int i = 0; i < 9; i++) sub_face[i] = face[i];
    lineSurfaceIntersect(point, direction, sub_face, intersection, dist, tag);

    if (tag == 0) {  // 有交点
      // 将交点与组成面的点降为2维
      double
          new_points[(n + 2) *
                     3];  // 包含组成face的点以及一个face[n]=face[0]点以及交点
      for (int i = 0; i < n * 3; i++) new_points[i] = face[i];
      new_points[n * 3] = face[0];
      new_points[n * 3 + 1] = face[1];
      new_points[n * 3 + 2] = face[2];
      new_points[(n + 1) * 3] = intersection[0];
      new_points[(n + 1) * 3 + 1] = intersection[1];
      new_points[(n + 1) * 3 + 2] = intersection[2];
      double result_points[(n + 2) * 2];
      reduceDim(new_points, n + 2, result_points);

      // 判断交点是否在face构成的多边形中
      double new_intersection[2] = {result_points[(n + 1) * 2],
                                    result_points[(n + 1) * 2 + 1]};
      if (wnPnpoly(new_intersection, result_points, n) != 0)
        return;
      else
        tag = -3;
    }
  }

  /**
   * @brief 判断点是否在网格中
   *
   * @param points 输入点坐标
   * @param n 点的数量
   * @return std::vector<int> 返回每个点所在网格单元，
   *                          -1代表在网格外，
   *                          >=0代表在网格中
   */
  std::vector<int> pointInGrid(const double* points, int n) {
#ifdef DEBUG_CHECK_ASSERTION
    assert(points != NULL);
    assert(n > 0);
#endif

    std::vector<int> result;
    for (int i = 0; i < n; i++) {
      int result_temp = -1;  // 记录当前点是否在该网格中
      double point[3];
      for (int d = 0; d < 3; d++) point[d] = points[i * 3 + d];
      // 第一步，通过包围盒和区间树大致判断点在那些单元之间，获取这些单元集合
      std::vector<int> bbox_intersect_cells;
      interval_tree_->getElementsListFromRange(point, point,
                                               bbox_intersect_cells);
      if (bbox_intersect_cells.size() == 0) {
        result.push_back(-1);
        continue;
      }

      //  第二步，遍历第一步得到的单元集合
      for (int c = 0; c < static_cast<int>(bbox_intersect_cells.size()); c++) {
        CellDescriptor cell = patch_cells_[bbox_intersect_cells[c]];

        // 采用cross number的方法来判断点是否在网格单元中有一个弊端，
        // 就是当点在单元边缘时,由于double误差关系，会使计算结果不准确，
        // 所以取点上下左右前后附近六个点进行测试，如果这7个点中有3个以
        // 上的点在该网格内，则认为该点在网格内
        bool point_result[7] = {false};
        const double around_points[21] = {point[0],
                                          point[1],
                                          point[2],
                                          point[0] * 1.0000000023879534,
                                          point[1] * 1.00000000568792345,
                                          point[2] * 0.9999999948626578,
                                          point[0] * 0.9999999978954623,
                                          point[1] * 0.9999999956321589,
                                          point[2] * 1.0000000098563214,
                                          point[0] * 0.9999999912356987,
                                          point[1] * 1.0000000045698753,
                                          point[2] * 1.0000000058746895,
                                          point[0] * 0.9999999932569874,
                                          point[1] * 1.0000000056987425,
                                          point[2] * 1.0000000074158963,
                                          point[0] * 0.9999999965874239,
                                          point[1] * 1.0000000065897412,
                                          point[2] * 0.9999999936985214,
                                          point[0] * 1.0000000074258635,
                                          point[1] * 1.0000000089632456,
                                          point[2] * 0.9999999923654789};
        for (int p = 0; p < 7; p++) {
          const double around_point[3] = {around_points[p * 3],
                                          around_points[p * 3 + 1],
                                          around_points[p * 3 + 2]};

          // 第三步，遍历该单元的所有face，从点出发取x方向的射线，利用rayIntersectFace计算是否与face相交
          // 记录交点个数，如果个数为偶数则该点在网格单元外，如果为奇数，则在网格单元内
          int cross_number = 0;
          double intersection[3];  // 存储交点
          double dist = 1.0e+100;  // 存储距离
          const double direction[3] = {1, 0, 0};
          for (int f = 0; f < static_cast<int>(cell.faces.size()); f++) {
            const std::vector<int>& face_nodes =
                patch_faces_[cell.faces[f]].nodes;
            double face_nodes_coord[face_nodes.size() * 3];
            for (int nnode = 0; nnode < static_cast<int>(face_nodes.size());
                 nnode++) {
              const hier::DoubleVector<3>& node_ref =
                  patch_nodes_[face_nodes[nnode]]->coord;
              face_nodes_coord[nnode * 3] = node_ref[0];
              face_nodes_coord[nnode * 3 + 1] = node_ref[1];
              face_nodes_coord[nnode * 3 + 2] = node_ref[2];
            }
            //  tag标识射线与面相交的状态，0代表有交点，-1代表平行，
            // -2代表与射线反方向相交，-3代表相交但交点不在face里面
            int tag = -4;
            double dist_temp = dist;
            rayIntersectFace(around_point, direction, face_nodes_coord,
                             face_nodes.size(), intersection, dist, tag);
            // 如果交点在face内，并且该交点没有在其他face出现过（对于棱的点或者顶点只统计一次）

            if (tag == 0 && abs(dist_temp - dist) > delta_) cross_number++;
          }  //  end loop cell face
          if (cross_number > 0 && (cross_number % 2 != 0))
            point_result[p] = true;
        }
        // 收集结果,目前的判断是当该点以及附近6个点中有3个点在该单元中，则认为该点在该单元中
        int result_true_sum = 0;
        for (int p = 0; p < 7; p++)
          if (point_result[p] == true) result_true_sum++;
        if (result_true_sum >= 3) {
          result_temp = cell.index;
          break;
        }
      }  //  end loop cell
      result.push_back(result_temp);
    }  //  end loop points
    return result;
  }

  /**
   * @brief 计算输入射线集合与网格的交点，输出交点坐标以及交点所在网格单元
   *
   * @param points 射线起点
   * @param directions 射线方向
   * @param n 射线数量
   * @param ids 每条射线与网格相交，交点所在网格单元索引，
   * 如果射线不相交,索引值为-1
   * @param intersections_coordinates 每条射线与网格相交，交点坐标，
   * 如果射线不相交，交点坐标为{-1.0e100, -1.0e100, -1.0e100}
   */
  void rayIntersectGrid(const double* points, const double* directions, int n,
                        std::vector<int>& ids,
                        double* intersections_coordinates) {
#ifdef DEBUG_CHECK_ASSERTION
    assert(points != NULL && directions != NULL &&
           intersections_coordinates != NULL);
#endif
    // 遍历每一条射线
    for (int r = 0; r < n; r++) {
      const double point[3] = {points[r * 3], points[r * 3 + 1],
                               points[r * 3 + 2]};
      const double direction[3] = {directions[r * 3], directions[r * 3 + 1],
                                   directions[r * 3 + 2]};
      double intersection[3] = {-1.0e100, -1.0e100, -1.0e100};
      double dist = 1.0e100;
      int intersect_cell_index = -1;
      // 遍历每一个外表面
      for (std::map<int, int>::iterator iter = patch_outer_faces_.begin();
           iter != patch_outer_faces_.end(); iter++) {
        double sub_intersection[3];
        double sub_dist;
        //  tag标识射线与面相交的状态，0代表有交点，-1代表平行，
        // -2代表与射线反方向相交，-3代表相交但交点不在face里面
        int tag;
        const std::vector<int>& face_nodes = patch_faces_[iter->first].nodes;
        double face_nodes_coord[face_nodes.size() * 3];
        for (int i = 0; i < static_cast<int>(face_nodes.size()); i++) {
#ifdef DEBUG_CHECK_ASSERTION
          assert(face_nodes[i] <= patch_faces_.size());
#endif
          const hier::DoubleVector<3>& node_ref =
              patch_nodes_[face_nodes[i]]->coord;
          face_nodes_coord[i * 3] = node_ref[0];
          face_nodes_coord[i * 3 + 1] = node_ref[1];
          face_nodes_coord[i * 3 + 2] = node_ref[2];
        }
        // 获取射线与face面的交点
        rayIntersectFace(point, direction, face_nodes_coord, face_nodes.size(),
                         sub_intersection, sub_dist, tag);
        if (tag == 0 && sub_dist < dist) {
          for (int i = 0; i < 3; i++) intersection[i] = sub_intersection[i];
          dist = sub_dist;
          intersect_cell_index = iter->second;
        }
      }  //  end loop patch outer face
      ids.push_back(intersect_cell_index);
      for (int i = 0; i < 3; i++)
        intersections_coordinates[r * 3 + i] = intersection[i];
    }  //  end loop ray
  }

  /**
   * @brief 根据输入目的网格单元索引，计算出与之相交源网格单元
   *
   * @param dest_patch 输入参数，目的网格
   * @param cells_index 输入参数，需要求交目的网格单元索引集合
   * @param intersect_index
   * 输出参数，与目的网格单元相交的源网格单元索引集合（已经进行去重处理）
   */
  void gridIntersectGrid(tbox::Pointer<hier::Patch<3> >& dest_patch,
                         std::vector<int>& cells_index,
                         std::set<int>& intersect_index) {
    int node_number = dest_patch->getNumberOfNodes(0);
    int cell_number = dest_patch->getNumberOfCells(0);
    tbox::Pointer<hier::PatchTopology<3> > topology =
        dest_patch->getPatchTopology();
    tbox::Pointer<hier::PatchGeometry<3> > patch_geometry =
        dest_patch->getPatchGeometry();

    // 获取节点坐标
    tbox::Pointer<pdat::NodeData<3, double> > nodes =
        patch_geometry->getNodeCoordinates();
    const double* nodes_coordinates = nodes->getPointer(0);

    //获取单元
    tbox::Array<int> cell_adj_nodes_extent;
    tbox::Array<int> cell_adj_nodes_indices;
    topology->getCellAdjacencyNodes(cell_adj_nodes_extent,
                                    cell_adj_nodes_indices);
    int focused_cells_size = cells_index.size();
    for (int cell = 0; cell < focused_cells_size; cell++) {
      int cell_index = cells_index[cell];
#ifdef DEBUG_CHECK_ASSERTION
      assert(cell_index >= 0 && cell_index < cell_number);
#endif

      //第一步，获取目的网格片每个网格单元的bbox
      double bbox[6] = {1.0e100,  1.0e100,  1.0e100,
                        -1.0e100, -1.0e100, -1.0e100};
      for (int j = cell_adj_nodes_extent[cell_index];
           j < cell_adj_nodes_extent[cell_index + 1]; j++) {
        int node_index = cell_adj_nodes_indices[j];
#ifdef DEBUG_CHECK_ASSERTION
        assert(node_index < node_number);
#endif
        for (int d = 0; d < 3; d++) {
          bbox[d] = std::min(bbox[d * 2], nodes_coordinates[node_index + d]);
          bbox[d + 3] =
              std::max(bbox[d * 2 + 1], nodes_coordinates[node_index + d]);
        }
      }
      //第二步，利用区间树求解与该网格单元bbox相交的源网格单元，并存放结果
      std::vector<int> bbox_intersect_cells;
      interval_tree_->getElementsListFromRange(bbox, bbox + 3,
                                               bbox_intersect_cells);
      for (int i = 0; i < static_cast<int>(bbox_intersect_cells.size()); i++)
        intersect_index.insert(bbox_intersect_cells[i]);
    }
  }

private:
  std::unordered_map<int,
                     NodeDescriptor*> patch_nodes_;  //  patch内节点集合
  std::unordered_map<int, FaceDescriptor> patch_faces_;  //  patch内面集合
  std::unordered_map<int, CellDescriptor> patch_cells_;  //  patch内体集合
  std::map<int, int> patch_outer_faces_;  //  patch外表面集合<face-id,cell-id>
  tbox::Pointer<hier::IntervalTree<3> > interval_tree_;  //  存储区间树
  double delta_;  //  用于double数值大小比较中的常数
};

//
//  指向GridIntersect的实现
//

GridIntersect::GridIntersect(tbox::Pointer<hier::Patch<3> >& source_patch)
    : impl_(new GridIntersectImpl(source_patch)) {}

GridIntersect::~GridIntersect() {}

std::vector<int> GridIntersect::pointInGrid(const double* points, int n) {
  return impl_->pointInGrid(points, n);
}

void GridIntersect::rayIntersectGrid(const double* start_points,
                                     const double* direction, int n,
                                     std::vector<int>& ids,
                                     double* intersection_coordinates) {
  impl_->rayIntersectGrid(start_points, direction, n, ids,
                          intersection_coordinates);
}

void GridIntersect::gridIntersectGrid(
    tbox::Pointer<hier::Patch<3> >& dest_patch, int& intersect_num,
    std::vector<int>& intersect_index) {
  std::vector<int> cells_index;
  std::set<int> intersect_indexs;
  int cell_number = dest_patch->getNumberOfCells(0);
  for (int i = 0; i < cell_number; i++) cells_index.push_back(i);

  impl_->gridIntersectGrid(dest_patch, cells_index, intersect_indexs);
  // 将set转化为vector
  for (std::set<int>::iterator iter = intersect_indexs.begin();
       iter != intersect_indexs.end(); iter++)
    intersect_index.push_back(*iter);
  intersect_num = intersect_indexs.size();
}

std::vector<int> GridIntersect::gridIntersectGrid(
    tbox::Pointer<hier::Patch<3> >& dest_patch, int focused_cell_index) {
  std::set<int> intersect_indexs;
  std::vector<int> cells_index(1, focused_cell_index);
  std::vector<int> intersect_index;

  impl_->gridIntersectGrid(dest_patch, cells_index, intersect_indexs);
  // 将set转化为vector
  for (std::set<int>::iterator iter = intersect_indexs.begin();
       iter != intersect_indexs.end(); iter++)
    intersect_index.push_back(*iter);
  return intersect_index;
}

#endif  //  GRID_INTERSECT_CPP