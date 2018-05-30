#ifndef GRID_INTERSECT_H
#define GRID_INTERSECT_H

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>

#include "Patch.h"
using namespace JAUMIN;

class GridIntersectImpl;

/**
 * @brief  GridIntersect
 * 该模块实现三个功能：1.查询输入点是否在输入网格内，如果在输出点所在网格单元
 * 2.给定网格外射线与网格，输出射线与网格外表面相加网格单元以及交点坐标
 * 3.给定源数据网格和目标单元网格，输出两个网格相交的源单元个数，以及源单元在网格片中的索引号
 *
 */

class GridIntersect : public boost::noncopyable {
public:
  /**
   * @brief 构造函数
   *
   * @param source_patch 用于点相交、射线相交、网格相交的目的网格
   */
  GridIntersect(tbox::Pointer<hier::Patch<3>>& source_patch);
  ~GridIntersect();

  /**
   * @brief 判断点是否在网格中
   *
   * @param points 输入点坐标
   * @param n 点的数量
   * @return std::vector<int> 返回每个点所在网格单元，
   *                          -1代表在网格外，
   *                          >=0代表在网格中
   */
  std::vector<int> pointInGrid(const double* points, int n);

  /**
   * @brief 点与网格外表面相交,输出相交网格单元编号与交点坐标
   *
   * @param start_points 射线起点集合
   * @param direction 射线方向集合
   * @param n 射线的数量
   * @param ids 输出参数，交点所在网格单元编号
   *        -1代表射线与网格外表面不相交，交点无意义
   *        >=0代表射线与网格外表面相交，交点所在网格单元索引
   * @param intersection 输出参数，交点坐标
   */
  void rayIntersectGrid(const double* start_points, const double* direction,
                        int n, std::vector<int>& ids,
                        double* intersection_coordinates);

  /**
   * @brief
   * 输入目的网格片，求解源网格片与目的网格片相交的单元个数以及源源单元索引号
   *
   * @param dest_patch 输入参数，目的网格片
   * @param intersect_num 输出参数，相交源网格单元个数
   * @param intersect_index 输出参数，相交源网格单元索引号
   */
  void gridIntersectGrid(tbox::Pointer<hier::Patch<3>>& dest_patch,
                         int& intersect_num, std::vector<int>& intersect_index);

  void gridIntersectGrid(
      tbox::Pointer<hier::Patch<3>>& dest_patch,
      std::vector<int>& focused_cells_index,
      std::vector<std::vector<int>>& focused_intersect_index);

private:
  boost::shared_ptr<GridIntersectImpl> impl_;
};

// #include "GridIntersect.cpp" 模板类加上
#endif  // GRID_INTERSECT_H
