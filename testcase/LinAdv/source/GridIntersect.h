//
// 文件名: GridIntersect.h
// 软件包: JAUMIN 
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 1 $
// 修改  : $Date: 2018-10-9 13:11:58 +0800
// 描述  : 网格相交算法类
//

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
 * 3.给定源数据网格和目的单元网格，输出两个网格相交的源单元个数，以及源单元在网格片中的索引号
 *
 */

class GridIntersect : public boost::noncopyable {
public:
  /**
   * @brief 构造函数
   *
   * @param source_patch 用于点相交、射线相交、网格相交的目的网格
   */
  GridIntersect(tbox::Pointer<hier::Patch<3> >& source_patch);
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
   * @brief 输入目的网格片，求解源网格片与目的网格片相交的单元个数
   * 以及源单元索引号。注：由于采用区间树算法，所以求出的结果不是
   * 所有单元都真正相交，真正相交的网格单元只是输出集合的子集
   *
   * @param dest_patch 输入参数，目的网格片
   * @param intersect_num 输出参数，相交源网格单元个数
   * @param intersect_index 输出参数，相交源网格单元索引号
   */
  void gridIntersectGrid(tbox::Pointer<hier::Patch<3> >& dest_patch,
                         int& intersect_num, std::vector<int>& intersect_index);

  /**
   * @brief 输入目的网格片，以及指定单元索引号，求解与该网格单元
   * 相交的源网格单元索引。注：由于采用区间树算法，所以求出的结果不是
   * 所有单元都真正相交，真正相交的网格单元只是输出集合的子集
   *
   * @param dest_patch 输入参数，目的网格片
   * @param focused_cell_index 输入参数，目的网格单元索引
   *
   * @return std::vector<int> 输出参数，与输入网格单元相交的源网格单元索引集合
   */
  std::vector<int> gridIntersectGrid(tbox::Pointer<hier::Patch<3> >& dest_patch,
                                     int focused_cell_index);

private:
  boost::shared_ptr<GridIntersectImpl> impl_;
};

// #include "GridIntersect.cpp" 模板类加上
#endif  // GRID_INTERSECT_H
