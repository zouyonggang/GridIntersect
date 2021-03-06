# GridIntersect
GridIntersect模块主要是为了解决点与网格相交、射线与网格相交、网格与网格相交三类问题

## 一、功能特性
1. 给定点以及三维非结构网格，计算点是否在网格中，如果点在网格中，输出网格单元编号
2. 给定三维非结构网格以及网格外的一条射线，计算射线与网格外表面的交点坐标以及交点所在单元编号
3. 给定两个三维非结构网格，计算两个网格中哪些单元相交，并求出交点所在单元编号

## 二、接口说明
```
  /**
   * @brief 构造函数
   *
   * @param source_patch 用于点相交、射线相交、网格相交的源网格片（jaumin）
   */
  GridIntersect(tbox::Pointer<hier::Patch<3>> source_patch);
```
```
  /**
   * @brief 判断点是否在网格中
   *
   * @param points 输入点坐标
   * @param n 点的数量
   * @return std::vector<int> 返回每个点所在网格单元，
   *                          -1代表在网格外，
   *                          >=0代表在网格中
   */
  std::vector<int> pointInGrid(const double *points, int n);
```
```
  /**
   * @brief 与网格外表面相交,输出相交网格单元编号与交点坐标
   *
   * @param start_points 射线起点集合
   * @param direction 射线方向集合
   * @param n 射线的数量
   * @param ids 输出参数，交点所在网格单元编号
   *        -1代表射线与网格外表面不相交，交点无意义
   *        >=0代表射线与网格外表面相交，交点所在网格单元索引
   * @param intersection 输出参数，交点坐标
   */
  void rayIntersectGrid(const double *start_points, const double *direction,
                        int n, std::vector<int> &ids,
                        double *intersection_coordinates);
```
```
  /**
   * @brief
   * 输入目的网格片，求解源网格片与目的网格片相交的单元个数以及源源单元索引号
   *
   * @param dest_patch 输入参数，目的网格片
   * @param intersect_num 输出参数，相交源网格单元个数
   * @param intersect_index 输出参数，相交源网格单元索引号
   */
  void gridIntersectGrid(tbox::Pointer<hier::Patch<3> >& dest_patch,
                         int& intersect_num, std::vector<int>& intersect_index);
```
```
  /**
   * @brief
   * 输入目的网格片，以及指定单元索引号，求解与该网格单元相交的源网格单元索引
   *
   * @param dest_patch 输入参数，目的网格片
   * @param focused_cell_index 输入参数，网格单元所以
   *
   * @return std::vector<int> 输出参数，与输入网格单元相交的源网格单元索引集合
   */
  std::vector<int> gridIntersectGrid(tbox::Pointer<hier::Patch<3>>& dest_patch,
                                     int focused_cell_index);
```

## 三、模块依赖
该模块因使用jaumin底层数据结构，以及区间树算法，所以依赖jaumin，无其他三方库依赖项.

## 四、编译和使用
1. 本模块只包含一个类,用户可以将两个源码文件放入基于jaumin框架研发的应用程序中直接使用
2. 本模块也提供cmake编译，在linux环境下执行以下命令生成动态库(jaumin版本大于1.9.0)
```bash
cd build
cmake -DJAUMIN_ROOT=path/to/jaumin ../source  (jaumin>1.9.0)
cmake -DJAUMIN=path/to/jaumin ../source       (jaumin<1.9.0)
make 
```

## 五、测试
本模块基本jaumin示例程序LinAdv设计了一份测试用例，网格来源使用“二维网格转三维网格”任务生成的网格。
###功能测试
1. 针对点与网格相交使用以下几类输入
- 内部点测试: 取网格所有内部单元体坐标作为内部点输入，输出应该为[0-网格单元总数]的递增数列
  * 结果：测试通过，但是某些特殊网格单元（网格单元过小）的中心点会被判到相邻网格单元
- 边界点测试：取网格顶点中x值最大的点作为输入，输出应该是该点所在的某个网格单元
  * 结果：测试通过
- 外部点测试：取影像区网格单元中心点作为输入（如果网格片没有影像区，结果无意义）,输出应该为-1
  * 结果：测试通过，但是某些特殊影像区网格单元中心点会被判到网格中
2. 针对射线与网格相交
- 取x值最大的节点作为原始点，并将该点x值加一，将新点作为起点，（-1，0，0）作为方向，输出射线与网格外表面交点坐标
  * 结果：测试通过
- 取x值最大的节点作为原始点，并将该点x值加一，将新点作为起点，（1，0，0）作为方向，输出射线与网格交点单元为-1，并且交点坐标无意义
  * 结果：测试通过
3. 针对网格与网格相交
- 取第一个patch作为源网格，让网格层中其他patch与之相交
  * 结果：测试通过

###性能测试
网格共94万个网格单元，2万个外表面，测试机配置为：8核Intel i7-3770 + 8G内存

| 数量       |  10e+4  |  10e+5  |
| --------   | :-----: | :----:  |
| 点         |  2.24s  |  33.02s |
| 射线       |  45.51s | 422.39s |
| 网格       |  19.01s | 133.19s |



## 六、算法设计
1. 点与网格相交
- a.获取网格片中的每个网格单元
- b.为网格单元构造包围盒
- c.通过区间树算法大概得到输入点属于那些网格单元，遍历这些网格单元
- d.获取组成网格单元组成的面，遍历每一个面
- e.以输入点作为起点，任意方向做一条射线，判断射线是否与面相交
- f.如果相交，判断交点是否在面中，如果在面中，进行计数
- g.遍历完所有面，查看交点总数，如果为偶数则该点在网格单元外，如果为奇数该点在网格单元中
- h.取该点上下左右前后六个点是否也在该网格单元中，如果有3个以上在，则认为该点在该网格单元中

2. 射线与网格外表面相交
- a.遍历该网格所有面，如果该面只被一个网格单元包围，则该面属于外表面，存储所有外表面
- b.遍历所有外表面，计算射线与该面的交点
- c.计算交点是否在该面内，如果是则记录射线起点到交点的距离
- d.重复步骤c，得到最小距离交点坐标，该坐标为射线与网格的交点坐标

3. 网格与网格相交
- a.遍历目的网格所有网格单元，并计算出网格的bbox
- b.在源网格每个网格单元bounding box构成的区间树中，求解出与第一步网格单元bounding box相交的网格单元
- c.输出与所有目的网格单元相交的源单元个数以及索引
