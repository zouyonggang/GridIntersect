//
// 文件名: LinAdv.h
// 软件包: JAUMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 138 $
// 修改  : $Date: 2007-09-28 13:11:58 +0800
// 描述  : 线性对流问题的网格片时间积分算法类.
//

#ifndef included_LinAdv
#define included_LinAdv

#include "Pointer.h"
#include "Array.h"
#include "RestartManager.h"
#include "JaVisDataWriter.h"
#include "DoubleVector.h"
#include "CellData.h"
#include "NodeData.h"
#include "StandardComponentPatchStrategy.h"
#include "NodeVariable.h"
#include "CellVariable.h"
#include "EdgeVariable.h"
#include "FaceVariable.h"
#include "MPI.h"
#include "Utilities.h"

using namespace std;
using namespace JAUMIN;

/**
 * @brief 该类从标准构件网格片策略类 algs::StandardComponentPatchStrategy 派生,
 * 实现求解线性对流方程的数值计算子程序.
 *
 * 该类需要从输入文件读取如下参数.
 *     - \b constant_x_velocity \n
 *       双精度浮点型, X方向的对流速度.
 */
class LinAdv : public algs::StandardComponentPatchStrategy<NDIM>,
               public tbox::Serializable {
public:
  /*!
   * @brief 构造函数.
   * @param object_name 输入参数, 字符串, 表示对象名称.
   * @param input_db    输入参数, 指针,   指向输入数据库.
   *
   * @note
   * 该函数主要完成以下操作:
   *  -# 初始化内部数据成员;
   *  -# 定义变量和数据片.
   */
  LinAdv(const string& object_name, tbox::Pointer<tbox::Database> input_db);

  /*!
   * @brief 析构函数.
   */
  virtual ~LinAdv();

  /// @name 重载基类 algs::StandardComponentPatchStrategy<NDIM> 的函数:
  // @{

  /*!
   * @brief 初始化指定的积分构件.
   *
   * 注册待填充的数据片或待调度内存空间的数据片到积分构件.
   *
   * @param component 输入参数, 指针, 指向待初始化的积分构件对象.
   */
  void initializeComponent(algs::IntegratorComponent<NDIM>* component) const;

  /**
   * @brief 初始化数据片（支持初值构件）.
   *
   * @param patch          输入参数, 网格片类,     表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示初始化的时刻.
   * @param initial_time   输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
   * @param component_name 输入参数, 字符串, 当前调用该函数的初值构件的名称.
   */
  void initializePatchData(hier::Patch<NDIM>& patch, const double time,
                           const bool initial_time,
                           const string& component_name);

  /*!
   * @brief 完成单个网格片上的数值计算（支持数值构件）.
   *
   * 该函数基于显式迎风格式，实现通量和守恒量的计算。
   *
   * @param patch          输入参数, 网格片类,     表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param initial_time   输入参数, 逻辑型,       表示当前是否为初始时刻.
   * @param component_name 输入参数, 字符串, 当前调用该函数的数值构件的名称.
   */
  void computeOnPatch(hier::Patch<NDIM>& patch, const double time,
                      const double dt, const bool initial_time,
                      const string& component_name);

  /*!
   * @brief 在单个网格片上计算稳定性时间步长（支持时间步长构件）.
   *
   * @param patch          输入参数, 网格片类,     表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示初始化的时刻.
   * @param initial_time   输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
   * @param flag_last_dt   输入参数, 整型,         表示前次积分返回的状态.
   * @param last_dt        输入参数, 双精度浮点型, 表示前次积分的时间步长.
   * @param component_name 输入参数, 字符串, 当前调用该函数的时间步长构件的名称.
   *
   * @return 双精度浮点型, 时间步长.
   */
  double getPatchDt(hier::Patch<NDIM>& patch, const double time,
                    const bool initial_time, const int flag_last_dt,
                    const double last_dt, const string& component_name);

  /**
   * @brief 根据边界条件, 填充物理边界影像区数据.
   *
   * @param patch          输入参数, 网格片类,     表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示求解时间步长.
   * @param initial_time   输入参数, 逻辑型,       表示当前是否为初始时刻.
   * @param component_name 输入参数, 字符串, 当前调用该函数的积分构件的名称.
   */
  void setPhysicalBoundaryOnPatch(hier::Patch<NDIM>& patch, const double time,
                                  const double dt, const bool initial_time,
                                  const string& component_name);

  //@}

  ///@name 重载基类tbox::Serializable的函数
  //@{
  /*!
   * @brief 将数据成员输出到重启动数据库.
   * @param db 输入参数, 指针, 指向重启动数据库.
   */
  void putToDatabase(tbox::Pointer<tbox::Database> db);

  //@}

  ///@name 自定义函数
  //@{

  /*!
   * @brief 注册绘图量.
   * @param javis_writer 输入参数, 指针, 表示 JaVis 数据输出器.
   */
  void registerPlotData(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer);

  //@}

private:
  /*!@brief 计算通量.  */
  void computeFluxOnPatch(hier::Patch<NDIM>& patch, const double time,
                          const double dt, const bool initial_time,
                          const string& component_name);

  /*!@brief 守恒差分.  */
  void conserDiffOnPatch(hier::Patch<NDIM>& patch, const double time,
                         const double dt, const bool initial_time,
                         const string& component_name);

  /*!@brief 从输入数据库读入数据.  */
  void getFromInput(tbox::Pointer<tbox::Database> db);

  /*!@brief 从重启动数据库读入数据.  */
  void getFromRestart();

  /*!@brief 注册变量和数据片.  */
  void registerModelVariables();

  string d_object_name; /*!@brief 对象名.  */
  int d_gw;             /*! 影像区宽度. */
  double d_dx;          /*! X方向的网格步长. */
  double d_x_velocity;  /*! @brief X方向速度常量 */

  int d_uval_current_id; /*!@brief <d_uval, d_current> 数据片的索引号 */
  int d_uval_new_id;     /*!@brief <d_uval, d_new> 数据片的索引号 */
  int d_uval_scratch_id; /*!@brief <d_uval, scratch> 数据片的索引号 */
  int d_flux_new_id;     /*!@brief <d_flux, new> 数据片的索引号 */

  int d_label_current_id; /*!@brief <d_label, d_current> 数据片的索引号 */

  int d_node_id;
};

#endif
