//
// 文件名: LinAdv.C
// 软件包: JAUMIN applications
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 146 $
// 修改  : $Date: 2011-08-26 17:09:34 +0800 (五, 2011-08-26) $
// 描述  : 线性对流问题的网格片时间积分算法类的实现.
//

#include "PIO.h"
#include "Patch.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "EntityUtilities.h"
#include "IntegratorComponent.h"
#include "EdgeData.h"
#include "CellData.h"
#include "FaceData.h"
#include "LinAdv.h"
#include "LinAdvFort.h"
#include "Database.h"

/*************************************************************************
 * 构造函数.
 *************************************************************************/
LinAdv::LinAdv(const string& object_name,
               tbox::Pointer<tbox::Database> input_db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!object_name.empty());
  TBOX_ASSERT(!input_db.isNull());
#endif

  d_object_name = object_name;

  // 从输入文件或重启动文件读入数据.
  bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
  if (is_from_restart) {
    getFromRestart();
  }

  getFromInput(input_db);

  // 注册变量和数据片.
  registerModelVariables();

  // 将当前类对象注册为重启动对象.
  tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
}

/*************************************************************************
 * 构造函数.
 ************************************************************************/
LinAdv::~LinAdv() {
  tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
};

/*************************************************************************
 * 注册变量和数据片.
 ************************************************************************/
void LinAdv::registerModelVariables() {
  // 获取变量数据库.
  hier::VariableDatabase<NDIM>* variable_db =
      hier::VariableDatabase<NDIM>::getDatabase();

  // 定义变量.
  // 守恒量。
  tbox::Pointer<pdat::CellVariable<NDIM, double> > d_uval =
      new pdat::CellVariable<NDIM, double>("uval", 1);  // 中心量，深度为1.

#if (NDIM == 2)
  // 二维情况下的通量。
  tbox::Pointer<pdat::EdgeVariable<NDIM, double> > d_flux =
      new pdat::EdgeVariable<NDIM, double>("flux", 1);  // 边心量，深度为1.
#else
  // 三维情况下的通量。
  tbox::Pointer<pdat::FaceVariable<NDIM, double> > d_flux =
      new pdat::FaceVariable<NDIM, double>("flux", 1);  // 边心量，深度为1.
#endif

  // 辅助变量，用于在后处理时标示网格单元索引.
  tbox::Pointer<pdat::CellVariable<NDIM, double> > d_label =
      new pdat::CellVariable<NDIM, double>("label", 1);  // 中心量，深度为1.
  tbox::Pointer<pdat::NodeVariable<NDIM, double> > d_node =
      new pdat::NodeVariable<NDIM, double>("node", 1);  // 中心量，深度为1.

  // 当前值上下文, 新值上下文, 演算上下文.
  tbox::Pointer<hier::VariableContext> d_current =
      variable_db->getContext("CURRENT");
  tbox::Pointer<hier::VariableContext> d_new = variable_db->getContext("NEW");
  tbox::Pointer<hier::VariableContext> d_scratch =
      variable_db->getContext("SCRATCH");

  d_gw = 1;  // 影像区宽度.

  // 数据片<uval,current>: 存储当前时刻的 u 值, 影像区宽度为0.
  d_node_id = variable_db->registerVariableAndContext(d_node, d_current);
  d_uval_current_id =
      variable_db->registerVariableAndContext(d_uval, d_current);

  // 数据片<uval,new>: 存储新时刻的 u 值, 影像区宽度为0.
  d_uval_new_id = variable_db->registerVariableAndContext(d_uval, d_new);

  // 数据片<flux,new>: 存储新时刻的 f 值, 影像区宽度为0.
  d_flux_new_id = variable_db->registerVariableAndContext(d_flux, d_new);

  // 数据片<uval,scratch>: 存储变量 u 的演算值, 影像区宽度为1.
  d_uval_scratch_id =
      variable_db->registerVariableAndContext(d_uval, d_scratch, d_gw);

  // 数据片<label,current>: 存储当前时刻的 label 值, 影像区宽度为0.
  d_label_current_id =
      variable_db->registerVariableAndContext(d_label, d_current);
}

/*************************************************************************
 * 注册绘图量.
 *************************************************************************/
void LinAdv::registerPlotData(
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!(javis_writer.isNull()));
#endif

  // 注册守恒量, 即: u 在当前时刻的值.
  javis_writer->registerPlotQuantity("UUU", "SCALAR", d_uval_current_id);
  // 注册辅助量.
  // javis_writer->registerPlotQuantity("LABEL", "SCALAR", d_label_current_id);
}

/*************************************************************************
 *  初始化指定的积分构件.
 *
 *  注册待填充的数据片或待调度内存空间的数据片到积分构件.
 ************************************************************************/
void LinAdv::initializeComponent(
    algs::IntegratorComponent<NDIM>* component) const {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(component);
#endif

  const string& component_name = component->getName();

  if (component_name == "ALLOC_NEW_DATA") {  // 内存构件
    component->registerPatchData(d_uval_new_id);
    component->registerPatchData(d_flux_new_id);

  } else if (component_name == "ALLOC_SCRATCH_DATA") {  // 内存构件
    component->registerPatchData(d_uval_scratch_id);

  } else if (component_name == "INIT_SET_VALUE") {  // 初值构件
    component->registerInitPatchData(d_node_id);
    component->registerInitPatchData(d_uval_current_id);
    component->registerInitPatchData(d_label_current_id);

  } else if (component_name == "STEP_SIZE") {  // 步长构件

  } else if (component_name == "COMPUTE_FLUX") {  // 数值构件, 计算通量
    component->registerCommunicationPatchData(d_uval_scratch_id,
                                              d_uval_current_id);

  } else if (component_name == "SET_BDRY") {  // 数值构件, 设置物理边界条件

  } else if (component_name == "CONSER_DIFF") {  // 数值构件, 计算守恒量

  } else if (component_name == "COPY_SOLUTION") {  // 复制构件
    component->registerCopyPatchData(d_uval_current_id, d_uval_new_id);
  } else {
    TBOX_ERROR("\n::initializeComponent() : component "
               << component_name << " is not matched. " << endl);
  }
}

/*************************************************************************
 *  初始化数据片（支持初值构件）.
 ************************************************************************/
void LinAdv::initializePatchData(hier::Patch<NDIM>& patch, const double time,
                                 const bool initial_time,
                                 const string& component_name) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(component_name == "INIT_SET_VALUE");
#endif

  NULL_USE(time); /**< 初始化中没有用到time */

  if (initial_time) {
    // 获取当前时刻的守恒量.
    tbox::Pointer<pdat::CellData<NDIM, double> > uval_current =
        patch.getPatchData(d_uval_current_id);

    // 获取网格标记量。
    tbox::Pointer<pdat::CellData<NDIM, double> > label_current =
        patch.getPatchData(d_label_current_id);

    // 获取当前网格片关联的网格片几何对象和网格片拓扑对象.
    tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
        patch.getPatchGeometry();
    tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
        patch.getPatchTopology();

    // 获取网格单元的中心坐标.
    tbox::Pointer<pdat::CellData<NDIM, double> > cell_coords =
        patch_geo->getCellCoordinates();
    // 获取网格结点坐标.
    tbox::Pointer<pdat::NodeData<NDIM, double> > node_coords =
        patch_geo->getNodeCoordinates();

    // node_coords->print(0);

    tbox::Pointer<pdat::NodeData<NDIM, double> > node_current =
        patch.getPatchData(d_node_id);
    // 获取当前网格片的内部单元数目.
    int num_cells = patch.getNumberOfCells();
    // 获取网格片内部网格结点数目.
    int num_nodes = patch.getNumberOfNodes();
    node_current->fill(0.0);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_cells <= cell_coords->getNumberOfEntity());
    TBOX_ASSERT(num_nodes <= node_coords->getNumberOfEntity());
#endif

    // 获取<单元, 结点>的拓扑关系.
    tbox::Array<int> cell_adj_node_extent, cell_adj_node_index;
    patch_top->getCellAdjacencyNodes(cell_adj_node_extent, cell_adj_node_index);

    // 调用Fortran程序初始化.
    initsetvalue_(NDIM, d_dx, num_cells, num_nodes, cell_coords->getPointer(),
                  node_coords->getPointer(), cell_adj_node_extent.getPointer(),
                  cell_adj_node_index.getPointer(), uval_current->getPointer(),
                  label_current->getPointer());

    label_current->fillAll(patch.getIndex());
  }
}

/*************************************************************************
 *  计算稳定性时间步长（支持步长构件）.
 ************************************************************************/
double LinAdv::getPatchDt(hier::Patch<NDIM>& patch, const double time,
                          const bool initial_time, const int flag_last_dt,
                          const double last_dt, const string& component_name) {
  // 返回时间步长.
  return (0.8 * d_dx / d_x_velocity);
}

/*************************************************************************
 * 完成单个网格片上的数值计算（支持数值构件）.
 * 该函数基于显式迎风格式，实现通量和守恒量的计算.
 ************************************************************************/
void LinAdv::computeOnPatch(hier::Patch<NDIM>& patch, const double time,
                            const double dt, const bool initial_time,
                            const string& component_name) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(component_name == "COMPUTE_FLUX" ||
              component_name == "CONSER_DIFF" || component_name == "SET_BDRY");
#endif

  // tbox::Pointer< pdat::NodeData<NDIM,double> > node_extend_coords =
  // patch.getPatchGeometry()->getNodeCoordinates();
  // node_extend_coords->print();

  if (component_name == "COMPUTE_FLUX") {
    computeFluxOnPatch(patch, time, dt, initial_time, component_name);
  } else if (component_name == "CONSER_DIFF") {
    conserDiffOnPatch(patch, time, dt, initial_time, component_name);
  } else if (component_name == "SET_BDRY") {
    setPhysicalBoundaryOnPatch(patch, time, dt, initial_time, component_name);
  }
}

/*************************************************************************
 * 该函数基于显式迎风格式，计算通量.
 ************************************************************************/
void LinAdv::computeFluxOnPatch(hier::Patch<NDIM>& patch, const double time,
                                const double dt, const bool initial_time,
                                const string& component_name) {
#ifdef DEBUG_CHECK_ASSERTIONS
#endif

  // 获取当前patch关联的网格片几何对象和网格片拓扑对象.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();
  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
      patch.getPatchTopology();

  // 获取网格片内部及前d_gw层影像区中的网格单元数目的和.
  int num_extend_cells = patch.getNumberOfCells(d_gw);

  int num_extend_nodes = patch.getNumberOfNodes(d_gw);
  // 获取网格结点坐标.
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_extend_coords =
      patch_geo->getNodeCoordinates();
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(node_extend_coords->getNumberOfEntity() == num_extend_nodes);
#endif

  // 获取网格单元中心坐标数组.
  tbox::Pointer<pdat::CellData<NDIM, double> > cell_extend_coords =
      patch_geo->getCellCoordinates();

#ifdef DEBUG_CHECK_ASSERTIONS

  TBOX_ASSERT(cell_extend_coords->getNumberOfEntity() == num_extend_cells);
#endif

  // 获取中心量演算数据片.
  tbox::Pointer<pdat::CellData<NDIM, double> > uval_scratch =
      patch.getPatchData(d_uval_scratch_id);

#if (NDIM == 2)
  // 获取网格片内部的边的数目.
  int num_edges = patch.getNumberOfEdges();

  // 获取边的中点坐标.
  tbox::Pointer<pdat::EdgeData<NDIM, double> > edge_coords =
      patch_geo->getEdgeCoordinates();

  // 获取<边, 单元>的拓扑关系.
  tbox::Array<int> edge_adj_cell_extent, edge_adj_cell_index;
  patch_top->getEdgeAdjacencyCells(edge_adj_cell_extent, edge_adj_cell_index);

  // 获取通量数据片(边心量).
  tbox::Pointer<pdat::EdgeData<NDIM, double> > flux_new =
      patch.getPatchData(d_flux_new_id);

  // 调用Fortran程序计算通量.
  advancepatch2d_(d_x_velocity, dt, num_edges, num_extend_cells,
                  edge_coords->getPointer(), cell_extend_coords->getPointer(),
                  edge_adj_cell_extent.getPointer(),
                  edge_adj_cell_index.getPointer(), uval_scratch->getPointer(),
                  flux_new->getPointer());

#else

  // 获取网格片内部的面的数目.
  int num_faces = patch.getNumberOfFaces();

  // 获取面的中心坐标.
  tbox::Pointer<pdat::FaceData<NDIM, double> > face_coords =
      patch_geo->getFaceCoordinates();

  // 获取<面, 单元>的拓扑关系.
  tbox::Array<int> face_adj_cell_extent, face_adj_cell_index;
  patch_top->getFaceAdjacencyCells(face_adj_cell_extent, face_adj_cell_index);

  // 获取通量数据片(面心量).
  tbox::Pointer<pdat::FaceData<NDIM, double> > flux_new =
      patch.getPatchData(d_flux_new_id);

  // 调用Fortran程序计算通量.
  advancepatch3d_(d_x_velocity, dt, num_faces, num_extend_cells,
                  face_coords->getPointer(), cell_extend_coords->getPointer(),
                  face_adj_cell_extent.getPointer(),
                  face_adj_cell_index.getPointer(), uval_scratch->getPointer(),
                  flux_new->getPointer());
#endif
}

/*************************************************************************
 * 该函数基于显式迎风格式，更新守恒量.
 ************************************************************************/
void LinAdv::conserDiffOnPatch(hier::Patch<NDIM>& patch, const double time,
                               const double dt, const bool initial_time,
                               const string& component_name) {
#ifdef DEBUG_CHECK_ASSERTIONS
#endif

  // 获取当前patch关联的网格片几何对象和网格片拓扑对象.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();
  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
      patch.getPatchTopology();

  // 获取网格片内部网格单元的数目.
  int num_cells = patch.getNumberOfCells();

  // 获取网格片内部单元的中心坐标.
  tbox::Pointer<pdat::CellData<NDIM, double> > cell_coords =
      patch_geo->getCellCoordinates();

  // 获取当前时刻守恒量数据片.
  tbox::Pointer<pdat::CellData<NDIM, double> > uval_current =
      patch.getPatchData(d_uval_current_id);

  // 获取新时刻的守恒量数据片.
  tbox::Pointer<pdat::CellData<NDIM, double> > uval_new =
      patch.getPatchData(d_uval_new_id);

#if (NDIM == 2)

  // 获取网格片内部边的数目.
  int num_edges = patch.getNumberOfEdges();

  // 获取边的中点坐标.
  tbox::Pointer<pdat::EdgeData<NDIM, double> > edge_coords =
      patch_geo->getEdgeCoordinates();

  // 获取<单元, 边>的拓扑关系.
  tbox::Array<int> cell_adj_edge_extent, cell_adj_edge_index;
  patch_top->getCellAdjacencyEdges(cell_adj_edge_extent, cell_adj_edge_index);

  // 获取下一时间步的通量数据片.
  tbox::Pointer<pdat::EdgeData<NDIM, double> > flux_new =
      patch.getPatchData(d_flux_new_id);

  // 调用Fortran程序, 更新守恒量.
  conspatch2d_(d_dx, num_cells, num_edges, cell_coords->getPointer(),
               edge_coords->getPointer(), cell_adj_edge_extent.getPointer(),
               cell_adj_edge_index.getPointer(), uval_current->getPointer(),
               uval_new->getPointer(), flux_new->getPointer());
#else

  // 获取网格片内部的面的数目.
  int num_faces = patch.getNumberOfFaces();

  // 获取面的中心坐标.
  tbox::Pointer<pdat::FaceData<NDIM, double> > face_coords =
      patch_geo->getFaceCoordinates();

  // 获取<单元, 面>的拓扑关系.
  tbox::Array<int> cell_adj_face_extent, cell_adj_face_index;
  patch_top->getCellAdjacencyFaces(cell_adj_face_extent, cell_adj_face_index);

  // 获取下一时间步的通量数据片.
  tbox::Pointer<pdat::FaceData<NDIM, double> > flux_new =
      patch.getPatchData(d_flux_new_id);

  // 调用Fortran程序, 更新守恒量.
  conspatch3d_(d_dx, num_cells, num_faces, cell_coords->getPointer(),
               face_coords->getPointer(), cell_adj_face_extent.getPointer(),
               cell_adj_face_index.getPointer(), uval_current->getPointer(),
               uval_new->getPointer(), flux_new->getPointer());
#endif
}

/*************************************************************************
 *  填充物理边界条件.
 ************************************************************************/
void LinAdv::setPhysicalBoundaryOnPatch(hier::Patch<NDIM>& patch,
                                        const double time, const double dt,
                                        const bool initial_time,
                                        const string& component_name) {
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(component_name == "SET_BDRY");
#endif

  // 获取当前patch关联的网格片几何对象和网格片拓扑对象.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();
  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
      patch.getPatchTopology();

  // 前处理中定义的网格实体集合编号和类型.
  int set_id;  // 暂时写死，将来从输入文件读入。
  hier::EntityUtilities::EntityType set_type;

#if (NDIM == 2)
  set_id = 1;  // 暂时写死，将来从输入文件读入。
  set_type = hier::EntityUtilities::EDGE;

  // 获取通量数据片, 边心量.
  tbox::Pointer<pdat::EdgeData<NDIM, double> > flux =
      patch.getPatchData(d_flux_new_id);

  // 获取当前网格片内部边的数目。
  int num_inter = patch.getNumberOfEdges();

#else  // 三维情况

  set_id = 1;
  set_type = hier::EntityUtilities::FACE;

  // 获取通量数据片, 面心量.
  tbox::Pointer<pdat::FaceData<NDIM, double> > flux =
      patch.getPatchData(d_flux_new_id);

  // 获取当前网格片内部面的数目。
  int num_inter = patch.getNumberOfFaces();

#endif

  if (patch.hasEntitySet(set_id, set_type)) {
    // 获取指定编号和类型的集合包含的网格实体的索引。
    const tbox::Array<int>& entity_idx =
        patch_geo->getEntityIndicesInSet(set_id, set_type);

    // 获取物理边界上的边或者面的数目.
    int size = entity_idx.getSize();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(size > 0);
#endif

    for (int i = 0; i < size; ++i) {
      if (entity_idx[i] <
          num_inter) {  // 只为定义在网格片内部的边或者面设置通量。
        flux->getPointer()[entity_idx[i]] = dt * 2.0 * d_x_velocity;
      }
    }
  }
}

/*************************************************************************
 *  从输入数据库读入数据.
 ************************************************************************/
void LinAdv::getFromInput(tbox::Pointer<tbox::Database> db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!db.isNull());
#endif

  if (db->keyExists("constant_x_velocity")) {
    d_x_velocity = db->getDouble("constant_x_velocity");
    TBOX_ASSERT(d_x_velocity > 0.0);
  } else {
    TBOX_ERROR(d_object_name << ": "
                             << " No key `constant_x_velocity' found in data."
                             << endl);
  }
}

/*************************************************************************
 *  输出数据成员到重启动数据库.
 ************************************************************************/
void LinAdv::putToDatabase(tbox::Pointer<tbox::Database> db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!db.isNull());
#endif

  db->putDouble("d_dx", d_dx);
}

/*************************************************************************
 *  从重启动数据库读取数据.
 ************************************************************************/
void LinAdv::getFromRestart() {
  tbox::Pointer<tbox::Database> root_db =
      tbox::RestartManager::getManager()->getRootDatabase();

  tbox::Pointer<tbox::Database> db;
  if (root_db->isDatabase(d_object_name)) {
    db = root_db->getDatabase(d_object_name);
  } else {
    TBOX_ERROR("Restart database corresponding to "
               << d_object_name << " not found in restart file.");
  }

  d_dx = db->getDouble("d_dx");
}
