//
// 文件名: LinAdvFort.h
// 软件包: JAUMIN application
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 95 $
// 修改  : $Date: 2011-05-30 09:57:32 +0800 (一, 2011-05-30) $
// 描述  : F77外部子程序接口.
//

extern "C" {

void initsetvalue_(const int &,     // DIM
                   double &,        // d_dx
                   const int &,     // num_cells
                   const int &,     // num_extend_nodes
                   const double *,  // cell_coords
                   const double *,  // node_extend_coords
                   const int *,     // cell_adj_node_extent
                   const int *,     // cell_adj_node_index
                   double *,        // cell_val_current
                   double *);       // cell_label_current

#if (NDIM == 2)
void advancepatch2d_(const double &,  // d_x_velocity
                     const double &,  // dt
                     const int &,     // num_edges
                     const int &,     // num_extend_cells
                     const double *,  // edge_coords
                     const double *,  // cell_extend_coords
                     const int *,     // edge_adj_cell_extent
                     const int *,     // edge_adj_cell_index
                     const double *,  // uval_sctatch
                     double *);       // flux_new

void conspatch2d_(const double &,   // d_dx
                  const int &,      // num_cells
                  const int &,      // num_edges
                  const double *,   // cell_coords
                  const double *,   // edge_coords
                  const int *,      // cell_adj_edge_extent
                  const int *,      // cell_adj_edge_index
                  double *,         // uval_current
                  double *,         // uval_new
                  const double *);  // flux_new

#else
//  void initsetvalue3d_(double&,          // d_dx
//                       const int&,       // num_cells
//                       const int&,       // num_extend_nodes
//                       const double*,    // cell_coords
//                       const double*,    // node_extend_coords
//                       const int*,       // cell_adj_node_extent
//                       const int*,       // cell_adj_node_index
//                       double*);         // cell_val_current

void advancepatch3d_(const double &,  // d_x_velocity
                     const double &,  // dt
                     const int &,     // num_faces
                     const int &,     // num_extend_cells
                     const double *,  // face_coords
                     const double *,  // cell_extend_coords
                     const int *,     // face_adj_cell_extent
                     const int *,     // face_adj_cell_index
                     const double *,  // uval_sctatch
                     double *);       // flux_new

void conspatch3d_(const double &,   // d_dx
                  const int &,      // num_cells
                  const int &,      // num_faces
                  const double *,   // cell_coords
                  const double *,   // face_coords
                  const int *,      // cell_adj_face_extent
                  const int *,      // cell_adj_face_index
                  double *,         // uval_current
                  double *,         // uval_new
                  const double *);  // flux_new

#endif
}
